!******************************************************************************
! Copyright c 2013, The Ames Laboratory, Iowa State University, and Rutgers
! University*.  All rights reserved.
!
! This software was authored by Yongxin Yao, Nicola Lanata*, Gabriel Kotliar*,
! Cai-Zhuang Wang, and Kai-Ming Ho, at The Ames Laboratory and
! Rutgers University and was supported by the U.S.
! Department of Energy (DOE), Office of Science,
! Basic Energy Sciences, Materials Science and Engineering Division.
! The Ames Laboratory is operated by Iowa State University for DOE
! under U.S. Government contract DE-AC02-07CH11358.
! The U.S. Government has the rights to use, reproduce, and
! distribute this software.
! NEITHER THE GOVERNMENT, THE AMES LABORATORY, IOWA STATE UNIVERSITY,
! NOR RUTGERS UNIVERSITY MAKES ANY WARRANTY,
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
! If software is modified to produce derivative works,
! such modified software should be clearly marked,
! so as to not confuse it with the version available from
! The Ames Laboratory and Rutgers University.
!
! Additionally, redistribution and use in source and binary forms,
! with or without modification,
! are permitted provided that the following conditions are met:
!
!     Redistribution of source code must retain the above copyright notice,
!     this list of conditions, and the following disclaimer.
!
!     Redistribution in binary form must reproduce the above copyright notice,
!     this list of conditions, and the following disclaimer
!     in the documentation and/or other materials provided with distribution.
!
!     Neither the name of The Ames Laboratory, Iowa State University,
!     Rutgers University, the U.S. Government, nor the names of
!     its contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE AMES LABORATORY, IOWA STATE UNIVERSITY,
! RUTGERS UNIVERSITY, AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
! THE IMPLIED WARRANTIES OF MERCHANTABILITY
! AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
! IN NO EVENT SHALL THE GOVERNMENT, THE AMES LABORATORY,
! IOWA STATE UNIVERSITY, RUTGERS UNIVERSITY, OR CONTRIBUTORS BE LIABLE
! FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
! HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
! WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
! OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
! EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!******************************************************************************

module bandstru
    use gprec, only: q
    use gmpi
    use gutil, only: anmxbmm,fermi_fun,gauss_fun,int_to_str,uhau,annxb,hermev
    use gconstant, only: const_z0=>z0,const_z1=>z1,const_iu_kgen=>iu_kgen,pi
    use ghdf5
    use localstore
    use gtime
    implicit none
    private

    type k_points
        integer::dim=-1,diml=-1,ismear,icor,ensemble=0,base=0
        real(q),pointer :: wt(:)=>null()
        real(q) twt,delta,cordeg
        character::file_name*512
    end type k_points

    type sym_info
        integer::nop,mode=0,ie
    end type sym_info

    type,public::bandstru_ob
        !< iso=2 => soc
        integer ispin_in,nspin_in,iso,nspin,rspo,riso,ispin
        integer n_mott,mode_hk ! mode_hk: 0 for original basis, otherwise sab.
        real(q) ehybrd,eband(2),ets2(2),eband_bare(2)
        !< ne(3). 1: total number of bands, 2-3: correlated bands interval
        integer,pointer :: ne(:,:,:)
        real(q),pointer :: ek(:,:,:) ! kohn-sham / gutz eigen values
        real(q),allocatable :: normso(:,:,:) ! list of spin-up/dn weights.
        complex(q),pointer :: r(:,:,:),la1(:,:,:),nc_phy(:,:,:),nks(:,:,:), &
                &d0(:,:,:),nrl(:,:,:) ! for onsite local 1pdm
        complex(q),pointer :: psi0_b(:,:,:,:) ! <bare psi0 | basis orbital>
        complex(q),pointer :: vk (:,:,:,:,:) ! eigen-vector
        complex(q),allocatable :: hk0(:,:,:,:,:) ! kohn-sham bare dispersion
        complex(q),allocatable :: kswt(:,:,:,:) ! ks band renormalized occup.
        integer nmax ! nmax: miximal number of bands
        integer nmaxin !< maximal number of bands inside the enrgy window.
        real(q),pointer :: ferwe(:,:,:)
        !< fermi-weight, convention including kpt%wt and spin degeneracy
        real(q) :: eflda,nelet,nelel,nelec,mup_dn=0._q
        real(q) :: ef=0._q,nelet_mott=0._q
        real(q) ebmin,ebmax,ebwidth
        complex(q),pointer::cwt(:,:,:)   ! occupation matrix for \rho
        type (sym_info) :: sym
        type (k_points) :: kpt

        contains
        procedure::init=>set_bnd_info
        procedure::read_hk0=>read_bare_hamiltonian
        procedure::rotate_hk0=>rotate_bare_hamiltonian
        procedure::corr_ebwidth=>calc_corr_ebwidth
        procedure::calc_all=>calc_band_all
        procedure::rm_h1e_from_hk0=>rm_h1e_from_bare_hamiltonian
        procedure::calc_fermi=>gutz_fermi
        procedure::calc_bz_integration_correction
        procedure::calc_nks
        procedure::write_hk0_sab=>write_revised_bare_hamiltonian
        procedure::modify_mott=>bnd_modify_mott
        procedure::calc_mup_dn
        procedure::calc_da0
        procedure::calc_rnrl
        procedure::calc_kswt
        procedure::write_kswt=>h5write_kswt
        procedure::write_results=>h5write_bnd_results
        procedure::write_bands=>h5write_bands
        procedure::calc_ebnd
        procedure::calc_num_electrons
    end type bandstru_ob
      
    contains


    subroutine set_bnd_info(this,gh5,mpi,loc)
    class(bandstru_ob) :: this
    class(hdf5_ob) :: gh5
    class(mpi_ob) :: mpi
    class(localstore_ob) :: loc

    integer iso,isp,err
    logical lexist

    this%nelet_mott=sum(loc%mott(:)%nelect)
    this%n_mott=sum(loc%mott(:)%nsorb)

    call gh5%fopen('GBareH.h5',1,"r")
    call gh5%gopen("/",1,1)
    call gh5%read(iso,'iso',1)
    call gh5%read(this%ispin_in,'ispin',1)
    call gh5%read(this%kpt%dim,'kptdim',1)
    call gh5%read(this%nmax,'nbmax',1)
    loc%ispin_in=this%ispin_in
    this%ispin=loc%ispin
    if(mpi%io>0)then
        write(mpi%io,'(" ispin_in=",i2," iso=",i2, &
                &" ispin=",i2," ispo=",i2)')this%ispin_in, &
                &loc%iso,loc%ispin,loc%ispo
    endif
    if(iso==2.and.loc%ispin==2)then
        allocate(this%normso(this%nmax,2,this%kpt%dim))
        call gh5%read(this%normso,this%nmax,2,this%kpt%dim,'/normso',1)
    endif
    ! modify ispin_in for iso=2
    this%nspin_in=max(1, this%ispin_in/iso)
    allocate(this%ne(3,this%kpt%dim,this%nspin_in))
    call gh5%exists(1,'/ispin_0/ne_list',lexist)
    if(lexist)then
        ! upper-case name: fortran convention
        do isp=1,this%nspin_in
            call gh5%read(this%ne(:,:,isp),3,this%kpt%dim,'/ispin_'//&
                    &trim(int_to_str(isp-1))//"/ne_list",1)
        enddo
        ! 0-base to 1-base
        this%ne(2:3,:,:)=this%ne(2:3,:,:)+1
    else
        this%ne(1,:,:)=this%nmax
        this%ne(2,:,:)=1
        this%ne(3,:,:)=this%nmax
    endif

    call gh5%exists(1,'ensemble',lexist)
    if(lexist)then
        call gh5%read(this%kpt%ensemble,'ensemble',1)
    endif

    call gh5%exists(1,'chempot',lexist)
    if(lexist)then
        call gh5%read(this%ef,'chempot',1)
    endif

    ! check shft_init_la1
    call gh5%exists(1,'shft_init_la1',lexist)
    if(lexist)then
        call gh5%read(loc%shft_init_la1,'shft_init_la1',1)
    endif

    allocate(this%kpt%wt(this%kpt%dim))
    call gh5%read(this%kpt%wt(1),'kptwt',1)
    call gh5%read(this%kpt%ismear,'ismear',1)
    call gh5%read(this%kpt%delta,'delta',1)
    call gh5%read(this%nelet,'nelectron',1)
    call gh5%exists(1,'kptfname',lexist)
    if(lexist)then
        call gh5%read(this%kpt%file_name,512,'kptfname',1)
    endif
    call gh5%read(this%sym%nop,'symnop',1)
    call gh5%read(this%sym%ie,'SYMIE',1)
    call gh5%gclose(1)
    call loc%read_h1e(gh5,this%nspin_in,1)
    call gh5%fclose(1)
    call loc%rotate_h1e()
    call loc%herm_matrices_pp(loc%h1e,'h1e',loc%hm,.false.,mpi%io,-1)
    this%kpt%twt=sum(this%kpt%wt)
    this%kpt%wt=this%kpt%wt/this%kpt%twt
    ! Consistence check
    if(loc%iso /= iso)then
        stop ' error in set_bnd_info: loc%iso /= this%iso!'
    endif
    this%iso=iso
    this%nspin=loc%nspin
    this%rspo=loc%rspo
    this%riso=loc%riso

    this%nmaxin=maxval(this%ne(3,:,:)-this%ne(2,:,:)+1)
    if(mpi%io>0)then
        write(mpi%io,'(" min/max(this%ne(1,:,:))=",2i8)') &
                &minval(this%ne(1,:,:)), &
                &maxval(this%ne(1,:,:))
        write(mpi%io,'(" min/max(this%ne(2,:,:))=",2i8)') &
                &minval(this%ne(2,:,:)), &
                &maxval(this%ne(2,:,:))
        write(mpi%io,'(" min/max(this%ne(3,:,:))=",2i8)') &
                &minval(this%ne(3,:,:)), &
                &maxval(this%ne(3,:,:))
    endif

    this%nelec=0
    do isp=1,this%nspin_in
        this%nelec=this%nelec+sum((this%ne(2,:,isp)-1)*this%kpt%wt)*this%rspo
    enddo
    this%nelel=this%nelet-this%nelec
    if(mpi%io>0)then
        write(mpi%io,'(" this%nmax=",i4)')this%nmax
        write(mpi%io,'(" valence electrons: total=",f8.1)')this%nelet
        write(mpi%io,'("         correlated block=",f8.1)')this%nelel
    endif
    call set_kpt_icor(this,mpi%io)
    call set_kpt_diml(this,mpi)
    call alloc_bnd(this,loc)
    return
      
    end subroutine set_bnd_info
    

    subroutine set_kpt_diml(this,mpi)
    class(bandstru_ob) :: this
    class(mpi_ob) :: mpi

    integer nk_residual

    associate(nk1=>this%kpt%diml, kbase=>this%kpt%base)
        nk1=this%kpt%dim/mpi%nprocs
        nk_residual=this%kpt%dim-nk1*mpi%nprocs
        if(mpi%myrank<nk_residual)then
            nk1=nk1+1
        endif

        kbase=nk1*mpi%myrank
        if(nk_residual>0.and.mpi%myrank>=nk_residual)then
            kbase=kbase+min(mpi%myrank,nk_residual)
        endif
        nk1=min(nk1,this%kpt%dim-kbase)
    end associate
    return
    
    end subroutine set_kpt_diml


    subroutine alloc_bnd(this,loc)
    class(bandstru_ob) :: this
    class(localstore_ob) :: loc

    integer i
    
    allocate(this%ferwe(this%nmax,this%kpt%dim,this%nspin))
    allocate(this%ek(this%nmax,this%kpt%dim,this%nspin))
    this%ek=0
    allocate(this%hk0(this%nmaxin,this%nmaxin,this%sym%nop, &
            &this%kpt%diml,this%nspin_in))
    allocate(this%vk (this%nmaxin,this%nmaxin,this%sym%nop, &
            &this%kpt%diml,this%nspin))
    allocate(this%psi0_b(this%nmaxin,this%nmaxin,this%kpt%diml, &
            &this%nspin_in))
    this%psi0_b=0
    allocate(this%r(loc%nasotot,loc%nasotot,this%nspin))
    this%r=0
    forall(i=1:loc%nasotot) &
            &this%r(i,i,:)=1._q
    allocate(this%d0 (loc%nasotot,loc%nasotot,this%nspin))
    allocate(this%la1(loc%nasotot,loc%nasotot,this%nspin))
    this%la1=0
    allocate(this%nc_phy(loc%nasotot,loc%nasotot,this%nspin))
    this%nc_phy=0
    allocate(this%nks(loc%nasotot,loc%nasotot,this%nspin))
    this%nks=0
    allocate(this%nrl(loc%nasotot,loc%nasotot,this%nspin))
    this%nrl=0
    return
      
    end subroutine alloc_bnd
    

    subroutine read_bare_hamiltonian(this,gh5)
    class(bandstru_ob) :: this
    class(hdf5_ob) :: gh5

    integer isp
      
    call gh5%fopen("GBareH.h5",1,"r")
    do isp=1,this%nspin_in
        call gh5%gopen("/",1,1)
        call gh5%read(this%mode_hk,'mode_hk',1)
        call gh5%gclose(1)
        call gh5%read(this%ek(:,:,isp),this%nmax,this%kpt%dim,"/ispin_"// &
                &trim(int_to_str(isp-1))//"/e_list",1)
        call gh5%dopen("/ispin_"//trim(int_to_str(isp-1))// &
                &"/U_PSIK0_TO_HK0_BASIS_LIST",1,1)
        call gh5%read(this%psi0_b(:,:,:,isp),this%nmaxin,this%nmaxin, &
                &this%kpt%diml,0,0,this%kpt%base,1)
        call gh5%dclose(1)
        call gh5%dopen("/ispin_"//trim(int_to_str(isp-1))//"/HK0_LIST",1,1)
        call gh5%read(this%hk0(:,:,:,:,isp),this%nmaxin,this%nmaxin, &
                &this%sym%nop,this%kpt%diml,0,0,0,this%kpt%base,1)
        call gh5%dclose(1)
    enddo
    call gh5%fclose(1)

    if(this%iso==1.and.this%nspin_in==1.and.this%nspin==2)then
        this%ek(:,:,2)=this%ek(:,:,1)
    endif
    return
      
    end subroutine read_bare_hamiltonian
   

    subroutine write_revised_bare_hamiltonian(this,gh5)
    class(bandstru_ob) :: this
    class(hdf5_ob) :: gh5

    integer isp

    call gh5%fopen("GBareH.h5",1,"rw")
    do isp=1,this%nspin_in
        call gh5%gopen("/",1,1)
        call gh5%awrite(1,'mode_hk',1)
        call gh5%gclose(1)
        call gh5%dopen("/ispin_"//trim(int_to_str(isp-1))// &
                &"/U_PSIK0_TO_HK0_BASIS_LIST",1,1)
        call gh5%dwrite(this%psi0_b(:,:,:,isp),this%nmaxin,this%nmaxin,&
                &this%kpt%diml,0,0,this%kpt%base,1)
        call gh5%dclose(1)
        call gh5%dopen("/ispin_"//trim(int_to_str(isp-1))//"/HK0_LIST",1,1)
        call gh5%dwrite(this%hk0(:,:,:,:,isp),this%nmaxin,this%nmaxin, &
                &this%sym%nop,this%kpt%diml,0,0,0,this%kpt%base,1)
        call gh5%dclose(1)
    enddo
    call gh5%fclose(1)
    return

    end subroutine write_revised_bare_hamiltonian


    ! rotate bare hamiltonian to symmetry-adapted basis.
    subroutine rotate_bare_hamiltonian(this,loc)
    class(bandstru_ob) :: this
    class(localstore_ob) :: loc

    integer ikp,ikpl,isym,nbands,isp
    complex(q),target::zbuf(this%nmaxin**2)
    complex(q),pointer::p_2d(:,:)

    ! check if necessary
    if(.not.associated(loc%db2sab).or.this%mode_hk/=0)return

    do ikpl=1,this%kpt%diml
        ikp=this%kpt%base+ikpl
        do isp=1,this%nspin_in
            nbands=this%ne(3,ikp,isp)-this%ne(2,ikp,isp)+1
            ! contiguous block
            p_2d(1:nbands,1:nbands)=>zbuf(1:nbands**2)
            p_2d=this%psi0_b(1:nbands,1:nbands,ikpl,isp)
            call loc%rtrans(p_2d,nbands,1)
            this%psi0_b(1:nbands,1:nbands,ikpl,isp)=p_2d
            do isym=1,this%sym%nop
                p_2d=this%hk0(1:nbands,1:nbands,isym,ikpl,isp)
                call loc%rtrans(p_2d,nbands,1)
                p_2d=transpose(p_2d)
                call loc%rtrans(p_2d,nbands,-1)
                this%hk0(1:nbands,1:nbands,isym,ikpl,isp)=transpose(p_2d)
            enddo
        enddo
    enddo
    nullify(p_2d)
    return
      
    end subroutine rotate_bare_hamiltonian


    subroutine merge_h1e(this, loc, h1e)
    class(bandstru_ob) :: this
    class(localstore_ob) :: loc
    complex(q),intent(out) :: h1e(loc%nasotot,loc%nasotot,this%nspin_in)

    integer i,nbase,nstep,naso,isp

    h1e=0
    nbase=1
    nstep=2/loc%iso
    do i=1,loc%num_imp
        naso=loc%co(i)%dimso
        do isp=1,this%nspin_in
            h1e(nbase:nbase+naso-1,nbase:nbase+naso-1,isp)= &
                    &loc%co(i)%h1e(isp::nstep,isp::nstep)
        enddo
        nbase=nbase+naso
    enddo
    return

    end subroutine merge_h1e


    subroutine rm_h1e_from_bare_hamiltonian(this,loc,mode)
    class(bandstru_ob) :: this
    class(localstore_ob) :: loc
    integer,intent(in) :: mode

    integer isym,ikp,ikpl,isp,nbands
    complex(q) h1e(loc%nasotot,loc%nasotot,this%nspin_in)

    call merge_h1e(this, loc, h1e)
    associate(nasot=>loc%nasotot)
    do ikpl=1,this%kpt%diml
        ikp=this%kpt%base+ikpl
        do isp=1,this%nspin_in
            nbands=this%ne(3,ikp,isp)-this%ne(2,ikp,isp)+1
            do isym=1,this%sym%nop
                this%hk0(1:nasot,1:nasot,isym,ikpl,isp)=&
                        &this%hk0(1:nasot,1:nasot,isym,ikpl,isp)+mode*&
                        &h1e(:,:,isp)
            enddo
        enddo
    enddo
    end associate
    return

    end subroutine rm_h1e_from_bare_hamiltonian

    !*************************************************************************  
    !< calc band energy
    !*************************************************************************
    subroutine calc_ebnd(this,mpi,bare_mode)
    class(bandstru_ob) :: this
    class(mpi_ob) :: mpi
    integer,optional,intent(in) :: bare_mode

    integer ikp,ikpl,isp,ispp,nemin,nemax,nbands,i,j
    complex(q),pointer::p_a(:,:),p_v(:,:),p_psi0b(:,:)
    complex(q),target::buf_a(this%nmaxin**2),buf_b(this%nmaxin**2), &
            &buf_v(this%nmaxin**2)
     
    this%eband=0
    if(this%ispin_in/=2.or.this%iso/=2)then
        this%eband=this%ets2
        do ikp=1,this%kpt%dim
            do isp=1,this%nspin
                ispp=min(isp,this%nspin_in)
                this%eband(isp)=this%eband(isp)+ &
                        &sum(this%ek(:,ikp,isp)*this%ferwe(:,ikp,isp))
            enddo
        enddo
    else
        ! this is the case of magnetic calculations with spin-orbit.
        ! we try to get spin-resolved band energy.
        ikpl=0
        do ikpl=1,this%kpt%diml
            ikp=this%kpt%base+ikpl
            nemin=this%ne(2,ikp,1)
            nemax=this%ne(3,ikp,1)
            nbands=nemax-nemin+1
            p_a(1:nbands,1:nbands)=>buf_a(1:nbands**2)
            p_psi0b(1:nbands,1:nbands)=>buf_b(1:nbands**2)
            p_v(1:nbands,1:nbands)=>buf_v(1:nbands**2)

            p_psi0b=this%psi0_b(1:nbands,1:nbands,ikpl,1)
            p_v=this%vk(1:nbands,1:nbands,this%sym%ie,ikpl,1)
            call zgemm('n','n',nbands,nbands,nbands,const_z1,p_psi0b, &
                    &nbands,p_v,nbands,const_z0,p_a,nbands)

            do isp=1,2
                if(nemin>1)then
                    this%eband(isp)=this%eband(isp)+sum(&
                            &this%ek(:nemin-1,ikp,1)* &
                            &this%ferwe(:nemin-1,ikp,1)* &
                            &this%normso(:nemin-1,isp,ikp))
                endif
                do i=1,nbands; do j=1,nbands
                    this%eband(isp)=this%eband(isp)+ &
                            &this%ek(nemin-1+i,ikp,1)*p_a(i,j)* &
                            &conjg(p_a(i,j))*this%ferwe(nemin-1+i,ikp,1)* &
                            &this%normso(nemin-1+j,isp,ikp)
                enddo; enddo
            enddo
        enddo
        call mpi%sum_all(this%eband,2)
        this%eband=this%eband+this%ets2/2
        nullify(p_a,p_v,p_psi0b)
    endif
    if(present(bare_mode))then
        this%eband_bare = this%eband
    endif
    return
      
    end subroutine calc_ebnd


    subroutine calc_e_hybrd(this,loc)
    class(bandstru_ob) :: this
    class(localstore_ob) :: loc

    this%ehybrd=real(sum(loc%r*(loc%d0)),q)
    return

    end subroutine calc_e_hybrd


    !*************************************************************************
    !< check magnetic moment for the case of spin-polarized calculations.
    !*************************************************************************
    subroutine chk_mag_moment(this,io)
    class(bandstru_ob) :: this
    integer,intent(in)::io

    integer isp,ikp
    real(q) nel(2)

    if(io<0)return
    if(this%nspin==1.or.this%iso==2)return
    nel=0
    do isp=1,this%nspin
        do ikp=1,this%kpt%dim
            nel(isp) = nel(isp) + sum(this%ferwe(:,ikp,isp))
        enddo
    enddo
    write(io,'(" total (quasi-particle) magnetic moment = ", f12.4)') &
            &nel(2)- nel(1)
    return

    end subroutine chk_mag_moment


    !*************************************************************************
    !> f_n <psi_n|a> r_{a,A} r^+_{B,b} <b|psi_n>
    !! =r^+_{B,b} <b|psi_n> f_n <psi_n|a> r_{a,A}
    !*************************************************************************
    subroutine calc_nabr_1k(this,loc,nabr,vk,ferwe,nbands,wtk,wtk0,isp)
    class(bandstru_ob) :: this
    class(localstore_ob) :: loc
    integer,intent(in)       :: nbands,isp
    real(q),intent(in)       :: wtk,wtk0,ferwe(nbands)
    complex(q),intent(in)    :: vk(nbands,nbands)
    complex(q),intent(inout) :: nabr(nbands,nbands)

    integer ia
    complex(q) r(loc%nasotot,loc%nasotot)

    !< r^+_{B,b}<b|psi>f<psi|a>
    call calc_nabr1_1k(this,loc,nabr,vk,ferwe,nbands,wtk,isp)
    r=this%r(:,:,isp)
    associate(nasot=>loc%nasotot)
        !< r^+_{B,b}<b|psi>f<psi|a>r(a,A)
        call anmxbmm('n',nabr(:,1:nasot),r,nbands,nasot)
        nabr=transpose(nabr) ! to {A,B}
        nabr(1:nasot,1:nasot)=nabr(1:nasot,1:nasot)+ &
                &(-this%nrl(:,:,isp)+this%nc_phy(:,:,isp))*wtk0
    end associate
    return

    end subroutine calc_nabr_1k


    !*************************************************************************
    !> r^+ applied to right side only, for d
    !! r^+ <b|psi>f<psi|a>. 
    !*************************************************************************
    subroutine calc_nabr1_1k(this,loc,nabr,vk,ferwe,nbands,wtk,isp)
    class(bandstru_ob) :: this
    class(localstore_ob) :: loc
    integer,intent(in)       :: nbands,isp
    real(q),intent(in)       :: wtk,ferwe(nbands)
    complex(q),intent(in)    :: vk(nbands,nbands)
    complex(q),intent(out) :: nabr(nbands,nbands)

    complex(q) zr(loc%nasotot,loc%nasotot),nabr_sub(loc%nasotot,nbands)

    call calc_nab_1k(nabr,vk,ferwe,nbands,nbands,wtk)
    zr=this%r(:,:,isp)
    nabr_sub=nabr(1:loc%nasotot,:)
    ! r^+ <b|psi>f<psi|a>
    call annxb('c',zr,nabr_sub,loc%nasotot,nbands,nabr(1:loc%nasotot,:))
    return

    end subroutine calc_nabr1_1k
     

    !*************************************************************************
    !D0^{t}_{A,a} 
    != f_n <psi_n|a> h_{A,B} r^+_{B,b} <b|psi_n>
    != h_{A,B} r^+_{B,b} <b|psi_n> f_n <psi_n|a>
    !*************************************************************************
    subroutine add_da_1k(this,loc,hk0,vk,ferwe,nbands,wtk,isp)
    class(bandstru_ob) :: this
    class(localstore_ob) :: loc
    integer,intent(in)    :: nbands,isp
    real(q),intent(in)    :: wtk,ferwe(nbands)
    complex(q),intent(in) :: hk0(nbands,nbands),vk(nbands,nbands)

    complex(q) nabr(nbands,nbands),zd(loc%nasotot,loc%nasotot), &
            &hk0_(loc%nasotot,nbands)

    ! r^+ <b|psi>f<psi|a>
    call calc_nabr1_1k(this,loc,nabr,vk,ferwe,nbands,wtk,isp)
    hk0_=hk0(1:loc%nasotot,:)
    call zgemm('n','n',loc%nasotot,loc%nasotot,nbands,const_z1, &
            &hk0_,loc%nasotot,nabr(:,1:loc%nasotot),nbands, &
            &const_z0,zd,loc%nasotot)
    ! D0 can be real or complex type .
    this%d0(:,:,isp)=this%d0(:,:,isp)+zd
    return

    end subroutine add_da_1k


    !< calculate n_{ab} at a single k-point. 
    subroutine calc_nab_1k(nab,vk,ferwe,nbasis,nbands,wtk)
    integer,intent(in)     :: nbasis,nbands
    real(q),intent(in)     :: wtk,ferwe(nbands)
    complex(q),intent(in)  :: vk(nbasis,nbands)
    complex(q),intent(out) :: nab(nbasis,nbasis)

    integer ib
    complex(q) vf(nbasis,nbands)

    do ib=1,nbands
        vf(:,ib)=vk(:,ib)*ferwe(ib)*wtk
    enddo
    call zgemm('n','c',nbasis,nbasis,nbands,const_z1,vf,nbasis,vk, &
            &nbasis,const_z0,nab,nbasis) ! <b|psi>f<psi|a>
    return

    end subroutine calc_nab_1k


    subroutine calc_da0(this,loc,mpi)
    class(bandstru_ob) :: this
    class(localstore_ob) :: loc
    class(mpi_ob) :: mpi

    integer isym,ikp,ikpl,isp,ispp,nbands,nemin,nemax
    real(q) wtk
    complex(q),target :: buf_vk(this%nmaxin**2),buf_hk0(this%nmaxin**2)
    complex(q),pointer :: p_vk(:,:),p_hk0(:,:)
    real(q),pointer :: ferwe(:)

    this%d0=0
    wtk=1._q/this%sym%nop/this%rspo
    do ikpl=1,this%kpt%diml
        ikp=this%kpt%base+ikpl
        do isp=1,this%nspin
            ispp=min(isp,this%nspin_in)
            nemin=this%ne(2,ikp,ispp)
            nemax=this%ne(3,ikp,ispp)
            nbands=nemax-nemin+1
            ferwe=>this%ferwe(nemin:nemax,ikp,isp)
            p_hk0(1:nbands,1:nbands)=>buf_hk0(1:nbands**2)
            p_vk(1:nbands,1:nbands)=>buf_vk(1:nbands**2)
            do isym=1,this%sym%nop
                p_hk0=this%hk0(1:nbands,1:nbands,isym,ikpl,ispp)
                p_vk=this%vk(1:nbands,1:nbands,isym,ikpl,isp)
                call add_da_1k(this,loc,p_hk0,p_vk,ferwe,nbands,wtk,isp)
            enddo                    
        enddo
    enddo
    call mpi%sum_all(this%d0(:,1,1),loc%nasotot*loc%nasotot*this%nspin)
    nullify(p_vk,p_hk0,ferwe)

    ! back to D0_{a,A}
    do isp=1,this%nspin
        this%d0(:,:,isp)=transpose(this%d0(:,:,isp))
    enddo
    call loc%map_bnd_matrix(loc%d0,this%d0,.true.)
    return

    end subroutine calc_da0


    subroutine calc_nks(this,loc,mpi,ferwes)
    class(bandstru_ob) :: this
    class(localstore_ob) :: loc
    class(mpi_ob) :: mpi
    real(q),optional,target::ferwes(this%nmax,this%kpt%dim,this%nspin)

    integer ikpl,ikp,isym,isp,ispp,nbands,nbase,&
            &n1,n2,i,naso,nemin,nemax
    real(q) wtk
    real(q),pointer :: p_ferwe(:),ferwes_(:,:,:)
    complex(q),pointer :: p_vk(:,:),p_vk_sub(:,:),p_nks(:,:)
    complex(q),target::nks(loc%nasomax*loc%nasomax)
    complex(q),target::vk_sub(loc%nasomax*this%nmaxin),buf_vk(this%nmaxin**2)

    if(present(ferwes))then
        ferwes_=>ferwes
    else
        ferwes_=>this%ferwe
    endif
    loc%nks=0

    wtk=1._q/this%sym%nop/this%rspo
    do ikpl=1,this%kpt%diml
        ikp=this%kpt%base+ikpl
        do isp=1,this%nspin
            ispp=min(isp,this%nspin_in)
            nemin=this%ne(2,ikp,ispp)
            nemax=this%ne(3,ikp,ispp)
            nbands=nemax-nemin+1
            p_ferwe=>ferwes_(nemin:nemax,ikp,isp)
            p_vk(1:nbands,1:nbands)=>buf_vk(1:nbands**2)
            do isym=1,this%sym%nop
                p_vk=this%vk(1:nbands,1:nbands,isym,ikpl,isp)
                nbase=0
                do i=1,loc%num_imp
                    naso=loc%co(i)%dimso
                    n1=1+(isp-1)*naso; n2=naso*isp
                    p_vk_sub(1:naso,1:nbands)=>vk_sub(1:naso*nbands)
                    p_vk_sub=p_vk(nbase+1:nbase+naso,1:nbands)
                    p_nks(1:naso,1:naso)=>nks(1:naso*naso)
                    call calc_nab_1k(p_nks,p_vk_sub,p_ferwe,naso, &
                            &nbands,wtk)
                    loc%co(i)%nks(n1:n2,n1:n2)=loc%co(i)%nks(n1:n2,n1:n2) +&
                            &transpose(p_nks)
                    nbase=nbase+naso
                enddo
            enddo
        enddo
    enddo
    nullify(p_ferwe,p_vk,p_vk_sub,p_nks)
    call mpi%sum_all(loc%nks,loc%na2112)
    do i=1,loc%num_imp
        call loc%co(i)%nks_patch_order(loc%ispo,loc%iso)
    enddo
    return

    end subroutine calc_nks


    !*************************************************************************
    !> f_n <psi_n|a> r_{a,A} r^+_{B,b} <b|psi_n>
    !! =r^+_{B,b} <b|psi_n> f_n <psi_n|a> r_{a,A}
    !*************************************************************************
    subroutine calc_rnrl(this,loc)
    class(bandstru_ob) :: this
    class(localstore_ob) :: loc

    integer isp
      
    associate(n=>loc%nasotot)
    do isp=1,this%nspin
        !< <b|psi_n> f_n <psi_n|a> r_{a,A}
        call zgemm('t','n',n,n,n,const_z1,this%nks(:,:,isp),n, &
                &this%r(:,:,isp),n,const_z0,this%nrl(:,:,isp),n)
        !< r^+_{B,b}<b|psi_n> f_n <psi_n|a> r_{a,A}
        call annxb('c',this%r(:,:,isp),this%nrl(:,:,isp),n,n)
        !< bring back to {A,B}
        this%nrl(:,:,isp)=transpose(this%nrl(:,:,isp))
    enddo
    end associate
    return
      
    end subroutine calc_rnrl
      

    !*************************************************************************
    subroutine calc_band_all(this,loc,time,mpi)
    class(bandstru_ob)::this
    class(localstore_ob)::loc
    class(time_ob)::time
    class(mpi_ob)::mpi

    integer isym,ikp,ikpl,nbands,isp,ispp
    complex(q),pointer :: hk(:,:)
    real(q),allocatable :: ek_list(:,:,:)
    complex(q),target :: buf_hk(this%nmaxin**2)
    complex(q) :: r(loc%nasotot,loc%nasotot),la1(loc%nasotot,loc%nasotot)
     
    call time%start(2)

    allocate(ek_list(this%nmax,this%kpt%dim,this%nspin))
    ek_list=0

    !$omp parallel do &
    !$omp &private(isym,isp,ispp,ikpl,ikp,nbands,hk,r,la1,buf_hk) &
    !$omp &schedule(static,1)
        do ikpl=1,this%kpt%diml
            ikp=this%kpt%base+ikpl
            do isym=1,this%sym%nop
                do isp=1,this%nspin
                    ispp=min(isp,this%nspin_in)
                    nbands=this%ne(3,ikp,ispp)-this%ne(2,ikp,ispp)+1
                    hk(1:nbands,1:nbands)=>buf_hk(1:nbands**2)
                    r=this%r(:,:,isp)
                    la1=this%la1(:,:,isp)
                    hk=this%hk0(1:nbands,1:nbands,isym,ikpl,ispp)
                    if(isym==1)then
                        ! save the energies outside of the energy window.
                        ek_list(:,ikp,isp)=this%ek(:,ikp,isp)
                        call calc_band_1k(hk,r,la1,nbands,loc%nasotot, &
                                &ek_list(this%ne(2,ikp,ispp): &
                                &this%ne(3,ikp,ispp),ikp,isp))
                    else
                        call calc_band_1k(hk,r,la1,nbands,loc%nasotot)   
                    endif
                    this%vk(1:nbands,1:nbands,isym,ikpl,isp)=hk
                enddo
            enddo !isym
        enddo ! ikpl
    !$omp end parallel do
    nullify(hk)
     
    call mpi%sum_all(ek_list(:,1,1),this%nmax*this%kpt%dim*this%nspin)
    this%ek=ek_list
    deallocate(ek_list)

    call time%finish(2)
    call time%print_usage('calc_band_all',2,mpi%io)
    return 

    end subroutine calc_band_all


    !*************************************************************************
    subroutine calc_band_1k(hk,r,la1,nbands,nasot,ek)
    integer,intent(in)::nbands,nasot
    complex(q),intent(in)::r(nasot,nasot),la1(nasot,nasot)
    complex(q),intent(inout)::hk(nbands,nbands)
    real(q),optional,intent(out)::ek(nbands)

    real(q) w(nbands)

    call calc_hamil_1k(hk,r,la1,nbands,nasot)
    call hermev('v','u',hk,w,nbands)
    if(present(ek))then
        ek=w
    endif
    return

    end subroutine calc_band_1k


    subroutine calc_kswt(this,loc,mpi,time)
    class(bandstru_ob) :: this
    class(localstore_ob) :: loc
    class(mpi_ob)::mpi
    class(time_ob)::time

    integer ikp,ikpl,nbands,isp,i,ispp,nsp,naso,nbase
    real(q) sumwt(2),maxoffdiag,fac_irspo
    complex(q),pointer :: p_vk(:,:),p_kswt(:,:),p_uk(:,:),p_utrans(:,:)
    complex(q),target::utrans(loc%nasomax**2),buf_kswt(this%nmaxin**2),&
            &buf_uk(this%nmaxin**2),buf_vk(this%nmaxin**2)
    real(q),pointer :: p_ferwe(:)
    
    allocate(this%kswt(this%nmaxin,this%nmaxin,this%kpt%diml,this%nspin_in))
    this%kswt=0
    call time%start(2)
    sumwt=0
    maxoffdiag=0
    fac_irspo=1._q/this%rspo

    !$omp parallel do &
    !$omp &private(isp,ispp,ikpl,ikp,i,nbands,nbase,nsp, &
    !$omp         &p_kswt,p_ferwe,p_vk,p_uk,utrans,p_utrans, &
    !$omp         &naso,buf_kswt,buf_uk,buf_vk) &
    !$omp &schedule(static,1) &
    !$omp &reduction(+:sumwt) &
    !$omp &reduction(max:maxoffdiag)
    do ikpl=1,this%kpt%diml
        ikp=this%kpt%base+ikpl
        do ispp=1,this%nspin_in
            if(this%nspin_in==1)then
                nsp=this%nspin
            else
                nsp=ispp
            endif
            nbands=this%ne(3,ikp,ispp)-this%ne(2,ikp,ispp)+1
            p_kswt(1:nbands,1:nbands)=>buf_kswt(1:nbands**2)
            p_uk(1:nbands,1:nbands)=>buf_uk(1:nbands**2)
            p_vk(1:nbands,1:nbands)=>buf_vk(1:nbands**2)
            p_kswt=0
            p_uk=this%psi0_b(1:nbands,1:nbands,ikpl,ispp)
            do isp=ispp,nsp
                p_ferwe=>this%ferwe(this%ne(2,ikp,ispp): &
                        &this%ne(3,ikp,ispp),ikp,isp)
                p_vk=this%vk(1:nbands,1:nbands,this%sym%ie,ikpl,isp)
                call calc_kswt_1k(this,loc,p_kswt,p_vk,p_uk,p_ferwe,nbands, &
                        &fac_irspo,this%kpt%wt(ikp),isp)
                sumwt(1)=sumwt(1)+sum(p_ferwe)
            enddo ! isp
            p_kswt=p_kswt*this%rspo
            this%kswt(1:nbands,1:nbands,ikpl,ispp)=p_kswt

            do i=1,nbands
                sumwt(2)=sumwt(2)+real(p_kswt(i,i),q)
                if(i==nbands)exit
                maxoffdiag=max(maxoffdiag, &
                        &maxval(abs(p_kswt(i+1:,i)))/this%kpt%wt(ikp))
            enddo
        enddo
    enddo ! ikpl
    !$omp end parallel do
    nullify(p_vk,p_uk,p_ferwe,p_kswt,p_utrans)
      
    call mpi%barrier()
    call mpi%sum_all(sumwt,2)
    call mpi%max_all(maxoffdiag)

    if(mpi%io>0)then
        write(mpi%io,'(" correlated subspace, sum_ferwt = ",f0.7, &
                &" sum_kswt = ",f0.7)')sumwt
        write(mpi%io,'("     max off diagonal of psik occ. mat. = ", &
                &f0.6)')maxoffdiag
        sumwt(1) = sumwt(1)-sumwt(2)
        if(abs(sumwt(1))>1.e-4_q)then
            write(0,'(" warning: too large sum_ferwt-sum_kswt = ", &
                    &f0.5,"!")')sumwt(1)
            write(mpi%io,'(" warning: too large sum_ferwt-sum_kswt = ", &
                    &f0.5,"!")')sumwt(1)
        endif
    endif

    call time%finish(2)
    call time%print_usage('calc_kswt',2,mpi%io)
    return 

    end subroutine calc_kswt


    subroutine calc_kswt_1k(this,loc,kswt,vk,uk,ferwe,nbands,wtk,wtk0,isp)
    class(bandstru_ob) :: this
    class(localstore_ob) :: loc
    integer,intent(in)::nbands,isp
    real(q),intent(in)::ferwe(nbands),wtk,wtk0
    complex(q),intent(in)::vk(nbands,nbands),uk(nbands,nbands)
    complex(q),intent(inout)::kswt(nbands,nbands)

    complex(q) nabr(nbands,nbands)

    ! ~ f<psi|a><b|psi>
    call calc_nabr_1k(this,loc,nabr,vk,ferwe,nbands,wtk,wtk0,isp)
    nabr=transpose(nabr) ! ~ form <a|psi>f<psi|b>
    call uhau(nabr,uk,nbands,nbands,trul='n',trur='c')
    kswt=kswt+nabr ! density matrix <a|psi>f_psi<psi|b>
    return

    end subroutine calc_kswt_1k


    subroutine h5write_kswt(this,loc,gh5)
    class(bandstru_ob) :: this
    class(localstore_ob) :: loc
    class(hdf5_ob) :: gh5

    real(8) e_gamma_dc,eband
    integer ikpl,isp

    call gh5%fopen('KSWT.h5',1,"w")
    call gh5%gopen("/",1,1)
    call gh5%awrite(this%ef,'efermi',1)
    call gh5%awrite(loc%egamma_dc,'egamma_dc',1)
    eband=sum(this%eband)
    call gh5%awrite(eband,'eband',1)
    if(this%ispin_in==2)then
        call gh5%awrite(this%eband(1),'eband_spin1',1)
        call gh5%awrite(this%eband(2),'eband_spin2',1)
    endif
    call gh5%gclose(1)
    call gh5%dwrite(this%kpt%wt,this%kpt%dim,'/kpt_wt',1)
    call gh5%dwrite(this%ne,3,this%kpt%dim,this%nspin_in,'/bnd_ne',1)
    call gh5%dcreate(4,(/this%nmaxin,this%nmaxin,this%kpt%dim,this%nspin_in/), &
            & "/KSWT",gh5%dcomplex_id,1,1)
    call gh5%dwrite(this%kswt,this%nmaxin,this%nmaxin,this%kpt%diml, &
            &this%nspin_in,0,0,this%kpt%base,0,1,serialio=.true.)
    call gh5%dclose(1)
    call gh5%fclose(1)
    return

    end subroutine h5write_kswt


    subroutine h5write_bnd_results(this,loc,gh5,i_g)
    class(bandstru_ob)::this
    class(localstore_ob)::loc
    class(hdf5_ob)::gh5
    integer,intent(in)::i_g

    complex(q) h1e(loc%nasotot,loc%nasotot,this%nspin_in)

    call gh5%awrite(this%ef,"efermi",i_g)
    call gh5%awrite(this%eband,2,"ebands",i_g)
    call gh5%awrite(this%eband_bare,2,"ebands_bare",i_g)
    call gh5%awrite(this%ets2,2,"ets2",i_g)
    call gh5%awrite(this%nelet,"nelectrons",i_g)
    call gh5%awrite(this%r, &
            &loc%nasotot,loc%nasotot,this%nspin,"RMAT",i_g)
    call gh5%awrite(this%la1, &
            &loc%nasotot,loc%nasotot,this%nspin,"LAMAT",i_g)
    call gh5%awrite(this%nrl, &
            &loc%nasotot,loc%nasotot,this%nspin,"NRLMAT",i_g)
    call gh5%awrite(this%nc_phy, &
            &loc%nasotot,loc%nasotot,this%nspin,"NPHYMAT",i_g)
    call merge_h1e(this,loc,h1e)
    call gh5%awrite(h1e, &
            &loc%nasotot,loc%nasotot,this%nspin_in,"H1E",i_g)
    return

    end subroutine h5write_bnd_results


    subroutine h5write_bands(this,loc,gh5)
    class(bandstru_ob)::this
    class(hdf5_ob)::gh5
    class(localstore_ob) :: loc

    call gh5%fopen("GBands.h5",1,"w")
    call gh5%gopen("/",1,1)
    call gh5%awrite(this%ef,"efermi",1)
    call gh5%awrite(this%sym%ie,"SYMIE",1)
    call gh5%gclose(1)
    call gh5%dwrite(this%ek,this%nmax,this%kpt%dim,&
            &this%nspin,"/e_list",1)
    call gh5%dcreate(5,(/this%nmaxin,this%nmaxin,&
            &this%sym%nop,this%kpt%dim, &
            &this%nspin/),"/V_LIST",gh5%dcomplex_id,1,1)
    call gh5%dwrite(this%vk,this%nmaxin,this%nmaxin,this%sym%nop, &
            &this%kpt%diml,this%nspin,0,0,0,this%kpt%base,0,1)
    call gh5%dclose(1)
    call gh5%dwrite(this%r,loc%nasotot,loc%nasotot,this%nspin,"/BND_R",1)
    call gh5%fclose(1)
    return

    end subroutine h5write_bands


    subroutine calc_hamil_1k(hk,r,la1,nbands,nasot)
    integer,intent(in)::nbands,nasot
    complex(q),intent(in)::r(nasot,nasot),la1(nasot,nasot)
    complex(q),intent(inout)::hk(nbands,nbands)

    integer ia
    complex(q) hk_sub(nasot,nbands)

    hk_sub = hk(1:nasot,:)
    call annxb('n',r,hk_sub,nasot,nbands,hk(1:nasot,:))
    call anmxbmm('c',hk(:,1:nasot),r,nbands,nasot) ! rhr^+
    hk(1:nasot,1:nasot)=hk(1:nasot,1:nasot)+la1
    return

    end subroutine calc_hamil_1k


    !*************************************************************************
    !< calc total magnetic moment
    !*************************************************************************
    subroutine calc_mup_dn(this,io)
    class(bandstru_ob)::this
    integer,intent(in)::io

    integer ik
      
    if(this%nspin==1)return
    this%mup_dn=0
    do ik=1,this%kpt%dim
        this%mup_dn=this%mup_dn+sum(this%ferwe(:,ik,1)-this%ferwe(:,ik,2))
    enddo
    if(io>0)then
        write(io,'(" total magnetic moment: ",f8.3)')this%mup_dn
    endif
    return
      
    end subroutine calc_mup_dn
   

    ! total number of electrons
    subroutine calc_num_electrons(this,io)
    class(bandstru_ob)::this
    integer,intent(in)::io

    real(q) ne

    ne=sum(this%ferwe)
    if(io>0)then
        write(io,'(" total number of electrons in simulation: ",f0.5)')ne
    endif
    if(abs(ne-this%nelet)>1.e-6_q.and.this%kpt%ensemble==0)then
        if(io>0)then
            write(io,'(" warning: electron number =",f0.2," vs ",f0.2)')&
                    &this%nelet,ne
        endif
    endif
    this%nelet=ne
    return

    end subroutine calc_num_electrons


    subroutine calc_corr_ebwidth(this,io)
    class(bandstru_ob)::this
    integer,intent(in)::io
      
    integer ik,isp
    real(q) emin,emax
      
    if(io<0)return
    emin=100._q; emax=-100._q
    do ik=1,this%kpt%dim
        do isp=1,this%nspin_in
            emin=min(emin,this%ek(this%ne(2,ik,isp),ik,isp))
            emax=max(emax,this%ek(this%ne(3,ik,isp),ik,isp))
        enddo
    enddo
    write(io,'(" correlated block: emin/emix=",2f10.4)')emin,emax
    this%ebmin=emin; this%ebmax=emax
    this%ebwidth=emax-emin
    return
      
    end subroutine calc_corr_ebwidth
   

    subroutine set_kpt_icor(this,io)
    class(bandstru_ob) :: this
    integer,intent(in)::io
      
    if(this%kpt%ismear/=-5)return
    this%kpt%icor=1
    if(abs(this%kpt%delta)>100._q)this%kpt%icor=0
    this%kpt%cordeg=-1.d-6
    if(this%kpt%delta>100._q)then
        this%kpt%cordeg=-this%kpt%delta+100._q
    elseif(this%kpt%delta>0._q)then
        this%kpt%cordeg=-this%kpt%delta
    elseif(this%kpt%delta<-100._q)then
        this%kpt%cordeg= this%kpt%delta-100._q
    elseif(this%kpt%delta<0._q)then
        this%kpt%cordeg= this%kpt%delta
    endif
    if(abs(this%kpt%cordeg)>.01_q)this%kpt%cordeg=-1.d-6
    if(io>0)write(io,'(" bz integration with tetra method, icor=",i2)')&
            &this%kpt%icor
    if(this%kpt%cordeg<0._q.and.io>0)write(io,'(" equal occupancy of degenerate &
        &states, tol=",e8.1)')-this%kpt%cordeg
    return
      
    end subroutine set_kpt_icor
    

    subroutine gutz_fermi(this,io)
    class(bandstru_ob)::this
    integer,intent(in)::io
     
    if(this%kpt%ismear==-5)then
        call gutz_fermi_tetra_w2k(this,io)
    elseif(this%kpt%ismear==0.or.this%kpt%ismear==-1)then
        call get_ef_fun(this)
    else
        write(0,'(" ismear=",i2)')this%kpt%ismear
        stop ' error: unsupported ismear!'
    endif
    call adjust_efermi_bandgap(this,io)
    if(io>0)write(io,'(" gutz fermi level=",f16.8)')this%ef
    return
      
    end subroutine gutz_fermi


    subroutine get_ef_fun(this)
    class(bandstru_ob) :: this

    real(q) emin,emax
    real(q),parameter::eps=1.e-10_q

    if(this%kpt%ensemble==0)then ! Canonical ensemble
        this%ferwe=0
        if(abs(this%nelet-this%nelet_mott)<1.e-6_q)then
            this%ef=0; return
        endif
        emin=minval(this%ek); emax=maxval(this%ek)
        this%ef=rtbis(this,emin,emax,eps)
    else
        call set_fermi_weight(this,this%ef)
    endif
    return

    end subroutine get_ef_fun


    subroutine calc_bz_integration_correction(this)
    class(bandstru_ob) :: this

    select case(this%kpt%ismear)
    case(-1)
        call calc_entropy_fermi(this)
    case(0)
        call calc_correction_gauss(this)
    end select
    return

    end subroutine calc_bz_integration_correction


    subroutine calc_entropy_fermi(this)
    class(bandstru_ob) :: this

    integer ikp,isp,ispp,ib
    real(q) entr,focc,f1,f2,eint

    this%ets2=0
    do isp=1,this%nspin
        entr=0
        ispp=min(isp,this%ispin_in)
        do ikp=1,this%kpt%dim
            do ib=1,this%ne(1,ikp,ispp)
                focc=this%ferwe(ib,ikp,isp)/this%rspo/this%kpt%wt(ikp)
                f1=focc; f2=1._q-focc
                if(f1>0._q.and.f2>0._q)then
                    eint=f1*log(f1)+f2*log(f2)
                    entr=entr+eint*this%kpt%wt(ikp)*this%rspo
                endif
            enddo
        enddo
        this%ets2(isp)=this%kpt%delta*entr
    enddo
    return

    end subroutine calc_entropy_fermi


    subroutine calc_correction_gauss(this)
    class(bandstru_ob) :: this

    integer ikp,isp,ispp,ib
    real(q) eta,de

    this%ets2=0
    do isp=1,this%nspin
        eta=0
        ispp=min(isp,this%nspin_in)
        do ikp=1,this%kpt%dim
            do ib=1,this%ne(1,ikp,ispp)
                de=(this%ek(ib,ikp,isp)-this%ef)/this%kpt%delta
                de=de*de
                if(de<15._q)eta=eta+0.5*this%kpt%delta*exp(-de)* &
                        &this%kpt%wt(ikp)*this%rspo
            enddo
        enddo
        eta=-eta*2._q/sqrt(pi)
        this%ets2(isp)=eta/2._q
    enddo
    return

    end subroutine calc_correction_gauss


    function dif_nele_fun(this,mu)
    class(bandstru_ob) :: this
    real(q),intent(in)::mu
    real(q) dif_nele_fun

    real(q) nelet

    call set_fermi_weight(this,mu)
    nelet=sum(this%ferwe)
    dif_nele_fun=nelet-(this%nelet-this%nelet_mott)
    return

    end function dif_nele_fun


    function rtbis(this,x1,x2,tol)
    class(bandstru_ob) :: this
    real(q),intent(in) :: x1,x2,tol

    real(q) dx,xmid,fmid,f,rtbis
    integer j
    integer,parameter::jmax=500

    fmid=dif_nele_fun(this,x2)
    f=dif_nele_fun(this,x1)
    if(f*fmid.ge.0._q)then
        if(abs(fmid).lt.tol)then
            rtbis=x2
            return
        else if(abs(f).lt.tol)then
            rtbis=x1
            return
        endif
        write(0,'(" fmid and f=",2f16.8)')fmid,f
        stop ' error in rtbis: root must be bracketed for bisection.'
    else
        if(f.lt.0.)then
            rtbis=x1
            dx=x2-x1
        else
            rtbis=x2
            dx=x1-x2
        endif
        do j=1,jmax
            dx=dx*.5_q
            xmid=rtbis+dx
            fmid=dif_nele_fun(this,xmid)
            if(fmid.le.0._q)rtbis=xmid
            if(abs(fmid)<tol)then
                rtbis=xmid
                return
            endif
        enddo
        stop ' error in rtbis: too many bisections in rtbis'
    endif

    end function rtbis


    subroutine adjust_efermi_bandgap(this,io)
    class(bandstru_ob) :: this
    integer,intent(in)::io

    integer ihomo,i,nelet,nelet_mott
    real(q) ecmin,evmax

    ! check band gap
    if(this%kpt%ensemble/=0)then
        return
    endif
    nelet=nint(this%nelet)
    if(abs(nelet-this%nelet)>1.e-6_q)then
        ! noninteger filling.
        return
    endif
    nelet_mott=nint(this%nelet_mott)
    if(abs(nelet_mott-this%nelet_mott)>1.e-6_q)then
        stop "error: noninteger mott localized electrons."
    endif
    nelet=nelet-nelet_mott
    if(mod(nelet,this%rspo)/=0)return
    ihomo=nelet/this%rspo
    if(minval(this%ne(1,:,:))<=ihomo)then
        ! minimal bands case
        return
    endif
    ! conduction band minimum
    ecmin=minval(this%ek(ihomo+1,:,:))
    ! valence band maximum
    evmax=maxval(this%ek(ihomo,:,:))
    if(evmax<ecmin)then
        if(io>0)then
            write(io,'(" band gap = ",f12.8)')ecmin-evmax
            write(io,'(" original ef = ",f12.6," adjusted ef = ",f12.6)') &
                    &this%ef,(evmax+ecmin)/2
        endif
        this%ef=(evmax+ecmin)/2
    endif
    return

    end subroutine adjust_efermi_bandgap
      
    !*************************************************************************
    subroutine gutz_fermi_tetra_w2k(this,io)
    class(bandstru_ob) :: this
    integer,intent(in) :: io

    integer,parameter::nw=250000
    integer :: iw(nw)
    real(q),allocatable::eb(:,:,:),e_(:,:),weight_(:)
    integer,allocatable::nehelp(:,:)
    real(q)elecn,ef
    integer nkpt,nemax,jspin,iso,nbmax
    integer ik,isp,ispp,n1,n2
      
    nkpt=this%kpt%dim; nemax=this%nmax; iso=this%iso
    allocate(nehelp(nkpt,2),eb(nemax,nkpt,2))
    eb=3._q
    do isp=1,2
        ispp=min(isp,this%nspin_in)
        nehelp(:,isp)=this%ne(1,:,ispp)
        ispp=min(isp,this%nspin)
        do ik=1,nkpt
            eb(1:nehelp(ik,isp),ik,isp)=this%ek(1:nehelp(ik,isp),ik,ispp)
        enddo
    enddo
    elecn=this%nelet-1.d-10-this%nelet_mott
    jspin=this%ispin
    allocate(e_(nemax*jspin,nkpt)); e_=0
    e_(1:nemax,:) = eb(:,:,1)
    if(jspin==2)then
        e_(1+nemax:,:) = eb(:,:,2)
    endif

    allocate(weight_(nemax*jspin*nkpt)); weight_=0
    open(const_iu_kgen,file=this%kpt%file_name,status='old')
    ! consistent check
    read(const_iu_kgen,*)n1
    rewind(const_iu_kgen)
    if(n1/=nkpt)then
        this%ef=0
        return
    endif

    ! fermi.f90: dos, eweigh
    call dos(nemax*jspin,nkpt,e_,weight_,elecn/2.d0*iso*jspin,ef,iw,nw,&
            &this%kpt%icor,io)
    call eweigh(ef,weight_,nemax,nkpt,jspin,nehelp,eb,nbmax,this%kpt%cordeg,io)
    close(const_iu_kgen)
    do isp=1,this%nspin
        do ik=1,nkpt
            n1=(ik-1)*nemax*jspin+(isp-1)*nemax+1
            n2=(ik-1)*nemax*jspin+ isp   *nemax
            this%ferwe(:,ik,isp)=weight_(n1:n2)*this%rspo  ! convention
        enddo
    enddo
    this%ets2=0; this%ef=ef
    deallocate(eb,e_,weight_,nehelp)
    return
      
    end subroutine gutz_fermi_tetra_w2k
      
    !*************************************************************************
    subroutine bnd_modify_mott(this)
    class(bandstru_ob) :: this
    integer ik,i,isp,ispp,ntop
      
    if(this%n_mott==0)return
    do ik=1,this%kpt%dim
        do isp=1,this%nspin
            ispp=min(isp,this%nspin_in)
            ntop=this%ne(3,ik,ispp)
            do i=ntop-this%n_mott*this%iso/2+1,ntop
                if(abs(this%ek(i,ik,isp)-30._q)>1.e-6_q)then
                    write(0,'(" fetal error in bnd_modify_mott (!=30): &
                            &fixed band level = ",2f12.4)')this%ek(i,ik,isp)
                    stop
                endif
                this%ek(i,ik,isp)=this%ef
                this%ferwe(i,ik,isp)=real(this%nelet_mott,q)/this%n_mott* &
                        &this%rspo*this%kpt%wt(ik)
            enddo
        enddo
    enddo
    return
      
    end subroutine bnd_modify_mott
    

    subroutine set_fermi_weight(this,mu)
    class(bandstru_ob) :: this
    real(q) mu
      
    integer isp,ispp,ikp,ib
    real(q) dt
      
    this%ferwe=0
    do isp=1,this%nspin
        ispp=min(isp,this%nspin_in)
        do ikp=1,this%kpt%dim
            do ib=1,this%ne(1,ikp,ispp)
                dt=(this%ek(ib,ikp,isp)-mu)/this%kpt%delta
                select case(this%kpt%ismear)
                case(-1)
                    this%ferwe(ib,ikp,isp)=fermi_fun(dt)
                case(0)
                    this%ferwe(ib,ikp,isp)=gauss_fun(dt)
                case default
                    stop ' error in nelect_fun: this%kpt%ismear not supported!'
                end select
                this%ferwe(ib,ikp,isp)=this%ferwe(ib,ikp,isp)*this%rspo* &
                        &this%kpt%wt(ikp)
            enddo
        enddo
    enddo
    return
      
    end subroutine set_fermi_weight


      
end module bandstru
