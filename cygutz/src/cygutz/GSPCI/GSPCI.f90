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

module gspci
    ! sparse-matrix based CI, with impurity-bath basis.
    use gprec
    use gprimme
    use sparse
    use gutil
    use ghdf5_sm
    use gtime
    implicit none

    type gspci_mem
        integer::norb,norb2,nstates=0,nval_bot,nval_top,norb_mott=0, &
                &nelect_mott=0,imp=1, &
                &mode=4,  & ! 1: no constraints, 2: sz; 3: s2; 4:jz
                &iup=0,idn=1
        type(dcoo_matrix),pointer::mcoo(:,:),nvcoo(:,:),npcoo(:,:)
        integer,pointer::bs(:),ibs(:),idx(:),bs_l(:),ibs_l(:),idx_l(:), &
                &iorb_mott(:)=>null(),m_struct(:,:), &
                &mem_binom(:,:)
        type(ivector),allocatable::i_phi(:)
        real(q),pointer::bs_sz(:)=>null(),bs_sz_l(:)=>null(), &
                &sz_orb(:)=>null()
#ifdef real_version
        real(q),  &
#else
        complex(q),  &
#endif
                &pointer::h1e(:,:),daalpha(:,:),lambdac(:,:), &
                &v2e(:,:,:,:),dm(:,:),v(:)=>null()
        complex(q),pointer::svec(:,:,:)=>null(),lvec(:,:,:)=>null(), &
                &jvec(:,:,:)=>null(),jn(:,:)=>null()
#ifdef real_version
        type(dcsr_matrix)::ucsr
#else
        type(zcsr_matrix)::ucsr
#endif
        type(dcsr_matrix)::s2op
        real(q)::lambda_j2=2.0_q,etot
    end type gspci_mem

    type(gspci_mem)::dmem
    type(hdf5_sm_ob)::gh5
    type(time_ob)::time
    private
    public::gspci_mem,dmem,gh5,gspci_gh5exe,setup_hdns_spci,av1_gspci

    contains


    subroutine gspci_gh5exe()
    integer :: j,ierr,nv,idx,iembdv=0
    character*32 str
    real(q) etot,de
    logical lexist,lexist2
    character(255) cmd,gpath

    call time%start(1)
    call get_command(cmd)
    write(0, '(" command: ", a)')trim(cmd)
    call get_kwarg(cmd,'-i',dmem%imp)
    write(0,'(" solving impurity ", i2)')dmem%imp
    call get_kwarg(cmd,'-m',dmem%mode)
    ! negative for analysis
    select case(abs(dmem%mode))
        case (1)
            write(0,'(" with {N} contraints. ")')
        case (2)
            write(0,'(" with {N, Sz} contraints. ")')
        case (3)
            write(0,'(" with {N, S2} contraints. ")')
        case (4)
            write(0,'(" with {N, Jz} contraints. ")')
        case (201)
            write(0,'(" diagonal fock solver. ")')
        case default
            write(0,'(" mode = ", i0, " not defined!")')dmem%mode
            stop
    end select

    if(dmem%mode==3)then
        call get_kwarg(cmd,'-l',dmem%lambda_j2)
        write(0,'(" lambda_s2 = ", f0.4)')dmem%lambda_j2
    endif

    ! analysis mode?
    if(dmem%mode<0)then
        write(0,'(" gspci in post-analysis mode.")')
        iembdv=1
    else
        call get_kwarg(cmd,'-e',iembdv)
    endif

    if(iembdv>0)then
        write(0,'(" try existing eigen-vector serves as the starting point.")')
    else
        write(0,'(" random vector serves as the starting point.")')
    endif

    call gh5%init()
    gpath="/impurity_"//trim(int_to_str(dmem%imp-1))
    call gh5%fopen('HEmbed.h5',1,"r",serialio=.true.)
    call gh5%gopen(trim(gpath),1,1)
    call gh5%read(dmem%norb,'na2',1)
    call gh5%read(dmem%nval_bot,'nval_bot',1)
    call gh5%read(dmem%nval_top,'nval_top',1)
    call gh5%read(dmem%norb_mott,'norb_mott',1)

    write(0,'(" valence range from HEmbed.h5: [", i0, ", ", i0, "]")') &
            &dmem%nval_bot, dmem%nval_top
    if(dmem%norb_mott>0)then
        call  gh5%read(dmem%nelect_mott,'nelect_mott',1)
        allocate(dmem%iorb_mott(dmem%norb_mott))
        call gh5%read(dmem%iorb_mott(1),'iorb_mott',1)
    endif
    call gh5%gclose(1)
    allocate(dmem%m_struct(dmem%norb,dmem%norb))
    call gh5%read(dmem%m_struct,dmem%norb,dmem%norb,trim(gpath)// &
            &'/m_struct',1,serialio=.true.)
    dmem%m_struct=transpose(dmem%m_struct)
    dmem%norb2=dmem%norb*2

    allocate(dmem%h1e(dmem%norb,dmem%norb), &
            &dmem%daalpha(dmem%norb,dmem%norb), &
            &dmem%lambdac(dmem%norb,dmem%norb), &
            &dmem%v2e(dmem%norb,dmem%norb,dmem%norb,dmem%norb), &
            &dmem%dm(dmem%norb2,dmem%norb2))
    call gh5%read(dmem%h1e,dmem%norb,dmem%norb,trim(gpath)//'/H1E',1, &
            &serialio=.true.)
    call gh5%read(dmem%daalpha,dmem%norb,dmem%norb,trim(gpath)//'/D',1, &
            &serialio=.true.)
    call gh5%read(dmem%lambdac,dmem%norb,dmem%norb,trim(gpath)//'/LAMBDA',1, &
            &serialio=.true.)
    call gh5%read(dmem%v2e,dmem%norb,dmem%norb,dmem%norb,dmem%norb, &
            &trim(gpath)//'/V2E',1,serialio=.true.)
    call gh5%fclose(1)

    ! get jz or sz for the spin-orbitals
    call gh5%fopen('GParam.h5',1,'r')
    if(abs(dmem%mode)==2.or.abs(dmem%mode)==3.or.abs(dmem%mode)==4)then
        allocate(dmem%svec(dmem%norb,dmem%norb,1))
        call gh5%read(dmem%svec(:,:,1),dmem%norb, &
                &dmem%norb,trim(gpath)//'/sz',1, &
                &serialio=.true.)
        if(abs(dmem%mode)==4)then
            allocate(dmem%lvec(dmem%norb,dmem%norb,1))
            call gh5%read(dmem%lvec(:,:,1),dmem%norb,dmem%norb, &
                    &trim(gpath)//'/lz',1,serialio=.true.)
        else
            if(real(dmem%svec(1,1,1),q)<real(dmem%svec(2,2,1),q))then
                write(0,'(" in spin dn-up convention.")')
                dmem%iup=1
                dmem%idn=0
            else
                write(0,'(" in spin up-dn convention.")')
            endif
        endif
    endif
    call gh5%fclose(1)

    ! sz_orb: all 0 if mode==1; sz if mode==2/3, or jz if mode==4.
    allocate(dmem%sz_orb(dmem%norb))
    dmem%sz_orb=0
    do j=1,dmem%norb
        if(associated(dmem%svec))then
            dmem%sz_orb(j)=dmem%svec(j,j,1)
            if(associated(dmem%lvec))then
                dmem%sz_orb(j)=dmem%sz_orb(j)+dmem%lvec(j,j,1)
            endif
        endif
    enddo

    call init_gspci_mem()

    ! Check previous solution vector
    call gh5%fopen('HEmbed.h5',1,"r",serialio=.true.)
    gpath=trim(gpath)//"/ans"
    call gh5%exists(1,trim(gpath),lexist)
    if(lexist.and.(iembdv>0))then
        call gh5%gopen(trim(gpath),1,1)
        call gh5%exists(1,"dimv",lexist2)
        if(lexist2)then
            call gh5%read(nv,"dimv",1)
            if(nv==dmem%nstates)then
                allocate(dmem%v(nv))
                call gh5%read(dmem%v,nv,trim(gpath)//'/evec',1,serialio=.true.)
            else
                write(0,'(" existing solution vec not match in dimension,")')
                write(0,'(" random vector to be initialized!")')
            endif
        endif
        call gh5%gclose(1)
    endif
    call gh5%fclose(1)

    if(dmem%mode>0.or..not.associated(dmem%v))then
        write(0,'(" calculating evec.")')
        ! find ground state solution 
        call solve_hembed_spci_drive()
        call calc_dm_spci()
        if(dmem%mode==3)then
            call chk_eval_s2(.true.)
        endif
   
        de=trace_a(dmem%lambdac,dmem%norb) 
        dmem%etot=dmem%etot-de
        write(0,'(" emol =", f0.5)')dmem%etot
    
        call gh5%fopen('HEmbedRes_'//trim(int_to_str(dmem%imp-1))//'.h5',1,"w", &
                &serialio=.true.)
        call gh5%gopen("/",1,1)
        call gh5%awrite(dmem%etot,'emol',1)
        call gh5%awrite(dmem%nstates,'dimv',1)
        call gh5%gclose(1)
        call gh5%dwrite(dmem%dm,dmem%norb2,dmem%norb2,'/DM',1,serialio=.true.)
        call gh5%dwrite(dmem%v,dmem%nstates,'/evec',1,serialio=.true.)
        call gh5%fclose(1)
    endif

    if(dmem%mode<0)then
        ! post analysis of \rho or equivalently \Phi^\dagger \Phi
        call calc_save_rho_cp_blks()
    endif
    call gh5%end()
    return

    end subroutine gspci_gh5exe


    subroutine init_gspci_mem()
    
    allocate(dmem%mem_binom(dmem%norb,0:dmem%norb))
    call set_binoms(dmem%norb,dmem%mem_binom)
    call set_fock_state_indices(dmem%norb,dmem%mem_binom,dmem%idx,dmem%bs, &
            &dmem%ibs,dmem%bs_sz,dmem%sz_orb)
    call set_full_fock_states_l_spci()
    call set_mncoo_spci()
    call calc_ucsr_spci()
    if(dmem%mode==3)then
        call set_s2_spci()
    endif
    return

    end subroutine init_gspci_mem


    subroutine setup_hdns_spci(a,n)
    integer,intent(in)::n
#ifdef real_version
    real(q),intent(out)::a(n,n)
#else
    complex(q),intent(out)::a(n,n)
#endif

    call setup_hdns_spci_dlh(a,n)
    if(dmem%mode==3)then
        call add_hdns_spci_s2(a,n)
    endif
    return

    end subroutine setup_hdns_spci


    subroutine av1_gspci(v1,v2)
#ifdef real_version
    real(q),intent(in)::v1(*)
    real(q),intent(out)::v2(*)
#else
    complex(q),intent(in)::v1(*)
    complex(q),intent(out)::v2(*)
#endif

    call av1_gspci_dlh(v1,v2)
    if(dmem%mode==3)then
        call av1_gspci_s2(v1,v2)
    endif
    return

    end subroutine av1_gspci   


end module gspci
