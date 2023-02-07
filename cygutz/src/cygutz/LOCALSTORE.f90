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

module localstore
    use gprec, only: q
    use gutil, only: uhau,output_matrices,chk_eigens_matrix_list, &
            &orbital_spin_trans,int_to_str,annxb,trace_a, &
            &get_hm_expand,anmxbmm
    use corrorb
    use ghdf5
    implicit none
    private

    type mott_t
        integer :: nsorb,nelect
        integer,pointer :: idx_orb(:)
    end type mott_t
      
    type,public::matrix_basis
        integer :: dimhst=-1,dimhsmax=0,dimhs1123
        !< dimension of the matrix basis for each atom
        integer,pointer :: dim_hs(:)
#ifdef real_version
        real(q), &
#else
        complex(q), &
#endif
                &pointer :: hs(:)=>null() ! matrix basis
        integer,pointer :: m_struct(:)=>null() ! self-energy structure
    end type matrix_basis

    type,public :: localstore_ob
        integer :: num_imp
        integer,allocatable :: & 
                &imap_list(:), & !< equiv. impurity map index. imap_list(i)<=i.
                &na2_list(:), & !< impurity spin-orbital dimension list
                &naso_imp(:), &
                &nval_bot_list(:), &
                &nval_top_list(:)
        integer :: na2112,na2max,nasomax,nasotot,na2tot
        integer :: iso,ispin,ispin_in,ispo,rspo,riso,nspin, &
                &ivext=1  !-1: included in h1e and energy, but not in init_la1
                          ! 1: included in h1e and energy and init_la1
                          ! 0: initial perturbation.
        type(corrorb_ob),allocatable::co(:)
        type(mott_t),pointer::mott(:)

        real(q) :: edcla1,symerr=0,egamma_dc=0,shft_init_la1=0
        integer :: r_factor= &
#ifdef real_version
                &1
#else
                &2
#endif
#ifdef real_version
        real(q), &
#else
        complex(q), &
#endif
                &pointer :: r(:),r0(:),d(:),d0(:),la1(:),la2(:),z(:), &
                &nks(:),nc_var(:),nc_phy(:),h1e(:),isimix(:), &
                &vext(:)=>null(),r_coef(:),d0_coef(:),d_coef(:)
        complex(q),pointer :: sx(:)=>null(),sy(:)=>null(),sz(:)=>null(), & 
                &lx(:)=>null(),ly(:)=>null(),lz(:)=>null(), &
                ! additional rotations for bare hamiltonian.
                &db2sab(:)=>null() ! complex spherical harmonics basis
        real(q),pointer :: nks_coef(:),la1_coef(:),la2_coef(:),ncv_coef(:)
        !> hermitian matrix basis, (hermitian) maxtrix basis with
        !! selective mott localization.
        type (matrix_basis) :: hm, hm_l, hm_r

        contains
        procedure::init=>init_localstore
        procedure::read_h1e=>gh5_read_loc_h1e
        procedure::rotate_h1e=>rotate_h1e_list
        procedure::herm_matrices_pp=>calc_herm_matrices_pp
        procedure::rtrans=>site_wise_rtrans
        procedure::nks_pp=>calc_nks_pp
        procedure::calc_nks_tot
        procedure::init_rln=>init_loc_rln
        procedure::symm_check
        procedure::modify_rl_mott
        procedure::map_bnd_matrix=>map_loc_bnd_matrix
        procedure::eval_sl_vec_all
        procedure::calc_lambdac_list
        procedure::dump_hembed_list=>h5prepare_hembed_list
        procedure::load_hembed_list=>h5read_hembed_list
        procedure::calc_isimix
        procedure::r01_pp=>calc_r01_pp
        procedure::la1_pp=>calc_la1_pp
        procedure::ncvar_pp=>calc_ncvar_pp
        procedure::ncphy_pp=>calc_ncphy_pp
        procedure::calc_et1_list
        procedure::calc_la1_hf_list
        procedure::calc_edcla1
        procedure::write_results=>h5write_loc_results
        procedure::calc_da
        procedure::symm_dm_across_atoms
        procedure,private::hm_expand_all_sym
        procedure,private::hm_expand_all_herm
        procedure,private::hm_expand_all_general
        procedure::gh5_create_impurity_groups
        generic,public::hm_expand_all => hm_expand_all_sym, &
                &hm_expand_all_herm, &
                &hm_expand_all_general
    end type localstore_ob

    public::report_maxerr

    contains


    subroutine init_localstore(this,gh5,io)
    class(localstore_ob)::this
    class(hdf5_ob)::gh5
    integer,intent(in)::io

    logical lexist
  
    call gh5%fopen('GParam.h5',1,"r")
    ! open root group
    call gh5%gopen('/',1,1)
    !< Get num_imp, imap_list and na2_list
    call set_impurities_index(this,gh5,1)
    call write_impurities_index(this,io)
    call gh5_read_loc_hs(this,gh5,'/','matbs_dim_list','/matrix_basis', &
            &'/symbol_matrix',this%hm,1,1)
    call set_loc_db2sab(this,gh5,io,1)
    call set_loc_sl_vec(this,gh5,1)
    allocate(this%co(this%num_imp))
    call gh5_set_v2e_list(this,gh5,io,1)
    ! close root group
    call gh5%gclose(1)
    ! check external potential vext
    call gh5%exists(1,'/vext',lexist)
    if(lexist)then
        call set_loc_vext(this,gh5,io,1)
        call gh5%gopen("/vext",1,1)
        call gh5%read(this%ivext,'givext',1)
        call gh5%gclose(1)
    endif
    ! check mott parameters
    call gh5%exists(1,'/mott',lexist)
    if(lexist)then
        call gh5%gopen('/mott',1,1)
    endif
    call set_loc_hs(this,gh5,'/mott','matbs_r_dim_list',"/matrix_basis_r", &
            &"/symbol_matrix_r",this%hm_r,this%hm,lexist,1,1,0)
    call set_loc_hs(this,gh5,'/mott','matbs_l_dim_list',"/matrix_basis_l", &
            &"/symbol_matrix_l",this%hm_l,this%hm,lexist,1,1,1)
    if(lexist)then
        call gh5%gclose(1)
    endif
    ! for the case of real version without mott group.
    if(this%hm_r%dimhst<0)then
        call gh5%gopen('/',1,1)
        call set_loc_hs(this,gh5,'/','matbs_r_dim_list',"/matrix_basis_r", &
                &"/symbol_matrix",this%hm_r,this%hm,.true.,1,1,1)
        call gh5%gclose(1)
    endif

    call set_mott_localized_orbitals(this,gh5,lexist,1)
    call gh5%fclose(1)
    call output_mott_localized_orbitals(this,io)
    call alloc_localstore(this)
    call link_co_localstore(this)
    call set_diagonal_r(this,1._q)
    return

    end subroutine init_localstore


    !< Impurities should be ordered according to types.
    subroutine gh5_set_v2e_list(this,gh5,io,i_f)
    class(localstore_ob)::this
    class(hdf5_ob)::gh5
    integer,intent(in)::io,i_f

    integer::i,imap,na2

    do i=1,this%num_imp
        imap=this%imap_list(i)
        if(imap==i)then
            na2=this%na2_list(i)
            allocate(this%co(i)%v2e(na2,na2,na2,na2))
            call gh5%read(this%co(i)%v2e,na2,na2,na2,na2,'/impurity_'// &
                    &trim(int_to_str(i-1))//'/V2E',i_f)
            if(io>0)then
                write(io,'(" imp = ",I2," v2e(1,1,1,1) = ",f0.2)')i, &
                        &real(this%co(i)%v2e(1,1,1,1))
            endif
            call this%co(i)%set_v_j2e(na2)
        else
            this%co(i)%v2e=>this%co(imap)%v2e
            this%co(i)%v_j2e=>this%co(imap)%v_j2e
        endif
    enddo
    return

    end subroutine gh5_set_v2e_list


    subroutine calc_la1_hf_list(this)
    class(localstore_ob)::this

    integer i

    do i=1,this%num_imp
        call this%co(i)%calc_la1_hf()
    enddo
    return

    end subroutine calc_la1_hf_list


    subroutine link_co_localstore(this)
    class(localstore_ob)::this

    integer i,na2,na22,naso,nasoo,dim_hs,ibase,jlbase,jrbase

    ibase=0; jlbase=0; jrbase=0
    do i=1,this%num_imp
        na2=this%na2_list(i)
        na22=na2*na2
        this%co(i)%dim2=na2; this%co(i)%dim=na2/2; this%co(i)%dim4=na2*2
        naso=na2*this%iso/2
        this%co(i)%dimso=naso
        nasoo=this%co(i)%dimso**2
        this%co(i)%dim_hs_l = this%hm_l%dim_hs(i)
        this%co(i)%dim_hs_r = this%hm_r%dim_hs(i)
        this%co(i)%nks   (1:na2,1:na2) => this%nks   (ibase+1:ibase+na22)
        this%co(i)%isimix(1:na2,1:na2) => this%isimix(ibase+1:ibase+na22)
        this%co(i)%nc_var(1:na2,1:na2) => this%nc_var(ibase+1:ibase+na22)
        this%co(i)%nc_phy(1:na2,1:na2) => this%nc_phy(ibase+1:ibase+na22)
        this%co(i)%r     (1:na2,1:na2) => this%r     (ibase+1:ibase+na22)
        this%co(i)%r0    (1:na2,1:na2) => this%r0    (ibase+1:ibase+na22)
        this%co(i)%z     (1:na2,1:na2) => this%z     (ibase+1:ibase+na22)
        this%co(i)%d0    (1:na2,1:na2) => this%d0    (ibase+1:ibase+na22)
        this%co(i)%d     (1:na2,1:na2) => this%d     (ibase+1:ibase+na22)
        this%co(i)%la1   (1:na2,1:na2) => this%la1   (ibase+1:ibase+na22)
        this%co(i)%la2   (1:na2,1:na2) => this%la2   (ibase+1:ibase+na22)
        this%co(i)%h1e   (1:na2,1:na2) => this%h1e   (ibase+1:ibase+na22)

        if(associated(this%vext))then
            this%co(i)%vext(1:na2,1:na2) => this%vext(ibase+1:ibase+na22)
        endif

        if(associated(this%sx))then
            this%co(i)%sx(1:na2,1:na2) => this%sx(ibase+1:ibase+na22)
            this%co(i)%sy(1:na2,1:na2) => this%sy(ibase+1:ibase+na22)
        endif
        if(associated(this%sz))then
            this%co(i)%sz(1:na2,1:na2) => this%sz(ibase+1:ibase+na22)
        endif
        if(associated(this%lx))then
            this%co(i)%lx(1:na2,1:na2) => this%lx(ibase+1:ibase+na22)
            this%co(i)%ly(1:na2,1:na2) => this%ly(ibase+1:ibase+na22)
        endif
        if(associated(this%lz))then
            this%co(i)%lz(1:na2,1:na2) => this%lz(ibase+1:ibase+na22)
        endif
        if(associated(this%db2sab))then
            this%co(i)%db2sab(1:na2,1:na2) => &
                    &this%db2sab(ibase+1:ibase+na22)
        endif

        this%co(i)%m_struct(1:na2,1:na2) => &
                &this%hm%m_struct(ibase+1:ibase+na22)

        dim_hs=this%hm_l%dim_hs(i)
        if(dim_hs>0)then
            this%co(i)%hs_l(1:na2,1:na2,1:dim_hs) => this%hm_l%hs(jlbase+1: &
                    &jlbase+na22*dim_hs)
        endif
        jlbase=jlbase+na22*dim_hs
        dim_hs=this%hm_r%dim_hs(i)
        if(dim_hs>0)then
            this%co(i)%hs_r(1:na2,1:na2,1:dim_hs) => this%hm_r%hs(jrbase+1: &
                    &jrbase+na22*dim_hs)
        endif
        jrbase=jrbase+na22*dim_hs
        ibase=ibase+na22
    enddo
    return

    end subroutine link_co_localstore


    ! \sum_ij \lambda_ij <c^\dagger_i c_j>
    subroutine calc_edcla1(this)
    class(localstore_ob)::this
    
    this%edcla1=real(sum(this%la1*this%nks),q)
    return

    end subroutine calc_edcla1


    subroutine set_loc_hs_dimtot(this,mb)
    class(localstore_ob)::this
    type(matrix_basis),intent(inout)::mb

    integer i

    mb%dimhsmax=maxval(mb%dim_hs)
    mb%dimhst=0
    mb%dimhs1123=0
    do i=1,this%num_imp
        mb%dimhs1123=mb%dimhs1123+this%na2_list(i)**2*mb%dim_hs(i)
        if(this%imap_list(i)==i)then
            mb%dimhst=mb%dimhst+mb%dim_hs(i)
        endif
    enddo
    return

    end subroutine set_loc_hs_dimtot


    subroutine set_impurities_index(this,gh5,i)
    class(localstore_ob)::this
    class(hdf5_ob)::gh5
    integer,intent(in)::i

    logical lexist

    ! read attributes 
    call gh5%read(this%iso, 'iso', i)
    call gh5%read(this%ispin, 'ispin', i)
    call gh5%read(this%num_imp, 'num_imp', i)
    allocate(this%imap_list(this%num_imp), &
            &this%na2_list(this%num_imp),this%naso_imp(this%num_imp))
    call gh5%read(this%imap_list(1), 'imap_list', i)
    call gh5%read(this%na2_list(1), 'na2_list', i)

    ! zero-base to one-base
    this%imap_list = this%imap_list + 1

    this%ispo=max(this%iso,this%ispin)
    this%rspo=3-this%ispo
    this%riso=3-this%iso
    this%nspin=max(1,this%ispin/this%iso)
    this%na2max=maxval(this%na2_list)
    this%nasomax=this%na2max*this%iso/2
    this%na2112=sum(this%na2_list**2)
    this%na2tot=sum(this%na2_list)
    this%nasotot=this%na2tot*this%iso/2
    this%naso_imp=this%na2_list*this%iso/2

    call gh5%exists(i,'nval_bot_list',lexist)
    if(lexist)then
        allocate(this%nval_bot_list(this%num_imp),  &
                &this%nval_top_list(this%num_imp))
        call gh5%read(this%nval_bot_list(1), 'nval_bot_list', i)
        call gh5%read(this%nval_top_list(1), 'nval_top_list', i)
    endif
    return

    end subroutine set_impurities_index


    subroutine write_impurities_index(this,io)
    class(localstore_ob)::this
    integer,intent(in)::io

    if(io<0)return
    write(io,'(" total number of impurities = ",I4)')this%num_imp
    write(io,'(" impurity imap indices:")')
    write(io,'(4x,10I4)')this%imap_list
    write(io,'(" impurity num_spin_orbitals:")')
    write(io,'(4x,10I4)')this%na2_list
    return

    end subroutine write_impurities_index


    subroutine alloc_localstore(this)
    class(localstore_ob)::this
      
    allocate(this%r(this%na2112), this%r0(this%na2112), &
            &this%z(this%na2112), this%d(this%na2112), &
            &this%d0(this%na2112), this%la1(this%na2112), &
            &this%la2(this%na2112), this%nks(this%na2112), &
            &this%isimix(this%na2112), this%nc_var(this%na2112), &
            &this%nc_phy(this%na2112), this%h1e(this%na2112), &
            &this%nks_coef(this%hm_l%dimhst), &
            &this%la1_coef(this%hm_l%dimhst), &
            &this%la2_coef(this%hm_l%dimhst), &
            &this%ncv_coef(this%hm_l%dimhst), &
            &this%r_coef(this%hm_r%dimhst), &
            &this%d_coef(this%hm_r%dimhst), &
            &this%d0_coef(this%hm_r%dimhst))
    this%r=0; this%r0=0; this%z=0; this%d=0; this%d0=0; this%la1=0; this%la2=0
    this%nks=0; this%isimix=0; this%nc_var=0; this%nc_phy=0
    return
      
    end subroutine alloc_localstore
    

    subroutine h5write_loc_results(this,gh5,i_g,i_f)
    class(localstore_ob)::this
    class(hdf5_ob)::gh5
    integer,intent(in)::i_g,i_f

    real(q)::rbuf(this%num_imp)

    call gh5_create_impurity_groups(this,gh5,0,i_f)
    call gh5_wrt_loc_matrix_list(this,gh5,'/NKS',this%nks,i_f)
    call gh5_wrt_loc_matrix_list(this,gh5,'/D0',this%d0,i_f)
    call gh5_wrt_loc_matrix_list(this,gh5,'/H1E',this%h1e,i_f)
    call gh5_wrt_loc_matrix_list(this,gh5,'/LA1',this%la1,i_f)
    call gh5_wrt_loc_matrix_list(this,gh5,'/D',this%d,i_f)
    call gh5_wrt_loc_matrix_list(this,gh5,'/R',this%r,i_f)
    call gh5_wrt_loc_matrix_list(this,gh5,'/NC_PHY',this%nc_phy,i_f)
    rbuf=this%co(:)%eu2
    call gh5%awrite(rbuf,this%num_imp,"interaction energies",i_g)
    rbuf=this%co(:)%et1
    call gh5%awrite(rbuf,this%num_imp,"one-body energies",i_g)
    call gh5%awrite(this%edcla1,"lambda-dc energy",i_g)
    call gh5%awrite(this%egamma_dc,"egamma_dc",i_g)
    return

    end subroutine h5write_loc_results


    subroutine init_loc_rln(this,gh5,io,r_default,iembeddiag)
    class(localstore_ob)::this
    class(hdf5_ob)::gh5
    integer,intent(in)::io,iembeddiag
    real(q),intent(in)::r_default

    logical lexist

    inquire(file='GLog.h5',exist=lexist)
    if(lexist)then
        call gh5_read_loc_rln(this,gh5,r_default,iembeddiag)
    else
        this%la1=this%h1e
        call shft_diagonal_la1(this,this%shft_init_la1)
        call set_diagonal_r(this,r_default)
        ! initial symmetry break term
        if(this%ivext>=0.and.associated(this%vext))then
            this%la1=this%la1+this%vext
        endif
    endif
    if(iembeddiag==10)then
        ! symmetry breaking
        if(this%ivext==0.and.associated(this%vext))then
            this%nks=this%nks-this%vext
        endif
        call calc_nks_pp(this,io)
    else
        call modify_rl_mott(this,0,30._q)
        call calc_r_pp(this,io,'r-inp')
        call calc_la1_pp(this,io,'la1-inp')
    endif
    return

    end subroutine init_loc_rln


    subroutine shft_diagonal_la1(this,delta)
    class(localstore_ob)::this
    real(q),intent(in)::delta

    integer i,j

    do i=1,this%num_imp
        do j=1,this%na2_list(i)
            this%co(i)%la1(j,j)=this%co(i)%la1(j,j)+delta
        enddo
    enddo
    return

    end subroutine shft_diagonal_la1


    subroutine set_diagonal_r(this,r_default)
    class(localstore_ob)::this
    real(q),intent(in)::r_default

    integer i,j

    do i=1,this%num_imp
        do j=1,this%na2_list(i)
            this%co(i)%r(j,j)=r_default
        enddo
    enddo
    return

    end subroutine set_diagonal_r


    !*************************************************************************
    subroutine gh5_wrt_loc_matrix_list(this,gh5,apath,a,i_f)
    class(localstore_ob)::this
    class(hdf5_ob)::gh5
    character(*),intent(in)::apath
#ifdef real_version
    real(q),  &
#else
    complex(q),  &
#endif
            &target,intent(in)::a(this%na2112)
    integer,intent(in)::i_f
      
    integer::i,na2,na22,ibase
#ifdef real_version
    real(q),  &
#else
    complex(q),  &
#endif
            &target::a_buf(this%na2max**2)
#ifdef real_version
    real(q),  &
#else
    complex(q),  &
#endif
            &pointer::p_a(:,:)
      
    ibase=0
    do i = 1,this%num_imp
        na2=this%na2_list(i)
        na22=na2**2
        p_a(1:na2,1:na2)=>a_buf(1:na2**2)
        a_buf=a(ibase+1:ibase+na22)
        if(apath(2:2)>="a".and.apath(2:2)<="z")then
            ! fortran-order to c-order transfer
            p_a=transpose(p_a)
        endif
        call gh5%dwrite(p_a,na2,na2,&
                &'/impurity_'//trim(int_to_str(i-1))//apath,i_f)
        ibase=ibase+na22
    enddo
    nullify(p_a)
    return
      
    end subroutine gh5_wrt_loc_matrix_list


    subroutine calc_isimix(this,io,mode)
    class(localstore_ob)::this
    integer,intent(in)::io,mode

    integer i,imap

    do i=1,this%num_imp
        imap=this%imap_list(i)
        if(imap==i)then
            call this%co(i)%calc_isimix(mode)
        else
            this%co(i)%isimix=this%co(imap)%isimix
        endif
    enddo
    call output_matrices('isimix-var',this%isimix,this%na2112, &
            &this%num_imp,this%na2_list,io,0)
    return

    end subroutine calc_isimix


    subroutine eval_sl_vec_all(this,mode,io)
    class(localstore_ob)::this
    integer,intent(in)::mode,io

    integer i
    character name_*12

    if(.not.associated(this%sx))return
    if(mode==1)then
        name_='var-hf-only'
    else
        name_='physical'
    endif
    do i=1,this%num_imp
        call this%co(i)%eval_sl_vec(mode)
        if(io>0)then
            write(io,'(" imp =",i4," s_xyz(",a12,") =",3f12.5)')i,name_, &
                &this%co(i)%s_val(:,mode)
            write(io,'(10x," l_xyz(",a12,") =",3f12.5)')name_, &
                &this%co(i)%l_val(:,mode)
        endif
    enddo
    return

    end subroutine eval_sl_vec_all


    subroutine gh5_read_loc_h1e(this,gh5,nspin_in,i_f)
    class(localstore_ob)::this
    class(hdf5_ob)::gh5
    integer,intent(in)::nspin_in,i_f

    integer i,naso,na22,ibase
    complex(q),allocatable::h1e(:,:,:,:)  ! the input is always complex

    this%h1e=0
    allocate(h1e(this%nasomax,this%nasomax,this%num_imp, &
            &max(nspin_in/this%iso,1)))
    call gh5%read(h1e(:,:,:,1),this%nasomax,this%nasomax,this%num_imp, &
            &"/ispin_0/H1E_LIST",i_f)
    if(this%iso==1.and.nspin_in==2)then
        call gh5%read(h1e(:,:,:,2),this%nasomax,this%nasomax,this%num_imp, &
                &"/ispin_1/H1E_LIST",i_f)
    endif

    do i=1,this%num_imp
        naso=this%naso_imp(i)
        this%co(i)%h1e(1:naso,1:naso)=h1e(1:naso,1:naso,i,1)
        if(this%iso==1)then
            if(nspin_in==2)then
                this%co(i)%h1e(1+naso:,1+naso:)=h1e(1:naso,1:naso,i,2)
            else
                this%co(i)%h1e(1+naso:,1+naso:)=h1e(1:naso,1:naso,i,1)
            endif
        endif
    enddo
    deallocate(h1e)
    return
      
    end subroutine gh5_read_loc_h1e


    subroutine rotate_h1e_list(this)
    class(localstore_ob)::this

    integer i,na2,na
    complex(q),pointer::p_h1e(:,:)
    complex(q),target::h1e(this%na2max**2)

    do i=1,this%num_imp
        na2=this%co(i)%dim2; na=this%co(i)%dim
        p_h1e(1:na2,1:na2) => h1e(1:na2**2)
        p_h1e=this%co(i)%h1e
        if(this%iso==1.and.(.not.associated(this%db2sab)))then
            ! {{orbs}_up, {orbs}_dn} -> {{orb_up, orb_dn}}
            this%co(i)%h1e(1::2,1::2)=p_h1e(1:na,1:na)
            this%co(i)%h1e(2::2,2::2)=p_h1e(1+na:,1+na:)
            this%co(i)%h1e(1::2,2::2)=0
            this%co(i)%h1e(2::2,1::2)=0
        else
            call uhau(p_h1e,this%co(i)%db2sab,na2,na2)
            this%co(i)%h1e=p_h1e
        endif
    enddo

    nullify(p_h1e)
    return

    end subroutine rotate_h1e_list


    subroutine gh5_read_loc_matrix_list(this,gh5,prefix,path,dim_imp,a,i_f)
    class(localstore_ob)::this
    class(hdf5_ob)::gh5
    character(*),intent(in)::prefix,path
    integer,intent(in)::dim_imp(*),i_f
#ifdef real_version
    real(q),  &
#else
    complex(q),  &
#endif
            &target,intent(out)::a(*)
      
    integer i,na2,na22,ibase
#ifdef real_version
    real(q),  &
#else
    complex(q),  &
#endif
            &pointer::p_a(:,:)
      
    ibase=0
    do i = 1,this%num_imp
        na2=dim_imp(i)
        na22=na2**2
        p_a(1:na2,1:na2)=>a(ibase+1:ibase+na22)
        call gh5%read(p_a,na2,na2,prefix//'/impurity_'// &
                &trim(int_to_str(i-1))//path,i_f)
        if(path(2:2)>="a".and.path(2:2)<="z")then
            ! path in lower case: c-order to fortran-order transfer.
            p_a(1:na2,1:na2)=transpose(p_a(1:na2,1:na2))
        endif
        ibase=ibase+na22
    enddo
    nullify(p_a)
    return
      
    end subroutine gh5_read_loc_matrix_list
  

    subroutine gh5_read_loc_zmatrix_list(this,gh5,prefix,path,dim_imp,a,i_f)
    class(localstore_ob)::this
    class(hdf5_ob)::gh5
    character(*),intent(in)::prefix,path
    integer,intent(in)::dim_imp(*),i_f
    complex(q),target,intent(out)::a(*)
      
    integer i,na2,na22,ibase
    complex(q),pointer::p_a(:,:)
      
    ibase=0
    do i = 1,this%num_imp
        na2=dim_imp(i)
        na22=na2**2
        p_a(1:na2,1:na2)=>a(ibase+1:ibase+na22)
        call gh5%read(p_a,na2,na2,prefix//'/impurity_'// &
                &trim(int_to_str(i-1))//path,i_f)
        if(path(2:2)>="a".and.path(2:2)<="z")then
            ! path in lower case: c-order to fortran-order transfer.
            p_a(1:na2,1:na2)=transpose(p_a(1:na2,1:na2))
        endif
        ibase=ibase+na22
    enddo
    nullify(p_a)
    return
      
    end subroutine gh5_read_loc_zmatrix_list
    

    !*************************************************************************
    subroutine gh5_read_loc_rln(this,gh5,r_default,iembeddiag)
    class(localstore_ob)::this
    class(hdf5_ob)::gh5
    real(q),intent(in)::r_default
    integer,intent(in)::iembeddiag

    logical lexist

    call gh5%fopen('GLog.h5',1,"r")
    if(iembeddiag==10)then
        ! hartree-fock
        call gh5_read_loc_matrix_list(this,gh5,'/','/NKS', &
                &this%na2_list,this%nks,1)
    else
        call gh5_read_loc_matrix_list(this,gh5,'/','/LA1', &
                &this%na2_list,this%la1,1)
    endif
    call gh5%exists(1,'/impurity_0/R',lexist)
    if(lexist.and.iembeddiag/=10)then
        call gh5_read_loc_matrix_list(this,gh5,'/','/R', &
                &this%na2_list,this%r,1)
    else
        call set_diagonal_r(this,r_default)
    endif
    call gh5%fclose(1)
    return
      
    end subroutine gh5_read_loc_rln
    

    !*************************************************************************
    ! read [s_x, s_y, s_z] and [l_x, l_y, l_z]
    !*************************************************************************
    subroutine set_loc_sl_vec(this,gh5,i_f)
    class(localstore_ob)::this
    class(hdf5_ob)::gh5
    integer,intent(in)::i_f

    logical lexist

    call gh5%exists(i_f,'/impurity_0/sz',lexist)
    if(lexist)then
        allocate(this%sz(this%na2112))
        call gh5_read_loc_zmatrix_list(this,gh5,'/','/sz', &
                &this%na2_list,this%sz,i_f)
    endif
    call gh5%exists(i_f,'/impurity_0/lz',lexist)
    if(lexist)then
        allocate(this%lz(this%na2112))
        call gh5_read_loc_zmatrix_list(this,gh5,'/','/lz', &
                &this%na2_list,this%lz,i_f)
    endif
    call gh5%exists(i_f,'/impurity_0/sx',lexist)
    if(lexist)then
        allocate(this%sx(this%na2112),this%sy(this%na2112), &
                &this%lx(this%na2112),this%ly(this%na2112))
        call gh5_read_loc_zmatrix_list(this,gh5,'/','/sx', &
                &this%na2_list,this%sx,i_f)
        call gh5_read_loc_zmatrix_list(this,gh5,'/','/sy', &
                &this%na2_list,this%sy,i_f)
        call gh5_read_loc_zmatrix_list(this,gh5,'/','/lx', &
                &this%na2_list,this%lx,i_f)
        call gh5_read_loc_zmatrix_list(this,gh5,'/','/ly', &
                &this%na2_list,this%ly,i_f)
    endif
    return
      
    end subroutine set_loc_sl_vec
    

    subroutine gh5_read_loc_hs(this,gh5,prefix,apath,hpath,spath,mb,i_g,i_f)
    class(localstore_ob)::this
    class(hdf5_ob)::gh5
    character(*),intent(in)::prefix,apath,hpath,spath
    type(matrix_basis),intent(inout)::mb
    integer,intent(in)::i_g,i_f

    integer na2,dim_hs,na223,i,j,ibase,jbase
#ifdef real_version
    real(q),  &
#else
    complex(q),  &
#endif
            &pointer::mb_hs(:,:,:)
    integer,pointer::mb_m_struct(:,:)

    allocate(mb%dim_hs(this%num_imp))
    call gh5%read(mb%dim_hs(1),apath,i_g)
    call set_loc_hs_dimtot(this,mb)
    allocate(mb%hs(mb%dimhs1123),mb%m_struct(this%na2112))
    mb%m_struct=0
    ibase=0; jbase=0
    do i=1,this%num_imp
        na2=this%na2_list(i)
        dim_hs=mb%dim_hs(i)
        na223=na2**2*dim_hs
        if(na223>0)then
            mb_hs(1:na2,1:na2,1:dim_hs)=>mb%hs(ibase+1:ibase+na223)
            call gh5%read(mb_hs,na2,na2,dim_hs,prefix//'/impurity_'// &
                    &trim(int_to_str(i-1))//hpath,i_f)
            ! hpath(1:1) is "/".
            if(hpath(2:2)>="a".and.hpath(2:2)<="z")then
                ! switch axis due to c-fortran convention.
                do j=1,dim_hs
                    mb_hs(:,:,j)=transpose(mb_hs(:,:,j))
                enddo
            endif
            ibase=ibase+na223
            mb_m_struct(1:na2,1:na2)=>mb%m_struct(jbase+1:jbase+na2**2)
            call gh5%read(mb_m_struct,na2,na2,prefix//'/impurity_'// &
                    &trim(int_to_str(i-1))//spath,i_f)
            if(spath(2:2)>="a".and.spath(2:2)<="z")then
                ! switch axis due to c-fortran convention.
                mb_m_struct=transpose(mb_m_struct)
            endif
        endif
        jbase=jbase+na2**2
    enddo
    nullify(mb_hs,mb_m_struct)
    return
      
    end subroutine gh5_read_loc_hs
     

    !*************************************************************************
    subroutine set_loc_hs(this,gh5,prefix,dpath,hpath,spath,mb,mbt, &
            &lopen,i_g,i_f,mode)
    class(localstore_ob)::this
    class(hdf5_ob)::gh5
    character(*),intent(in)::prefix,dpath,hpath,spath
    type(matrix_basis),intent(inout)::mb
    type(matrix_basis),intent(in)::mbt
    logical,intent(in)::lopen
    integer,intent(in)::i_g,i_f,mode

    logical lexist
     
    if(lopen)then
        call gh5%exists(i_f,dpath,lexist)
    else
        lexist=.false.
    endif
    if(lexist)then
        call gh5_read_loc_hs(this,gh5,prefix,dpath,hpath,spath,mb,i_g,i_f)
    elseif(mode>0)then
        mb%dim_hs=>mbt%dim_hs
        mb%hs=>mbt%hs
        mb%dimhsmax=mbt%dimhsmax
        mb%dimhst=mbt%dimhst
        mb%dimhs1123=mbt%dimhs1123
    endif
    return
      
    end subroutine set_loc_hs


    ! Convention: 
    ! {{complext spherical harmonics (CSH)}_dn,up}
    ! to symmetry adapted basis (SAB)
    subroutine set_loc_db2sab(this,gh5,io,i_f)
    class(localstore_ob)::this
    class(hdf5_ob)::gh5
    integer,intent(in)::io,i_f

    logical lexist
   
    call gh5%exists(i_f,'/impurity_0/db2sab',lexist)
    if(lexist)then
        allocate(this%db2sab(this%na2112))
        call gh5_read_loc_zmatrix_list(this,gh5,'/','/db2sab', &
                &this%na2_list,this%db2sab,i_f)
        if(io>0)then
            write(io,'(" transformation db_to_sab read in.")')
        endif
    endif
    return
      
    end subroutine set_loc_db2sab
     

    ! read in a list of local external potentials.
    subroutine set_loc_vext(this,gh5,io,i_f)
    class(localstore_ob)::this
    class(hdf5_ob)::gh5
    integer,intent(in)::io,i_f

    allocate(this%vext(this%na2112))
    call gh5_read_loc_matrix_list(this,gh5,'/vext','/v', &
            &this%na2_list,this%vext,i_f)
    call output_matrices('v_ext',this%vext,this%na2112,this%num_imp, &
            &this%na2_list,io,0)
    return
 
    end subroutine set_loc_vext


    !*************************************************************************
    ! Mott localized spin-orbital information
    !*************************************************************************
    subroutine set_mott_localized_orbitals(this,gh5,lopen,i_f)
    class(localstore_ob)::this
    class(hdf5_ob)::gh5
    logical,intent(in)::lopen
    integer,intent(in)::i_f

    integer i
    logical lexist

    allocate(this%mott(this%num_imp))
    if(lopen)then
        call gh5%exists(i_f,'/mott/impurity_0',lexist)
    else
        lexist=.false.
    endif
    if(.not.lexist)then
        this%mott(:)%nsorb=0; this%mott(:)%nelect=0
    else
        do i=1,this%num_imp
            call gh5%gopen('/mott/impurity_'//trim(int_to_str(i-1)),i_f,1)
            call gh5%read(this%mott(i)%nsorb,'num_mott_orbitals',1)
            call gh5%read(this%mott(i)%nelect,'num_mott_electrons',1)
            if(this%mott(i)%nsorb>0)then
                allocate(this%mott(i)%idx_orb(this%mott(i)%nsorb))
                call gh5%read(this%mott(i)%idx_orb(1),'mott_orbital_indices',1)
                this%mott(i)%idx_orb=this%mott(i)%idx_orb+1
            endif
            call gh5%gclose(1)
        enddo
    endif
    return
      
    end subroutine set_mott_localized_orbitals
    

    !*************************************************************************
    subroutine output_mott_localized_orbitals(this,io)
    class(localstore_ob)::this
    integer,intent(in)::io
      
    integer i
      
    if(io<0)return
    write(io,'(" mott localized orbital info:")')
    do i=1,this%num_imp
        write(io,'(" i=",i3," nsorb=",i2," nelect=",i2)')i, &
                &this%mott(i)%nsorb,this%mott(i)%nelect
        if(this%mott(i)%nsorb>0)then
            write(io,'("    idx_orb=",14i3)')this%mott(i)%idx_orb
        endif
    enddo
    return
      
    end subroutine output_mott_localized_orbitals
    

    !*************************************************************************
    subroutine modify_rl_mott(this,mode,val)
    class(localstore_ob)::this
    integer,intent(in)::mode
    real(q),intent(in)::val

    integer i,j,j_
      
    do i=1,this%num_imp
        do j=1,this%mott(i)%nsorb
            j_=this%mott(i)%idx_orb(j)
            this%co(i)%r(j_,:)=0
            this%co(i)%la1(j_,:)=0; this%co(i)%la1(:,j_)=0
            if(mode>0)then
                this%co(i)%la1(j_,j_)=val
            endif
        enddo
    enddo
    return
      
    end subroutine modify_rl_mott


    !*************************************************************************
    !< post processing: print and symmetizations.
    !! nks could be f(H) or f'(H), with later the derivative of fermi fun.
    !*************************************************************************
    subroutine calc_nks_pp(this,io,lchk,obj)
    class(localstore_ob) :: this
    integer,intent(in)::io
    logical,optional,intent(in)::lchk
    character(3),optional,intent(in)::obj

    real(q) maxerr
#ifdef real_version
    real(q)  &
#else
    complex(q)  &
#endif
            &copy(this%na2112)
    character(3)::ob="nks"

    if(present(obj))then
        ob=obj
    endif
    call output_matrices(ob//'-unsym',this%nks,this%na2112,this%num_imp, &
            &this%na2_list,io,1)
    call chk_eigens_matrix_list('nks-unsym',this%nks,this%na2112,this%num_imp, &
            &this%na2_list,io,0)
    copy=this%nks
    call symm_dm_across_atoms(this,this%nks)
    ! Get expansion coefficients only
    call hm_expand_all_herm(this,this%nks,this%nks_coef,this%hm_l,1,.true.)
    ! Symmetrization (note non-zero diagonals)
    call hm_expand_all_sym(this,this%nks,this%hm,.true.,.true.)
    maxerr=maxval(abs(this%nks-copy))
    this%symerr=max(this%symerr,maxerr)
    call output_matrices(ob//'-sym',this%nks,this%na2112,this%num_imp, &
            &this%na2_list,io,1)
    call report_maxerr(maxerr,ob,io)
    if(ob/="nks")then
        return
    endif
    call chk_eigens_matrix_list('nks-sym',this%nks,this%na2112,this%num_imp, &
            &this%na2_list,io,0)
    call calc_nks_tot(this,io)
    if(present(lchk))then
        call chk_nks_tiny(this)
    endif
    return

    end subroutine calc_nks_pp


    subroutine symm_dm_across_atoms(this,dm)
    class(localstore_ob) :: this
#ifdef real_version
    real(q),  &
#else
    complex(q),  &
#endif
            &intent(inout)::dm(this%na2112)

    integer i,j,isum,nbase,na22
#ifdef real_version
    real(q)  &
#else
    complex(q)  &
#endif
            &dms(this%na2max*this%na2max)

    dms=0; isum=0; nbase=0
    do i=1,this%num_imp
        na22=this%na2_list(i)**2
        nbase=nbase+na22
        if(isum==0)then
            dms(1:na22)=dm(nbase-na22+1:nbase)
        else
            dms(1:na22)=dms(1:na22)+dm(nbase-na22+1:nbase)
        endif
        isum=isum+1
        if(i<this%num_imp)then
            if(this%imap_list(i)==this%imap_list(i+1))cycle
        endif
        if(isum>1)then
            dms=dms/isum
            do j=1,isum
                dm(nbase-j*na22+1:nbase-(j-1)*na22)=dms(1:na22)
            enddo
        endif
        isum=0; dms=0
    enddo

    end subroutine symm_dm_across_atoms


    subroutine calc_nks_tot(this,io)
    class(localstore_ob) :: this
    integer,intent(in)::io

    integer i
    real(q) res

    do i=1,this%num_imp
        if(this%ispin_in==1)then
            ! case of lda
            res = trace_a(this%co(i)%nks,this%co(i)%dim2)
            this%co(i)%net=res/2
        else
            ! case of lsda
            call this%co(i)%calc_net()
        endif
    enddo
    if(io>0)then
        write(io,'(" nele_loc total:")')
        write(io,'(4x,2(2f10.5,3x))')(this%co(i)%net,i=1,this%num_imp)
    endif
    return

    end subroutine calc_nks_tot


    subroutine calc_ncvar_pp(this,io)
    class(localstore_ob) :: this
    integer,intent(in)::io

    real(q) maxerr
#ifdef real_version
    real(q)  &
#else
    complex(q)  &
#endif
            &copy(this%na2112)

    call output_matrices('ncv-unsym',this%nc_var,this%na2112,this%num_imp, &
            &this%na2_list,io,1)
    copy=this%nc_var
    ! Get expansion coefficients only
    call hm_expand_all_herm(this,this%nc_var,this%ncv_coef,this%hm_l,1,.true.)
    ! Symmetrization (note non-zero diagonals)
    call hm_expand_all_sym(this,this%nc_var,this%hm,.true.,.true.)
    maxerr=maxval(abs(this%nc_var-copy))
    this%symerr=max(this%symerr,maxerr)
    call output_matrices('ncv-sym',this%nc_var,this%na2112,this%num_imp, &
            &this%na2_list,io,1)
    call report_maxerr(maxerr,'ncv',io)
    return

    end subroutine calc_ncvar_pp


    subroutine chk_nks_tiny(this)
    class(localstore_ob) :: this

    integer i,i1
    real(q) res
    real(q),parameter::small=1.e-6_q

    do i=1,this%num_imp
        do i1=1,this%co(i)%dim2
            res=real(this%co(i)%nks(i1,i1),q)
            if(res<small)then
                write(0,'(" warning in chk_nks_tiny: too small diagonal &
                        &elements",2f16.5)')this%co(i)%nks(i1,i1)
            endif
            this%co(i)%nks(i1,i1)=max(real(this%co(i)%nks(i1,i1),q),small)
        enddo
    enddo
    return

    end subroutine chk_nks_tiny


    subroutine calc_ncphy_pp(this,io)
    class(localstore_ob) :: this
    integer,intent(in)::io

    integer i
    real(q) maxerr
#ifdef real_version
    real(q)  &
#else
    complex(q)  &
#endif
            &copy(this%na2112)

    call output_matrices('ncp-unsym',this%nc_phy,this%na2112,this%num_imp, &
            &this%na2_list,io,1)
    copy=this%nc_phy
    call hm_expand_all_sym(this,this%nc_phy,this%hm,.true.,.true.)
    maxerr=maxval(abs(this%nc_phy-copy))
    this%symerr=max(this%symerr,maxerr)
    call output_matrices('ncp-sym',this%nc_phy,this%na2112,this%num_imp, &
            &this%na2_list,io,1)
    call report_maxerr(maxerr,'ncp',io)
    call renorm_ncphy(this,io)
    call output_matrices('ncp-renorm',this%nc_phy,this%na2112,this%num_imp, &
            &this%na2_list,io,1)
    return

    end subroutine calc_ncphy_pp


    ! for this%h1e, etc.
    subroutine calc_herm_matrices_pp(this,matrices,sname,mb,ltrans,io,mode)
    class(localstore_ob)::this
#ifdef real_version
    real(q),  &
#else
    complex(q),  &
#endif
            &intent(inout)::matrices(this%na2112)
    character,intent(in)::sname*3
    type(matrix_basis),intent(in)::mb
    logical,intent(in)::ltrans
    integer,intent(in)::io,mode

    real(q) maxerr
#ifdef real_version
    real(q)  &
#else
    complex(q)  &
#endif
            & copy(this%na2112)

    call output_matrices(sname//'-unsym',matrices,this%na2112, &
            &this%num_imp,this%na2_list,io,mode)
    call chk_eigens_matrix_list(sname//'-unsym',matrices,this%na2112, &
            &this%num_imp,this%na2_list,io,0)
    copy=matrices
    call symm_dm_across_atoms(this,matrices)
    call hm_expand_all_sym(this,matrices,mb,ltrans,.true.)
    maxerr=maxval(abs(matrices-copy))
    this%symerr=max(this%symerr,maxerr)
    call output_matrices(sname//'-sym',matrices,this%na2112,this%num_imp, &
            &this%na2_list,io,mode)
    call chk_eigens_matrix_list(sname//'-sym',matrices,this%na2112, &
                    &this%num_imp,this%na2_list,io,0)
    call report_maxerr(maxerr,sname,io)
    return

    end subroutine calc_herm_matrices_pp


    subroutine calc_r01_pp(this,io)
    class(localstore_ob)::this
    integer,intent(in)::io

    integer i

    ! r0 -> r
    forall(i=1:this%num_imp)this%co(i)%r=matmul(transpose(this%co(i)%isimix), &
            &this%co(i)%r0)

    call output_matrices('r0-out',this%r0,this%na2112,this%num_imp, &
            &this%na2_list,io,0)
    call calc_r_pp(this,io,'r-out')

    ! r-> z=r^\dagger r
    do i=1,this%num_imp
        call annxb('c',this%co(i)%r,this%co(i)%r,this%co(i)%z, &
                &this%na2_list(i),this%na2_list(i))
    enddo

    call output_matrices('z-out-sym',this%z,this%na2112,this%num_imp, &
            &this%na2_list,io,0)
    call chk_eigens_matrix_list('z',this%z,this%na2112,this%num_imp, &
            &this%na2_list,io,0)
    return

    end subroutine calc_r01_pp


    subroutine calc_r_pp(this,io,sname)
    class(localstore_ob)::this
    integer,intent(in)::io
    character,intent(in)::sname*5

    integer i
    real(q) maxerr
#ifdef real_version
    real(q)  &
#else
    complex(q)  &
#endif
            &copy(this%na2112)

    call output_matrices(sname,this%r,this%na2112,this%num_imp, &
            &this%na2_list,io,0)
    copy=this%r
    call hm_expand_all_general(this,this%r,this%r_coef,this%hm_r,0, &
            &.false.,.false.)
    maxerr=maxval(abs(this%r-copy))
    this%symerr=max(this%symerr,maxerr)
    call output_matrices(sname//'-sym',this%r,this%na2112,this%num_imp, &
            &this%na2_list,io,0)
    call report_maxerr(maxerr,sname,io)
    return

    end subroutine calc_r_pp


    subroutine calc_la1_pp(this,io,sname)
    class(localstore_ob)::this
    integer,intent(in)::io
    character,intent(in)::sname*7

    integer i
    real(q) maxerr
#ifdef real_version
    real(q)  &
#else
    complex(q)  &
#endif
            &copy(this%na2112)

    call output_matrices(sname,this%la1,this%na2112,this%num_imp, &
            &this%na2_list,io,0)
    copy=this%la1
    call hm_expand_all_herm(this,this%la1,this%la1_coef,this%hm_l,0,.false.)
    maxerr=maxval(abs(this%la1-copy))
    this%symerr=max(this%symerr,maxerr) 
    call output_matrices(sname//'-sym',this%la1,this%na2112,this%num_imp, &
            &this%na2_list,io,0)
    call report_maxerr(maxerr,sname,io)
    return

    end subroutine calc_la1_pp


    subroutine calc_lambdac_list(this,io)
    class(localstore_ob)::this
    integer,intent(in) :: io

    integer i

    do i=1,this%num_imp
        call this%co(i)%calc_lambdac()
    enddo
    call hm_expand_all_herm(this,this%la2,this%la2_coef,this%hm_l,0,.false.)
    call output_matrices('la2',this%la2,this%na2112,this%num_imp, &
            &this%na2_list,io,0)
    return

    end subroutine calc_lambdac_list


    subroutine calc_da0_pp(this,io)
    class(localstore_ob)::this
    integer,intent(in)::io

    integer i
    real(q) maxerr
#ifdef real_version
    real(q)  &
#else
    complex(q)  &
#endif
            &copy(this%na2112)

    call output_matrices('d0-unsym',this%d0,this%na2112,this%num_imp, &
            &this%na2_list,io,0)
    copy=this%d0
    call symm_dm_across_atoms(this,this%d0)
    call hm_expand_all_general(this,this%d0,this%d0_coef,this%hm_r,0, &
            &.false.,.false.)
    maxerr=maxval(abs(this%d0-copy))
    this%symerr=max(this%symerr,maxerr)
    call output_matrices('d0-sym',this%d0,this%na2112,this%num_imp, &
            &this%na2_list,io,0)
    call report_maxerr(maxerr,'d0',io)
    return

    end subroutine calc_da0_pp


    subroutine calc_da(this,io)
    class(localstore_ob)::this
    integer,intent(in)::io

    call calc_da0_pp(this,io)
    call d0_to_d(this)
    call hm_expand_all_general(this,this%d,this%d_coef,this%hm_r,0, &
            &.false.,.false.)
    call output_matrices('d-sym',this%d,this%na2112,this%num_imp, &
         &this%na2_list,io,0)
    return

    end subroutine calc_da


    subroutine report_maxerr(err,dname,io)
    integer,intent(in)::io
    character(*),intent(in)::dname
    real(q),intent(in)::err

    if(io<0)return
    write(io,'(" max error due to symm in ",a," = ",f12.6)')dname,err
    return

    end subroutine report_maxerr


    !*************************************************************************
    !< mode>0: get c; mode<0: get a; mode=0: symmetrization (get c then a)
    !*************************************************************************
    subroutine hm_expand_all_herm(this,a,c,mb,mode,ltrans)
    class(localstore_ob) :: this
    type(matrix_basis),intent(in)::mb
#ifdef real_version
    real(q),  &
#else
    complex(q),  &
#endif
            &target,intent(inout)::a(this%na2112)
    real(q),target,intent(inout)::c(mb%dimhst)
    integer,intent(in)::mode
    logical,intent(in)::ltrans

#ifdef real_version
    real(q),  &
#else
    complex(q),  &
#endif
            &pointer :: c_buf(:)

#ifdef real_version
    c_buf=>c
#else
    allocate(c_buf(mb%dimhst))
    c_buf=c
#endif
    call hm_expand_all_general(this,a,c_buf,mb,mode,ltrans,.true.)
#ifndef real_version
    c=c_buf
    deallocate(c_buf)
#endif
    nullify(c_buf)
    return

    end subroutine hm_expand_all_herm


    subroutine hm_expand_all_general(this,a,c,mb,mode,ltrans,lherm)
    class(localstore_ob)::this
    type(matrix_basis),intent(in)::mb
#ifdef real_version
    real(q),  &
#else
    complex(q),  &
#endif
            &target,intent(inout)::a(this%na2112)
#ifdef real_version
    real(q),  &
#else
    complex(q),  &
#endif
            &intent(inout)::c(mb%dimhst)
    integer,intent(in)::mode
    logical,intent(in)::ltrans,lherm

    integer i,imap,ibase,abase,hsbase,nbase,dim_hs,na2,na22
#ifdef real_version
    real(q),  &
#else
    complex(q),  &
#endif
            &pointer::p_a(:,:),p_hs(:,:,:)

    ibase=0
    abase=0
    hsbase=0
    do i=1,this%num_imp
        na2=this%na2_list(i)
        na22=na2**2
        dim_hs=mb%dim_hs(i)
        if(dim_hs>0)then
            p_a(1:na2,1:na2)=>a(abase+1:abase+na22)
            p_hs(1:na2,1:na2,1:dim_hs)=>mb%hs(hsbase+1:hsbase+na22*dim_hs)
            imap=this%imap_list(i)
            if(imap==i)then
                call get_hm_expand(p_a,p_hs,na2,dim_hs, &
                        &c(ibase+1:ibase+dim_hs),mode,ltrans,lherm)
                ibase=ibase+dim_hs
            elseif(mode<=0)then
                nbase=sum(this%na2_list(1:imap-1)**2)
                a(abase+1:abase+na22)=a(nbase+1:nbase+na22)
            endif
        elseif(mode<0)then
            a(abase+1:abase+na22)=0
        endif
        abase=abase+na22
        hsbase=hsbase+na22*dim_hs
    enddo
    return

    end subroutine hm_expand_all_general


    subroutine hm_expand_all_sym(this,a,mb,ltrans,lherm)
    class(localstore_ob) :: this
#ifdef real_version
    real(q),  &
#else
    complex(q),  &
#endif
            &target,intent(inout)::a(this%na2112)
    type(matrix_basis),intent(in)::mb
    logical,intent(in)::ltrans,lherm

    integer i,abase,hsbase,dim_hs,na2,na22
#ifdef real_version
    real(q),  &
#else
    complex(q),  &
#endif
            &pointer::p_a(:,:),p_hs(:,:,:)

    abase=0
    hsbase=0
    do i=1,this%num_imp
        na2=this%na2_list(i)
        na22=na2**2
        dim_hs=mb%dim_hs(i)
        if(dim_hs>0)then
            p_a(1:na2,1:na2)=>a(abase+1:abase+na22)
            p_hs(1:na2,1:na2,1:dim_hs)=>mb%hs(hsbase+1:hsbase+na22*dim_hs)
            call get_hm_expand(p_a,p_hs,na2,dim_hs,ltrans,lherm)
        endif
        abase=abase+na22
        hsbase=hsbase+na22*dim_hs
    enddo
    return

    end subroutine hm_expand_all_sym


    subroutine site_wise_rtrans(this,a,n,mode)
    ! for qp hamiltonian, always complex
    class(localstore_ob) :: this
    integer,intent(in)::n,mode
    complex(q),intent(inout)::a(n,n)

    integer nbase,i,naso
    complex(q),pointer::p_trans(:,:)
    complex(q),target::utrans(this%nasomax**2)

    nbase=0
    do i=1,this%num_imp
        naso=this%co(i)%dimso
        p_trans(1:naso,1:naso) => utrans(1:naso**2)
        if(this%iso==1)then
            p_trans=this%co(i)%db2sab(1:naso,1::2)
        else
            p_trans = this%co(i)%db2sab
        endif
        if(mode<0)then
            p_trans=conjg(p_trans)
        endif
        call anmxbmm("n",a(:,nbase+1:nbase+naso), &
                &p_trans,n,naso)
        nbase=nbase+naso
    enddo
    nullify(p_trans)
    return

    end subroutine site_wise_rtrans


    subroutine map_loc_bnd_matrix(this,x,y,lback)
    class(localstore_ob) :: this
#ifdef real_version
    real(q),  &
#else
    complex(q),  &
#endif
            &intent(inout)::x(this%na2112)
    ! band array is always complex
    complex(q),intent(inout)::y(this%nasotot,this%nasotot,this%nspin)
    logical,intent(in)::lback

    integer i,isp
    integer na2,naso,nbase,ibase
    complex(q),target::buf(this%na2max**2)
    complex(q),pointer::p_buf(:,:)

    nbase=0
    ibase=0
    do i=1,this%num_imp
        na2=this%co(i)%dim2
        naso=this%co(i)%dimso
        p_buf(1:na2,1:na2)=>buf(1:na2**2)
        if(.not.lback)then
            buf(1:na2**2)=x(ibase+1:ibase+na2*na2)
            if(this%iso==1)call orbital_spin_trans(p_buf,na2,.false.)
            do isp=1,this%nspin
                y(nbase+1:nbase+naso,nbase+1:nbase+naso,isp)= &
                        &p_buf((isp-1)*naso+1:isp*naso,(isp-1)*naso+1:isp*naso)
            enddo
        else
            buf=0
            do isp=1,this%nspin
                p_buf((isp-1)*naso+1:isp*naso,(isp-1)*naso+1:isp*naso)= &
                        &y(nbase+1:nbase+naso,nbase+1:nbase+naso,isp)
            enddo
            if(this%ispo==1)p_buf(1+naso:na2,1+naso:na2)=p_buf(1:naso,1:naso)
            if(this%iso==1)call orbital_spin_trans(p_buf,na2,.true.)
            x(ibase+1:ibase+na2*na2)=buf(1:na2**2)
        endif
        nbase=nbase+naso
        ibase=ibase+na2*na2
    enddo
    return

    end subroutine map_loc_bnd_matrix


    !*************************************************************************
    ! renormalize nc_phy
    !*************************************************************************
    subroutine renorm_ncphy(this,io)
    class(localstore_ob) :: this
    integer,intent(in)::io

    integer i
    real(q) sum_ks,sum_phy,sum_var

    do i=1,this%num_imp
        sum_ks =trace_a(this%co(i)%nks   , this%co(i)%dim2)
        sum_phy=trace_a(this%co(i)%nc_phy, this%co(i)%dim2)
        sum_var=trace_a(this%co(i)%nc_var, this%co(i)%dim2)
        if(io>0)then
            write(io,'(" imp=",i3," sum_phy-sum_var=",f16.8)') &
                    &i,sum_phy-sum_var
            write(io,'(7x       ," sum_phy-sum_ks =",f16.8, &
                    &" would be renormalized!")')sum_phy-sum_ks
        endif
        if(sum_phy<1.e-16_q)then
            this%co(i)%nc_phy = 0._q
        else
            this%co(i)%nc_phy=this%co(i)%nc_phy/sum_phy*sum_ks ! renormalized
        endif
    enddo
    return

    end subroutine renorm_ncphy


    !*************************************************************************
    subroutine d0_to_d(this)
    class(localstore_ob) :: this
    integer i

    do i=1,this%num_imp
        call this%co(i)%d0_to_d()
    enddo
    return

    end subroutine d0_to_d


    subroutine gh5_create_impurity_groups(this,gh5,mode,i_f)
    class(localstore_ob) :: this
    class(hdf5_ob)::gh5
    integer,intent(in)::mode,i_f

    integer i

    do i=1,this%num_imp
        if(mode<=0.or.this%imap_list(i)==i)then
            call gh5%gcreate("/impurity_"//trim(int_to_str(i-1)),i_f)
        endif
    enddo
    return

    end subroutine gh5_create_impurity_groups


    subroutine h5prepare_hembed_list(this,gh5,mode)
    class(localstore_ob) :: this
    class(hdf5_ob)::gh5
    integer,intent(in)::mode

    logical lexist

    inquire(file="HEmbed.h5",exist=lexist)
    if(lexist)then
        call h5update_hembed_list(this,gh5,mode)
    else
        call h5write_hembed_list(this,gh5,mode)
    endif
    return

    end subroutine h5prepare_hembed_list


    subroutine h5write_hembed_list(this,gh5,mode)
    class(localstore_ob) :: this
    class(hdf5_ob)::gh5
    integer,intent(in)::mode

    integer i

    call gh5%fopen('HEmbed.h5',1,"w")
    call gh5_create_impurity_groups(this,gh5,1,1)
    do i=1,this%num_imp
        if(this%imap_list(i)==i)then
            call h5write_hembed(this,gh5,i,this%na2_list(i),1,mode)
        endif
    enddo
    call h5write_ndlc_coefs(this,gh5,1)
    call gh5%fclose(1)
    return

    end subroutine h5write_hembed_list


    subroutine h5write_hembed(this,gh5,i,na2,i_f,mode)
    class(localstore_ob) :: this
    class(hdf5_ob)::gh5
    integer,intent(in)::i,na2,i_f,mode

#ifdef real_version
    real(q)  &
#else
    complex(q)  &
#endif
            &h1e2(na2,na2)
    character(39) root

    h1e2=this%co(i)%h1e
    ! adding (orbital dependent) DC term to h1e
    if(associated(this%co(i)%vext).and.this%ivext/=0)then
        h1e2(:na2,:na2)=h1e2(:na2,:na2)+this%co(i)%vext
    endif
    root="/impurity_"//int_to_str(i-1)

    call gh5%dwrite(h1e2,na2,na2,trim(root)//'/H1E',i_f)
    call gh5%dwrite(this%co(i)%v2e,na2,na2,na2,na2,trim(root)//'/V2E',i_f)
    if(mode/=10)then
        ! not hartree-fock calculation.
        call gh5%dwrite(this%co(i)%d,na2,na2,trim(root)//'/D',i_f)
        call gh5%dwrite(this%co(i)%la2,na2,na2,trim(root)//'/LAMBDA',i_f)
    endif
    call gh5%dwrite(this%co(i)%m_struct,na2,na2,trim(root)//'/m_struct',i_f)
    call gh5%gopen(trim(root)//"/",i_f,1)
    call gh5%awrite(na2,'na2',1)
    call gh5%awrite(this%nval_bot_list(i),'nval_bot',1)
    call gh5%awrite(this%nval_top_list(i),'nval_top',1)
    call gh5%awrite(this%mott(i)%nsorb,'norb_mott',1)
    if(this%mott(i)%nsorb>0)then
        call gh5%awrite(this%mott(i)%nelect,'nelect_mott',1)
        call gh5%awrite(this%mott(i)%idx_orb,this%mott(i)%nsorb,'iorb_mott',1)
    endif
    call gh5%gclose(1)
    call gh5%dwrite(this%co(i)%nks,na2,na2,trim(root)//'/NKS',i_f)
    return

    end subroutine h5write_hembed


    subroutine h5write_ndlc_coefs(this,gh5,i_f)
    class(localstore_ob) :: this
    class(hdf5_ob)::gh5
    integer,intent(in)::i_f

    if(this%hm_l%dimhst>0)then
        call gh5%dwrite(this%la2_coef,this%hm_l%dimhst,"/lc_coef",i_f)
        call gh5%dwrite(this%nks_coef,this%hm_l%dimhst,"/nks_coef",i_f)
    endif
    if(this%hm_r%dimhst>0)then
        call gh5%dwrite(this%d_coef,this%hm_r%dimhst,"/d_coef",i_f)
        call gh5%dwrite(this%d0_coef,this%hm_r%dimhst,"/d0_coef",i_f)
    endif
    return

    end subroutine h5write_ndlc_coefs


    subroutine h5update_hembed_list(this,gh5,mode)
    class(localstore_ob) :: this
    class(hdf5_ob)::gh5
    integer,intent(in)::mode

    integer i

    call gh5%fopen('HEmbed.h5',1,"rw")
    do i=1,this%num_imp
        if(this%imap_list(i)==i)then
            call this%co(i)%update_hembed(gh5,"/impurity_"//&
                    &trim(int_to_str(i-1)),1,mode)
        endif
    enddo
    call h5write_ndlc_coefs(this,gh5,1)
    call gh5%fclose(1)
    return

    end subroutine h5update_hembed_list


    subroutine h5read_hembed_list(this,gh5,io,mode)
    class(localstore_ob) :: this
    class(hdf5_ob)::gh5
    integer,intent(in) :: io,mode

    integer i,imap

    call gh5%fopen('HEmbed.h5',1,"r")
    do i=1,this%num_imp
        imap=this%imap_list(i)
        if(imap==i)then
            call this%co(i)%read_hembed(gh5,"/impurity_"// &
                    &trim(int_to_str(i-1)),io,1,mode)
        else
            if(mode/=10)then
                this%co(i)%la2=this%co(imap)%la2
                this%co(i)%d=this%co(imap)%d
                this%co(i)%nc_phy=this%co(imap)%nc_phy
                this%co(i)%nc_var=this%co(imap)%nc_var
                this%co(i)%r0=this%co(imap)%r0
            endif
            this%co(i)%nks=this%co(imap)%nks
            this%co(i)%h1e=this%co(imap)%h1e
            this%co(i)%eu2=this%co(imap)%eu2
        endif
    enddo
    call gh5%fclose(1)
    return

    end subroutine h5read_hembed_list


    subroutine calc_et1_list(this)
    class(localstore_ob) :: this

    integer i

    do i=1,this%num_imp
        this%co(i)%et1=sum(this%co(i)%h1e*this%co(i)%nc_phy)
    enddo
    return

    end subroutine calc_et1_list


    subroutine symm_check(this,io)
    class(localstore_ob)::this
    integer,intent(in)::io


    if(io<0)return
    write(io,'(" maximal local symmetrization error = ", f10.5)')this%symerr
    if(this%symerr>0.05_q)then
        write(io,'(&
                &" WARNING: LARGE LOCAL SYMM ERROR!",/, &
                &"   YOU ARE USING TOO HIGH SYMMETRY FOR THIS SYSTEM!")')
    endif
    return

    end subroutine symm_check



end module localstore
