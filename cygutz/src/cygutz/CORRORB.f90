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

module corrorb
    use gprec
    use gutil, only:trace_a,get_hf_pot,dsimix,dpsimix,pfa_pa,uhau,atofa, &
            &orbital_spin_trans,get_coul_exchange
    use ghdf5
    implicit none
    private
      
    type,public :: corrorb_ob
        !< number of correlated orbital, spin-orbital, dim*ispin
        integer dim,dim2,dim4,dimso
        real(q) :: net(2),eu2,et1
        integer :: dim_hs_l,dim_hs_r
#ifdef real_version
        real(q),  &
#else
        complex(q),  &
#endif
                &pointer :: nks(:,:),nc_var(:,:),nc_phy(:,:),isimix(:,:), &
                &r(:,:),r0(:,:),z(:,:),d(:,:),d0(:,:), &
                &la1(:,:),la2(:,:),h1e(:,:),vext(:,:)=>null(), &
                &hs_l(:,:,:),hs_r(:,:,:),v2e(:,:,:,:), &
                &v_j2e(:,:,:,:), &
                &pncvplc(:,:,:)=>null(),pncvpd(:,:,:)=>null(), &
                &pr0plc(:,:,:)=>null(),pr0pd(:,:,:)=>null()
        complex(q),pointer :: sx(:,:),sy(:,:),sz(:,:),  &
                &lx(:,:),ly(:,:),lz(:,:),db2sab(:,:)=>null()
        real(q) :: s_val(3,2),l_val(3,2)  !< one-body operator
        integer,pointer :: m_struct(:,:)=>null()

        contains
        procedure::eval_sl_vec=>eval_co_sl_vec
        procedure::calc_lambdac=>calc_co_lambdac
        procedure::calc_net=>calc_co_net
        procedure::calc_isimix=>calc_co_isimix
        procedure::d0_to_d=>co_d0_to_d
        procedure::nks_patch_order=>co_nks_patch_order
        procedure::calc_la1_hf=>calc_co_la1_hf
        procedure::update_hembed=>h5update_hembed
        procedure::read_hembed=>h5read_hembed
        procedure::set_v_j2e=>set_co_v_j2e
    end type corrorb_ob
      
    contains
    

    subroutine eval_co_sl_vec(this,mode)
    class (corrorb_ob) this
    integer mode
      
    integer i
#ifdef real_version
    real(q)  &
#else
    complex(q)  &
#endif
            &n_(this%dim2,this%dim2)
      
    if(mode==1)then
        n_=this%nks
    else
        n_=this%nc_phy
    endif
    ! Tr(\rho A) = \sum_{ij}{\rho_ij*A_ji} = \sum_{ij}{\rho^{t}_ji*A_ji}
    this%s_val(1,mode)=sum(this%sx*n_)
    this%s_val(2,mode)=sum(this%sy*n_)
    this%s_val(3,mode)=sum(this%sz*n_)
    this%l_val(1,mode)=sum(this%lx*n_)
    this%l_val(2,mode)=sum(this%ly*n_)
    this%l_val(3,mode)=sum(this%lz*n_)
    return
      
    end subroutine eval_co_sl_vec
      

    subroutine calc_co_lambdac(this)
    class (corrorb_ob) this
      
    integer i
    real(q) res
#ifdef real_version
    real(q)  &
#else
    complex(q)  &
#endif
            &rd(this%dim2,this%dim2),pn(this%dim2,this%dim2), &
            &h(this%dim2,this%dim2)

    rd=matmul(this%r,transpose(this%d))
    do i=1,this%dim_hs_l
        ! nks is expanded in terms of h.t
        h=transpose(this%hs_l(:,:,i))
        call pfa_pa(this%nks,pn,h,this%dim2,dsimix,dpsimix) ! p f / p d_n
        res=-real(sum(pn*rd),q)*2
        this%la2=this%la2+this%hs_l(:,:,i)*res
    enddo
    this%la2=this%la2-this%la1
    return
      
    end subroutine calc_co_lambdac


    subroutine calc_co_net(this)
    class (corrorb_ob) this

    integer i
    complex(q) zbuf(this%dim2,this%dim2),sab2db(this%dim2,this%dim2)

    this%net=0
    if(associated(this%db2sab))then
        ! transposition is needed for nks.
        zbuf=transpose(this%nks)
        sab2db=conjg(transpose(this%db2sab))
        call uhau(zbuf,sab2db,this%dim2,this%dim2)
        ! orbital-index faster
        do i=1,this%dim
            this%net(1)=this%net(1)+zbuf(i,i)
            this%net(2)=this%net(2)+zbuf(i+this%dim,i+this%dim)
        enddo
    else
        ! spin-index faster
        do i=1,this%dim
            this%net(1)=this%net(1)+this%nks(2*i-1,2*i-1)
            this%net(2)=this%net(2)+this%nks(2*i,2*i)
        enddo
    endif
    return

    end subroutine calc_co_net

    
    subroutine set_co_v_j2e(this,na2)
    class (corrorb_ob) this
    integer,intent(in) :: na2

    allocate(this%v_j2e(na2,na2,na2,na2))
    call get_coul_exchange(na2,this%v2e,this%v_j2e)
    return

    end subroutine set_co_v_j2e


    !*************************************************************************
    ! isimix = (nks(1-nks))^(-1/2)
    !*************************************************************************
    subroutine calc_co_isimix(this,mode)
    class (corrorb_ob) this
    integer,intent(in)::mode
#ifdef real_version
    real(q)  &
#else
    complex(q)  &
#endif
            &xn(this%dim2,this%dim2)
      
    if(mode==1)then
        xn=this%nks
    else
        xn=this%nc_var
    endif
    xn=xn-matmul(xn,xn)
    call atofa(xn,this%isimix,this%dim2,-12,1.d0,.true.)
    return
      
    end subroutine calc_co_isimix
      
    !*************************************************************************
    ! d = (nks(1-nks))^(-1/2) d0
    !*************************************************************************
    subroutine co_d0_to_d(this)
    class (corrorb_ob) this
      
    this%d=matmul(this%isimix,this%d0)
    return
      
    end subroutine co_d0_to_d


    subroutine co_nks_patch_order(this,ispo,iso)
    class (corrorb_ob) this
    integer,intent(in)::ispo,iso

    if(ispo==1)then
        this%nks(1+this%dimso:,1+this%dimso:)= &
                &this%nks(1:this%dimso,1:this%dimso)
    endif
    if(iso==1)then
        call orbital_spin_trans(this%nks,this%dim2,.true.)
    endif
    return

    end subroutine co_nks_patch_order


    subroutine calc_co_la1_hf(this)
    class (corrorb_ob) this

    this%la1=0
    call get_hf_pot(this%dim2,this%nks,this%v_j2e,this%la1)
    this%la1=this%la1+this%h1e
    return

    end subroutine calc_co_la1_hf


    subroutine calc_co_eu2_hf(this)
    class (corrorb_ob) this
#ifdef real_version
    real(q)  &
#else
    complex(q)  &
#endif
            &vhf(this%dim2,this%dim2)

    vhf=0
    call get_hf_pot(this%dim2,this%nks,this%v_j2e,vhf)
    this%eu2=sum(vhf*this%nks)/2
    return

    end subroutine calc_co_eu2_hf


    subroutine h5update_hembed(this,gh5,root,i_f,mode)
    class(corrorb_ob) :: this
    class(hdf5_ob)::gh5
    character(*),intent(in)::root
    integer,intent(in)::i_f,mode ! f_id index

    if(mode/=10)then
        call gh5%dwrite(this%d,this%dim2,this%dim2,root//'/D',i_f)
        call gh5%dwrite(this%la2,this%dim2,this%dim2,root//'/LAMBDA',i_f)
    endif
    call gh5%dwrite(this%nks,this%dim2,this%dim2,root//'/NKS',i_f)
    return

    end subroutine h5update_hembed


    subroutine h5read_hembed(this,gh5,root,io,i_f,mode)
    class(corrorb_ob) :: this
    class(hdf5_ob)::gh5
    character(*),intent(in)::root
    integer,intent(in)::io,i_f,mode  ! f_id index

    real(q) de,etot
#ifdef real_version
    real(q)  &
#else
    complex(q)  &
#endif
            &h1e4(this%dim4,this%dim4), dm(this%dim4,this%dim4)

    if(mode/=10)then
        ! not hartree-fock calculation.
        call gh5%gopen(root//"/ans/",i_f,1)
        call gh5%read(etot,'emol',1)
        call gh5%gclose(1)
        call gh5%read(dm,this%dim4,this%dim4,root//'/ans/DM',i_f)
        call gh5%read(this%la2,this%dim2,this%dim2,root//'/LAMBDA',i_f)
        call gh5%read(this%d,this%dim2,this%dim2,root//'/D',i_f)
    endif
    call gh5%read(this%nks,this%dim2,this%dim2,root//'/NKS',i_f)
    call gh5%read(this%h1e,this%dim2,this%dim2,root//'/H1E',i_f)

    if(mode/=10)then
        de=trace_a(this%la2,this%dim2)
        if(io>0)then
            write(io,'(" e_mol = ",f0.7," e_tot = ",f0.7)')etot,etot+de
        endif

        ! setup h1e4
        h1e4(1:this%dim2,1:this%dim2)=this%h1e
        h1e4(1+this%dim2:,1:this%dim2)= &
#ifdef real_version
                &this%d
#else
                &conjg(this%d)
#endif
        h1e4(1:this%dim2,1+this%dim2:)=transpose(this%d)
        h1e4(1+this%dim2:,1+this%dim2:)=-this%la2
        this%eu2=etot-sum(h1e4*dm)
        call set_rn_from_embeddm(this,dm)
    else
        call calc_co_eu2_hf(this)
    endif
    return

    end subroutine h5read_hembed


    !> Get the impurity
    !! physical density matrix (ni_{AB} = <c_A^\dagger c_B>),
    !! bath density matrix (nb_{ab} = <f_b f_a^\dagger>
    !!                              = \delta_{ab}-<f_a^\dagger f_b>)
    !! and r0_{aA} = <c_A^\dagger f_a>.
    !! @param dm complete density matrix <{c_A f_b}^\dagger {c_B f_b}>
    !< Complex version.
    subroutine set_rn_from_embeddm(this,dm)
    class(corrorb_ob) :: this
#ifdef real_version
    real(q)  &
#else
    complex(q)  &
#endif
            &:: dm(this%dim4,this%dim4)

    integer i

    !< Get physical density matrix.
    this%nc_phy=dm(:this%dim2,:this%dim2)
    !< Bath auxiliary density matrix
    this%nc_var=-dm(this%dim2+1:,this%dim2+1:)
    forall(i=1:this%dim2) this%nc_var(i,i)=this%nc_var(i,i)+1
    !< r0_{aA} = <c_A^\dagger f_a>, with indices switch.
    this%r0=transpose(dm(:this%dim2,this%dim2+1:))
    return

    end subroutine set_rn_from_embeddm

      
end module corrorb
