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

module gkernel
    use gprec
    use ghdf5
    use gutil, only: output_matrices,get_kwarg
    use bandstru
    use localstore
    use gmpi
    use dcstd
    use gtime
    implicit none
    private
    
    type,public::gkernel_ob
        integer::iembeddiag,iupdaterho=0,iwrite=0,ideri=0
        integer::nx,nx1,nx2,iter=0
        real(q)::etot,r_default=0.99_q
        real(q),pointer::x(:)

        contains
        procedure::set_params=>set_params
        procedure::init_x=>init_newton_x
        procedure::dump_x=>h5write_gk_x
        procedure::load_x=>h5read_gk_x
        procedure::calc_band
        procedure::save_verr=>h5save_verr
        procedure::calc_total_energy
        procedure::save_results=>h5save_results
    end type gkernel_ob

    contains


    subroutine set_params(this,gh5)
    class(gkernel_ob)::this
    class(hdf5_ob)::gh5

    character(255) cmd

    call gh5%fopen("GParam.h5",1,"r")
    call gh5%gopen("/",1,1)
    call gh5%read(this%iembeddiag,'giembeddiag',1)
    call gh5%gclose(1)
    call gh5%fclose(1)

    if(this%iembeddiag==10)then
        this%r_default=1._q
    endif
    call get_command(cmd)
    call get_kwarg(cmd,'-r',this%iupdaterho)
    call get_kwarg(cmd,'-w',this%iwrite)
    call get_kwarg(cmd,'-d',this%ideri)
    return

    end subroutine set_params


    subroutine init_newton_x(this,loc)
    class(gkernel_ob)::this
    class(localstore_ob)::loc

    if(this%iembeddiag==10)then
        !! Hartree-Fock or DFT+U
        this%nx1=loc%hm%dimhst
        this%nx=this%nx1
        allocate(this%x(this%nx))
        this%x=loc%nks_coef
    else
        !! initialize solution {R,\lambda} vector
        this%nx1=loc%hm_r%dimhst*loc%r_factor
        this%nx2=loc%hm_l%dimhst
        this%nx=this%nx1+this%nx2
        allocate(this%x(this%nx))
        this%x(1:this%nx1)=transfer(loc%r_coef,this%x(1:this%nx1))
        this%x(1+this%nx1:)=loc%la1_coef
    endif
    return

    end subroutine init_newton_x


    subroutine h5write_gk_x(this,gh5)
    class(gkernel_ob)::this
    class(hdf5_ob)::gh5

    call gh5%fopen("GIter.h5",1,"w")
    call gh5%gopen("/",1,1)
    call gh5%awrite(this%nx,'x_dim',1)
    call gh5%awrite(this%iter,'iter',1)
    call gh5%awrite(this%nx1,'x_dim1',1)
    call gh5%gclose(1)
    if(this%nx>0)then
        call gh5%dwrite(this%x,this%nx,'/x_val',1)
    endif
    call gh5%fclose(1)
    return

    end subroutine h5write_gk_x


    subroutine h5read_gk_x(this,gh5)
    class(gkernel_ob)::this
    class(hdf5_ob)::gh5

    call gh5%fopen("GIter.h5",1,"r")
    call gh5%gopen("/",1,1)
    call gh5%read(this%nx,'x_dim',1)
    call gh5%read(this%nx1,'x_dim1',1)
    call gh5%read(this%iter,'iter',1)
    call gh5%gclose(1)
    this%nx2=this%nx-this%nx1
    allocate(this%x(this%nx))
    if(this%nx>0)then
        call gh5%read(this%x,this%nx,'/x_val',1)
    endif
    call gh5%fclose(1)
    return

    end subroutine h5read_gk_x


    subroutine calc_band(this,bnd,loc,dc,mpi,time)
    class(gkernel_ob)::this
    class(bandstru_ob)::bnd
    class(localstore_ob)::loc
    class(dcstd_ob)::dc
    class(mpi_ob)::mpi
    class(time_ob)::time

    if(this%iembeddiag==10)then
        call calc_band_n(this,bnd,loc,dc,mpi,time)
    else
        call calc_band_rl(this,bnd,loc,dc,mpi,time)
    endif
    return

    end subroutine calc_band


    subroutine calc_band_rl(this,bnd,loc,dc,mpi,time)
    class(gkernel_ob)::this
    class(bandstru_ob)::bnd
    class(localstore_ob)::loc
    class(dcstd_ob)::dc
    class(mpi_ob)::mpi
    class(time_ob)::time

    ! x to {r_coef, la1_coef}
    loc%r_coef=transfer(this%x(1:this%nx1),loc%r_coef)
    loc%la1_coef=this%x(1+this%nx1:)
    call loc%hm_expand_all(loc%r,loc%r_coef,loc%hm_r,-1, &
            &.false.,.false.)
    call loc%hm_expand_all(loc%la1,loc%la1_coef,loc%hm_l,-1,.false.)
    call loc%modify_rl_mott(1,30._q)
    call output_matrices('r-in',loc%r,loc%na2112,loc%num_imp, &
            &loc%na2_list,mpi%io,0)
    call output_matrices('la1-in',loc%la1,loc%na2112,loc%num_imp, &
            &loc%na2_list,mpi%io,0)
    call loc%map_bnd_matrix(loc%r,bnd%r,.false.)
    call loc%map_bnd_matrix(loc%la1,bnd%la1,.false.)
    call bnd%calc_all(loc,time,mpi)
    call bnd%calc_fermi(mpi%io)
    call bnd%modify_mott()
    call loc%modify_rl_mott(1,bnd%ef)
    call bnd%calc_bz_integration_correction()
    call bnd%calc_mup_dn(mpi%io)
    call bnd%calc_nks(loc,mpi)
    call loc%nks_pp(mpi%io)
    call loc%eval_sl_vec_all(1,mpi%io)
    call loc%map_bnd_matrix(loc%nks,bnd%nks,.false.)
    call bnd%calc_da0(loc,mpi)
    call loc%calc_isimix(mpi%io,1)
    call loc%calc_da(mpi%io)
    call dc%calc_vlist(loc)
    call loc%calc_lambdac_list(mpi%io)
    return

    end subroutine calc_band_rl


    subroutine calc_band_n(this,bnd,loc,dc,mpi,time)
    class(gkernel_ob)::this
    class(bandstru_ob)::bnd
    class(localstore_ob)::loc
    class(dcstd_ob) :: dc
    class(mpi_ob) :: mpi
    class(time_ob)::time

    ! x to nks_coef
    loc%nks_coef=this%x
    call loc%hm_expand_all(loc%nks,loc%nks_coef,loc%hm,-1,.false.)
    call loc%calc_nks_tot(mpi%io)
    call output_matrices('nks-in',loc%nks,loc%na2112,loc%num_imp, &
            &loc%na2_list,mpi%io,0)
    call loc%calc_la1_hf_list()
    call dc%add_to_la1_list(loc)
    if(associated(loc%vext).and.loc%ivext>0)then
        loc%la1=loc%la1+loc%vext
    endif
    call loc%la1_pp(mpi%io,'la1-inp')
    call loc%map_bnd_matrix(loc%la1,bnd%la1,.false.)
    call bnd%calc_all(loc,time,mpi)
    call bnd%calc_fermi(mpi%io)
    call bnd%calc_mup_dn(mpi%io)
    call bnd%calc_nks(loc,mpi)
    call loc%nks_pp(mpi%io)
    !! True for HF calculation.
    loc%nc_var=loc%nks
    loc%nc_phy=loc%nks
    call loc%eval_sl_vec_all(1,mpi%io)
    call loc%map_bnd_matrix(loc%nks,bnd%nks,.false.)
    return

    end subroutine calc_band_n


    subroutine h5save_verr(this,gh5,loc,io)
    class(gkernel_ob)::this
    class(hdf5_ob)::gh5
    class(localstore_ob)::loc
    integer,intent(in)::io

    real(q) fvec(this%nx),maxerr,fvec_(max(1,this%nx))

    if(this%iembeddiag==10)then
        fvec=loc%nks_coef-this%x
    else
        fvec(1:this%nx1)=transfer(loc%r_coef,fvec(1:this%nx1)) &
                &-this%x(1:this%nx1)
        fvec(1+this%nx1:)=loc%nks_coef-loc%ncv_coef
    endif

    if(this%nx>0)then
        fvec_=fvec
    else
        fvec_=0
    endif
    if(io>0)then
        write(io,'(" maxerr = ",f12.6)')maxval(abs(fvec_))
    endif
    call gh5%fopen("GIter.h5",1,"rw")
    call gh5%dwrite(fvec_,max(this%nx,1),"/v_err",1)
    call gh5%fclose(1)
    return

    end subroutine h5save_verr


    subroutine calc_total_energy(this,bnd,loc,dc,mpi)
    class(gkernel_ob)::this
    class(bandstru_ob)::bnd
    class(localstore_ob)::loc
    class(dcstd_ob) :: dc
    class(mpi_ob)::mpi

    call loc%calc_et1_list()
    call bnd%calc_ebnd(mpi)
    call loc%calc_edcla1()
    call dc%calc_edc_list(loc)
    loc%egamma_dc=sum(loc%co(:)%eu2)+sum(loc%co(:)%et1)-loc%edcla1-sum(dc%e)
    this%etot=sum(bnd%eband)+loc%egamma_dc
    if(mpi%io>0)then
        write(mpi%io,'(" e_total_model = ",f0.7)')this%etot
    endif
    return

    end subroutine calc_total_energy


    subroutine h5save_results(this,bnd,loc,dc,gh5)
    class(gkernel_ob)::this
    class(hdf5_ob)::gh5
    class(bandstru_ob)::bnd
    class(localstore_ob)::loc
    class(dcstd_ob) :: dc

    call gh5%fopen("GLog.h5",1,"w")
    call gh5%gopen("/",1,1)
    call gh5%awrite(this%etot,"etot_model",1)
    call bnd%write_results(loc,gh5,1)
    call loc%write_results(gh5,1,1)
    call dc%write_results(loc,gh5,1)
    call gh5%gclose(1)
    call gh5%fclose(1)
    return

    end subroutine h5save_results


end module gkernel
