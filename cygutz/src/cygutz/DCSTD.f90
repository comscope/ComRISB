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

module dcstd
    use gprec, only: q
    use gutil, only: uhau
    use gconstant
    use localstore
    use ghdf5
    implicit none
    private

    type,public::dcstd_ob
        integer::mode=1
        real(q),allocatable::e(:),vpot(:,:),u_avg(:),j_avg(:),nelf(:,:)

        contains
        procedure::init=>init_dc_std
        procedure::calc_vlist=>calc_vdc_list
        procedure::calc_edc_list=>calc_edc_list
        procedure::write_results=>h5write_results
        procedure::add_to_la1_list=>add_vdc_to_la1_list
        procedure::chk_err=>chk_nelf_err
    end type dcstd_ob

    contains


    subroutine init_dc_std(this,loc,gh5,io)
    class(dcstd_ob)::this
    class(localstore_ob)::loc
    class(hdf5_ob)::gh5
    integer,intent(in)::io

    allocate(this%vpot(2,loc%num_imp), &
            &this%u_avg(loc%num_imp), &
            &this%j_avg(loc%num_imp), &
            &this%nelf(2,loc%num_imp), &
            &this%e(loc%num_imp))
    this%u_avg=0; this%j_avg=0; this%nelf=0
    call set_nelf_list(this,loc,gh5,io)
    return

    end subroutine init_dc_std


    subroutine set_nelf_list(this,loc,gh5,io)
    class(dcstd_ob)::this
    class(localstore_ob)::loc
    class(hdf5_ob)::gh5
    integer,intent(in)::io

    integer i,err
    logical lexist

    call gh5%fopen('GParam.h5',1,"rw")
    call gh5%gopen("/",1,1)
    call gh5%read(this%mode,'dc_mode',1)
    if(this%mode>0)then
        call gh5%read(this%u_avg(1),'dc_u_avg_list',1)
        call gh5%read(this%j_avg(1),'dc_j_avg_list',1)
    endif
    if(this%mode>1)then
        call gh5%exists(1,'dc_nelf_list',lexist)
        if(lexist)then
            ! pass address only.
            call gh5%read(this%nelf(1,1),'dc_nelf_list',1)
        else
            do i=1,loc%num_imp
                this%nelf(:,i)=loc%co(i)%net
            enddo
            call gh5%awrite(this%nelf,2,loc%num_imp,"dc_nelf_list",1)
        endif
    endif
    call gh5%gclose(1)
    call gh5%fclose(1)

    if(io>0)then
        select case(this%mode)
        case(0)
            write(io,'(" no double counting.")')
        case(1)
            write(io,'(" standard fully localized limit (FLL) dc.")')
        case(2)
            write(io,'(" fixed double counting.")')
        case(12)
            write(io,'(" standard FLL dc with dc only updated at the &
                    &electron density cycles.")')
        case default
            stop ' error: dc mode not defined!'
        end select

        if(this%mode>1)then
            write(io,'(" input nelf:")')
            write(io,'(8x,5(f8.3,2x))')(sum(this%nelf(:,i)), i=1,loc%num_imp)
        endif

        if(this%mode>0)then
            write(io,'(" average hubbard u list:")')
            write(io,'(8x,5(f8.3,2x))')this%u_avg
            write(io,'(" average hund j list:")')
            write(io,'(8x,5(f8.3,2x))')this%j_avg
        endif
    endif
    return

    end subroutine set_nelf_list


    subroutine chk_nelf_err(this,loc,io)
    class(dcstd_ob)::this
    class(localstore_ob)::loc
    integer,intent(in)::io

    integer idx(1),i
    real(q) ndiff(loc%num_imp)

    ! no double counting.
    if(this%mode==0)return 

    do i=1,loc%num_imp
        ndiff(i)=sum(loc%co(i)%net-this%nelf(:,i))
    enddo
    idx=maxloc(abs(ndiff))

    if(io>0)then
        write(io,'(" max nelf diff=",f14.8)')ndiff(idx(1))
    endif
    return

    end subroutine chk_nelf_err


    subroutine calc_vdc_list(this,loc)
    class(dcstd_ob)::this
    class(localstore_ob)::loc

    integer i,j,na2,na
    real(q) ntot,ntotp
    complex(q),target::zbuf(loc%na2max**2)
    complex(q),pointer::p_la2(:,:)


    ntot=1._q; ntotp=1._q
    do i=1,loc%num_imp
        if(this%mode/=1)then
            ntot=sum(this%nelf(:,i))
            ntotp=sum(loc%co(i)%net)
        endif
        this%nelf(:,i)=loc%co(i)%net*ntot/ntotp
    enddo
    do i=1,loc%num_imp
        ntot=sum(this%nelf(:,i))
        do j=1,2
            this%vpot(j,i)=this%u_avg(i)*(ntot-.5_q)- &
                    &this%j_avg(i)*(this%nelf(j,i)-0.5_q)
        enddo
    enddo
    ! add to co%la2
    loc%la2=0
    do i=1,loc%num_imp
        na2=loc%na2_list(i); na=na2/2
        if(.not.associated(loc%db2sab))then
            do j=1,na
                ! spin-index-faster convention
                loc%co(i)%la2(2*j-1,2*j-1)=-this%vpot(1,i)
                loc%co(i)%la2(2*j,2*j)=-this%vpot(2,i)
            enddo
        else
            do j=1,na
                ! orbital-index-faster convention
                loc%co(i)%la2(j,j)=-this%vpot(1,i)
                loc%co(i)%la2(j+na,j+na)=-this%vpot(2,i)
            enddo
            p_la2(1:na2,1:na2)=>zbuf(1:na2*na2)
            p_la2=loc%co(i)%la2
            ! subsequent transformation to the symmetry-adapted basis
            call uhau(p_la2,loc%co(i)%db2sab,na2,na2)
            loc%co(i)%la2=p_la2
            if(maxval(abs(loc%co(i)%la2-p_la2))>1.e-6_q)then
                stop "error: real version with complex lambdac!"
            endif
        endif
    enddo
    return

    end subroutine calc_vdc_list


    subroutine add_vdc_to_la1_list(this,loc)
    class(dcstd_ob)::this
    class(localstore_ob)::loc

    call calc_vdc_list(this,loc)
    loc%la1=loc%la1+loc%la2
    return

    end subroutine add_vdc_to_la1_list


    subroutine calc_edc_list(this,loc)
    class(dcstd_ob)::this
    class(localstore_ob)::loc

    integer i,isp
    real(q) ntot

    do i=1,loc%num_imp
        this%e(i)=0
        if(this%mode==2)then
            ! nominal dc
            do isp=1,2
                this%e(i)=this%e(i)+this%vpot(isp,i)*loc%co(i)%net(isp)
            enddo
        else
            ! full-localized limit dc
            ntot=sum(loc%co(i)%net)
            this%e(i)=this%e(i)+this%u_avg(i)*ntot*(ntot-1)/2
            do isp=1,2
                this%e(i)=this%e(i)-this%j_avg(i)*loc%co(i)%net(isp)* &
                        &(loc%co(i)%net(isp)-1)/2
            enddo
        endif
    enddo
    return
    
    end subroutine calc_edc_list


    subroutine h5write_results(this,loc,gh5,i_g)
    class(dcstd_ob)::this
    class(hdf5_ob)::gh5
    class(localstore_ob)::loc
    integer,intent(in)::i_g

    integer i
    real(q)::nelf_list(2,loc%num_imp)

    call gh5%awrite(this%e,loc%num_imp,"interaction dc energies",i_g)
    call gh5%awrite(this%nelf,2,loc%num_imp,'dc_nelf_list_inp',i_g)
    do i=1,loc%num_imp
        nelf_list(:,i)=loc%co(i)%net
    enddo
    call gh5%awrite(nelf_list,2,loc%num_imp,'dc_nelf_list_out',i_g)
    return

    end subroutine h5write_results



end module dcstd
