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

module gkerder
    use gprec
    use ghdf5
    use localstder
    use gkernel
    implicit none
    private
    
    type,extends(gkernel_ob),public::gkerder_ob
        contains
        procedure::save_jacobian=>h5save_jacobian
    end type gkerder_ob

    contains


    subroutine h5save_jacobian(this,gh5,loc,io)
    class(gkerder_ob)::this
    class(hdf5_ob)::gh5
    class(localstder_ob)::loc
    integer,intent(in)::io

    integer i
    real(q) pfpx(this%nx,this%nx)

    if(this%iembeddiag==10)then
        return
    endif
    pfpx=0
    do i=1,this%nx1
        pfpx(1:this%nx1,i)=transfer(loc%pr_coefpr(:,i),pfpx(1:this%nx1,i))
        pfpx(i,i)=pfpx(i,i)-1
        pfpx(1+this%nx1:,i)=loc%pdm_coefpr(:,i)-loc%pncv_coefpr(:,i)
    enddo
    do i=1,this%nx2
        pfpx(1:this%nx1,this%nx1+i)=transfer(loc%pr_coefpl(:,i),  &
                &pfpx(1:this%nx1,this%nx1+i))
        pfpx(1+this%nx1:,this%nx1+i)=loc%pdm_coefpl(:,i)-loc%pncv_coefpl(:,i)
    enddo
    call gh5%fopen("GIter.h5",1,"rw")
    ! fortran to c order
    pfpx=transpose(pfpx)
    call gh5%dwrite(pfpx,this%nx,this%nx,"/pfpx",1)
    call gh5%fclose(1)
    return

    end subroutine h5save_jacobian



end module gkerder
