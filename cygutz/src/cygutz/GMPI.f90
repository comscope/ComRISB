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

module gmpi
    use gprec
    implicit none
    include "mpif.h"
    private

    type,public::mpi_ob
        integer::comm=0
        integer::master=0
        integer::nprocs=1
        integer::myrank=0
        integer::io=-1

        contains
        procedure::init=>init_gmpi
        procedure::finalize=>finalize_gmpi
        procedure::barrier=>gmpi_barrier
        procedure,private::imax1_master=>imax1_master_mpi
        procedure,private::imax_master=>imax_master_mpi
        procedure,private::dmax1_master=>dmax1_master_mpi
        procedure,private::dmax_master=>dmax_master_mpi
        procedure,private::imin1_master=>imin1_master_mpi
        procedure,private::imin_master=>imin_master_mpi
        procedure,private::isum1_master=>isum1_master_mpi
        procedure,private::isum_master=>isum_master_mpi
        procedure,private::dsum1_master=>dsum1_master_mpi
        procedure,private::dsum_master=>dsum_master_mpi
        procedure,private::zsum1_master=>zsum1_master_mpi
        procedure,private::zsum_master=>zsum_master_mpi
        procedure,private::imax1_all=>imax1_all_mpi
        procedure,private::imax_all=>imax_all_mpi
        procedure,private::dmax1_all=>dmax1_all_mpi
        procedure,private::dmax_all=>dmax_all_mpi
        procedure,private::imin1_all=>imin1_all_mpi
        procedure,private::imin_all=>imin_all_mpi
        procedure,private::isum1_all=>isum1_all_mpi
        procedure,private::isum_all=>isum_all_mpi
        procedure,private::dsum1_all=>dsum1_all_mpi
        procedure,private::dsum_all=>dsum_all_mpi
        procedure,private::zsum1_all=>zsum1_all_mpi
        procedure,private::zsum_all=>zsum_all_mpi
        procedure,private::zbcast=>zbcast_mpi
        procedure,private::zbcast1=>zbcast1_mpi
        procedure,private::zbcast2=>zbcast2_mpi
        procedure,private::dbcast=>dbcast_mpi
        procedure,private::dbcast1=>dbcast1_mpi
        procedure,private::dbcast2=>dbcast2_mpi
        procedure,private::ibcast=>ibcast_mpi
        procedure,private::ibcast1=>ibcast1_mpi
        procedure,private::ibcast2=>ibcast2_mpi
        procedure,private::lbcast1=>lbcast1_mpi
        procedure,private::sbcast=>sbcast_mpi
        generic::max_master=>imax1_master,imax_master,&
                &dmax1_master,dmax_master
        generic::min_master=>imin1_master,imin_master
        generic::sum_master=>isum1_master,isum_master, &
                &dsum1_master,dsum_master,zsum1_master,zsum_master
        generic::max_all=>imax1_all,imax_all,dmax1_all,dmax_all
        generic::min_all=>imin1_all,imin_all
        generic::sum_all=>isum1_all,isum_all,dsum1_all, &
                &dsum_all,zsum1_all,zsum_all
        generic::bcast=>zbcast,zbcast1,zbcast2,dbcast,dbcast1,dbcast2,&
                &ibcast,ibcast1,ibcast2,lbcast1,sbcast
    end type mpi_ob

    contains


    subroutine init_gmpi(this)
    class(mpi_ob)::this

    integer iflag,ierr

    call mpi_initialized(iflag,ierr)
    if(iflag==0)then
        call mpi_init(ierr)
    endif
    this%comm=mpi_comm_world
    call mpi_comm_rank(this%comm,this%myrank,ierr)
    call mpi_comm_size(this%comm,this%nprocs,ierr)
    if(this%myrank==this%master)then
        this%io=6
    endif
    return

    end subroutine init_gmpi


    subroutine finalize_gmpi(this)
    class(mpi_ob)::this
    integer iflag,ierr

    call mpi_finalized(iflag,ierr)
    if(iflag==0)then
        call mpi_finalize(ierr)
    endif
    return

    end subroutine finalize_gmpi


    subroutine gmpi_barrier(this)
    class(mpi_ob)::this

    integer ierr
      
    call mpi_barrier(this%comm,ierr)
      
    end subroutine gmpi_barrier
    

    subroutine imax1_master_mpi(this,n)
    class(mpi_ob)::this
    integer n
      
    integer i(1)
      
    i(1)=n; call imax_master_mpi(this,i,1); n=i(1)
    return
      
    end subroutine imax1_master_mpi
    

    subroutine imax_master_mpi(this,i,n)
    class(mpi_ob)::this
    integer n,i(n)
      
    integer ierr,maxi(n)
      
    call mpi_reduce(i,maxi,n,mpi_integer,mpi_max,this%master,this%comm,ierr)
    if(this%myrank.eq.this%master)i=maxi
    return
      
    end subroutine imax_master_mpi
    

    subroutine dmax1_master_mpi(this,a)
    class(mpi_ob)::this
    real(q) a
      
    real(q) b(1)
      
    b(1)=a; call dmax_master_mpi(this,b,1); a=b(1)
    return
      
    end subroutine dmax1_master_mpi
    

    subroutine dmax_master_mpi(this,a,n)
    class(mpi_ob)::this
    integer n
    real(q) a(n)
      
    integer ierr
    real(q) maxa(n)
      
    call mpi_reduce(a,maxa,n,mpi_double_precision,mpi_max,this%master, &
            &this%comm,ierr)
    if(this%myrank.eq.this%master)a=maxa
    return
      
    end subroutine dmax_master_mpi
   

    subroutine imin1_master_mpi(this,n)
    class(mpi_ob)::this
    integer n
      
    integer m(1)
      
    m(1)=n; call imin_master_mpi(this,m,1); n=m(1)
    return
      
    end subroutine imin1_master_mpi
      

    subroutine imin_master_mpi(this,i,n)
    class(mpi_ob)::this
    integer n,i(n)
      
    integer ierr,mini(n)
      
    call mpi_reduce(i,mini,n,mpi_integer,mpi_min,this%master,this%comm, &
            &ierr)
    if(this%myrank.eq.this%master)i=mini
    return
      
    end subroutine imin_master_mpi
      

    subroutine isum1_master_mpi(this,i)
    class(mpi_ob)::this
    integer i
      
    integer j(1)
      
    j(1)=i; call isum_master_mpi(this,j,1); i=j(1)
    return
      
    end subroutine isum1_master_mpi
    

    subroutine isum_master_mpi(this,i,n)
    class(mpi_ob)::this
    integer n,i(n)
      
    integer ierr
    integer,allocatable::j(:)
      
    allocate(j(n)); j=0
    call mpi_reduce(i,j,n,mpi_integer,mpi_sum,this%master,this%comm,ierr)
    if(this%myrank.eq.this%master)i=j
    deallocate(j)
    return
      
    end subroutine isum_master_mpi
    

    subroutine dsum1_master_mpi(this,a)
    class(mpi_ob)::this
    real(q) a
      
    real(q) b(1)
      
    b(1)=a; call dsum_master_mpi(this,b,1); a=b(1)
    return
      
    end subroutine dsum1_master_mpi
    

    subroutine dsum_master_mpi(this,a,n)
    class(mpi_ob)::this
    integer n
    real(q) a(n)
      
    integer ierr
    real(q),allocatable::b(:)
      
    allocate(b(n)); b=0
    call mpi_reduce(this,a,b,n,mpi_double_precision,mpi_sum,this%master, &
            &this%comm,ierr)
    if(this%myrank.eq.this%master)a=b
    deallocate(b)
    return
      
    end subroutine dsum_master_mpi
    

    subroutine zsum1_master_mpi(this,a)
    class(mpi_ob)::this
    complex(q) a
      
    complex(q) b(1)
      
    b(1)=a; call zsum_master_mpi(this,b,1); a=b(1)
    return
      
    end subroutine zsum1_master_mpi
    

    subroutine zsum_master_mpi(this,a,n)
    class(mpi_ob)::this
    integer n
    complex(q) a(n)
      
    integer ierr
    complex(q),allocatable::b(:)
      
    allocate(b(n)); b=0
    call mpi_reduce(a,b,n,mpi_double_complex,mpi_sum,this%master, &
            &this%comm,ierr)
    if(this%myrank.eq.this%master)a=b
    deallocate(b)
    return
      
    end subroutine zsum_master_mpi
    

    subroutine imax1_all_mpi(this,i)
    class(mpi_ob)::this
    integer i
      
    integer j(1)
      
    j(1)=i; call imax_all_mpi(this,j,1); i=j(1)
    return
      
    end subroutine imax1_all_mpi
    

    subroutine imax_all_mpi(this,i,n)
    class(mpi_ob)::this
    integer n,i(n)
      
    integer ierr,maxi(n)
      
    call mpi_allreduce(i,maxi,n,mpi_integer,mpi_max,this%comm,ierr)
    i=maxi
    return
      
    end subroutine imax_all_mpi
    

    subroutine imin1_all_mpi(this,i)
    class(mpi_ob)::this
    integer i
      
    integer j(1)
      
    j(1)=i; call imin_all_mpi(this,j,1); i=j(1)
    return
      
    end subroutine imin1_all_mpi
    

    subroutine imin_all_mpi(this,i,n)
    class(mpi_ob)::this
    integer n,i(n)
      
    integer ierr,mini(n)
      
    call mpi_allreduce(i,mini,n,mpi_integer,mpi_min,this%comm,ierr)
    i=mini
    return
      
    end subroutine imin_all_mpi
    

    subroutine dmax1_all_mpi(this,a)
    class(mpi_ob)::this
    real(q) a
      
    real(q) b(1)
      
    b(1)=a; call dmax_all_mpi(this,b,1); a=b(1)
    return
      
    end subroutine dmax1_all_mpi
    

    subroutine dmax_all_mpi(this,a,n)
    class(mpi_ob)::this
    integer n
    real(q) a(n)
      
    integer ierr
    real(q) maxa(n)
      
    call mpi_allreduce(a,maxa,n,mpi_double_precision,mpi_max, &
            &this%comm,ierr)
    a=maxa
    return
      
    end subroutine dmax_all_mpi
    

    subroutine isum1_all_mpi(this,i)
    class(mpi_ob)::this
    integer i
      
    integer j(1)
      
    j(1)=i; call isum_all_mpi(this,j,1); i=j(1)
    return
      
    end subroutine isum1_all_mpi
    

    subroutine isum_all_mpi(this,i,n)
    class(mpi_ob)::this
    integer n,i(n)
      
    integer ierr
    integer,allocatable::j(:)
      
    allocate(j(n)); j=0
    call mpi_allreduce(i,j,n,mpi_integer,mpi_sum,this%comm,ierr)
    i=j; deallocate(j)
    return
      
    end subroutine isum_all_mpi
    

    subroutine dsum1_all_mpi(this,a)
    class(mpi_ob)::this
    real(q) a
      
    real(q) b(1)
      
    b(1)=a; call dsum_all_mpi(this,b,1); a=b(1)
    return
      
    end subroutine dsum1_all_mpi
    

    subroutine dsum_all_mpi(this,a,n)
    class(mpi_ob)::this
    integer n
    real(q) a(n)
      
    integer ierr
    real(q),allocatable::b(:)
      
    allocate(b(n)); b=0
    call mpi_allreduce(a,b,n,mpi_double_precision,mpi_sum,this%comm,ierr)
    a=b; deallocate(b)
    return
      
    end subroutine dsum_all_mpi
    

    subroutine zsum1_all_mpi(this,a)
    class(mpi_ob)::this
    complex(q) a
      
    complex(q) b(1)
      
    b(1)=a; call zsum_all_mpi(this,b,1); a=b(1)
    return
      
    end subroutine zsum1_all_mpi
    

    subroutine zsum_all_mpi(this,a,n)
    class(mpi_ob)::this
    integer n
    complex(q) a(n)
      
    integer ierr
    complex(q),allocatable::b(:)
      
    allocate(b(n)); b=0
    call mpi_allreduce(a,b,n,mpi_double_complex,mpi_sum,this%comm,ierr)
    a=b; deallocate(b)
    return
      
    end subroutine zsum_all_mpi


    subroutine dbcast_mpi(this,a,n)
    class(mpi_ob)::this
    integer,intent(in)::n
    real(q),intent(inout)::a(n)

    integer ierr

    call mpi_bcast(a, n, mpi_double_precision, this%master, &
            &mpi_comm_world, ierr)
    return

    end subroutine dbcast_mpi


    subroutine dbcast1_mpi(this,a)
    class(mpi_ob)::this
    real(q),intent(inout)::a

    integer ierr

    call mpi_bcast(a, 1, mpi_double_precision, this%master, &
            &mpi_comm_world, ierr)
    return

    end subroutine dbcast1_mpi


    !! useful to bcast multidimensional array
    !! by passing the first element(address).
    subroutine dbcast2_mpi(this,a,n)
    class(mpi_ob)::this
    integer,intent(in)::n
    real(q),intent(inout)::a

    integer ierr

    call mpi_bcast(a, n, mpi_double_precision, this%master, &
            &mpi_comm_world, ierr)
    return

    end subroutine dbcast2_mpi


    subroutine zbcast_mpi(this,a,n)
    class(mpi_ob)::this
    integer,intent(in)::n
    complex(q),intent(inout)::a(n)

    integer ierr

    call mpi_bcast(a, n, mpi_double_complex, this%master, &
            &mpi_comm_world, ierr)
    return

    end subroutine zbcast_mpi


    subroutine zbcast1_mpi(this,a)
    class(mpi_ob)::this
    complex(q),intent(inout)::a

    integer ierr

    call mpi_bcast(a, 1, mpi_double_complex, this%master, &
            &mpi_comm_world, ierr)
    return

    end subroutine zbcast1_mpi


    subroutine zbcast2_mpi(this,a,n)
    class(mpi_ob)::this
    integer,intent(in)::n
    complex(q),intent(inout)::a

    integer ierr

    call mpi_bcast(a, n, mpi_double_complex, this%master, &
            &mpi_comm_world, ierr)
    return

    end subroutine zbcast2_mpi


    subroutine ibcast_mpi(this,a,n)
    class(mpi_ob)::this
    integer,intent(in)::n
    integer,intent(inout)::a(n)
    
    integer ierr

    call mpi_bcast(a, n, mpi_integer, this%master, &
            &mpi_comm_world, ierr)
    return

    end subroutine ibcast_mpi


    subroutine sbcast_mpi(this,a,n)
    class(mpi_ob)::this
    integer,intent(in)::n
    character(n),intent(inout)::a

    integer ierr

    call mpi_bcast(a, n, mpi_character, this%master, &
            &mpi_comm_world, ierr)
    return

    end subroutine sbcast_mpi


    subroutine ibcast1_mpi(this,a)
    class(mpi_ob)::this
    integer,intent(inout)::a

    integer ierr

    call mpi_bcast(a, 1, mpi_integer, this%master, &
            &mpi_comm_world, ierr)
    return

    end subroutine ibcast1_mpi


    subroutine ibcast2_mpi(this,a,n)
    class(mpi_ob)::this
    integer,intent(in)::n
    integer,intent(inout)::a
    
    integer ierr

    call mpi_bcast(a, n, mpi_integer, this%master, &
            &mpi_comm_world, ierr)
    return

    end subroutine ibcast2_mpi


    subroutine lbcast1_mpi(this,a)
    class(mpi_ob)::this
    logical,intent(inout)::a

    integer ierr

    call mpi_bcast(a, 1, mpi_logical, this%master, &
            &mpi_comm_world, ierr)
    return

    end subroutine lbcast1_mpi

      
end module gmpi
