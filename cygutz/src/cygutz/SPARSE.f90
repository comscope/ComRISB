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

module sparse
    use gprec
    implicit none
    private

    type,public::sp_matrix
        integer :: nrow=0,ncol=0
        integer,allocatable :: i(:),j(:)
    end type sp_matrix

    type,public,extends(sp_matrix) :: dsp_matrix
        real(q),allocatable :: a(:)
    end type dsp_matrix

    type,public,extends(dsp_matrix) :: dcsr_matrix
    end type dcsr_matrix

    type,public,extends(dsp_matrix) :: dcoo_matrix
        integer :: nnz=0
    end type dcoo_matrix

    type,public,extends(sp_matrix) :: zsp_matrix
        complex(q),allocatable :: a(:)
    end type zsp_matrix

    type,public,extends(zsp_matrix) :: zcsr_matrix
    end type zcsr_matrix

    type,public::zvector
        complex(q),allocatable :: a(:)
    end type zvector

    type,public::ivector
        integer :: imin=0,imax=0
        integer,allocatable :: i(:)
    end type ivector

    public::csr_syamux,sp_amux,alloc_sp_matrix,dealloc_sp_matrix, &
            &coohmux,vh_sp_v,vh_cooh_v,vh_sycsr_v

    interface csr_syamux
        module procedure dcsr_syamuzx_sk, dcsr_sylamuzx_sk, &
                &dcsr_sylamux_sk, zcsr_syamux_sk
    end interface csr_syamux

    interface sp_amux
        module procedure zcsr_amux_sk, zs_dcoomux, ds_coomux
    end interface sp_amux

    interface coohmux
        module procedure zs_dcoohmux, ds_dcoohmux
    end interface coohmux

    interface vh_sp_v
        module procedure zvh_dcoo_zv, dvh_coo_v, zv2h_sdcoo_zv1,  &
                &dv2h_scoo_v1, zvh_dcsr_zv, dvh_csr_v, zvh_csr_v
    end interface vh_sp_v

    interface vh_cooh_v
        module procedure zv2h_sdcooh_zv1, dv2h_scooh_v1
    end interface vh_cooh_v

    interface vh_sycsr_v
        module procedure zvh_sydcsr_zv, dvh_sycsr_v
    end interface vh_sycsr_v

    contains

    subroutine dcsr_syamuzx_sk(a,x,ax)
    type(dcsr_matrix),intent(in)::a
    complex(q),intent(in)::x(*)
    complex(q),intent(inout)::ax(*)
      
    integer i,k,j
      
    do i=1,a%nrow
        do k=a%i(i),a%i(i+1)-1
            j=a%j(k)
            ax(i)=ax(i)+a%a(k)*x(j)
            ! lower triangular part
            if(j==i)cycle
            ax(j)=ax(j)+a%a(k)*x(i)
        enddo
    enddo
    return
      
    end subroutine dcsr_syamuzx_sk


    subroutine dcsr_sylamuzx_sk(coef,a,x,ax)
    real(q),intent(in)::coef
    type(dcsr_matrix),intent(in)::a
    complex(q),intent(in)::x(*)
    complex(q),intent(inout)::ax(*)
      
    integer i,k,j
      
    do i=1,a%nrow
        do k=a%i(i),a%i(i+1)-1
            j=a%j(k)
            ax(i)=ax(i)+a%a(k)*x(j)*coef
            ! lower triangular part
            if(j==i)cycle
            ax(j)=ax(j)+a%a(k)*x(i)*coef
        enddo
    enddo
    return
      
    end subroutine dcsr_sylamuzx_sk


    subroutine dcsr_sylamux_sk(coef,a,x,ax)
    real(q),intent(in)::coef
    type(dcsr_matrix),intent(in)::a
    real(q),intent(in)::x(*)
    real(q),intent(inout)::ax(*)
      
    integer i,k,j
      
    do i=1,a%nrow
        do k=a%i(i),a%i(i+1)-1
            j=a%j(k)
            ax(i)=ax(i)+a%a(k)*x(j)*coef
            ! lower triangular part
            if(j==i)cycle
            ax(j)=ax(j)+a%a(k)*x(i)*coef
        enddo
    enddo
    return
      
    end subroutine dcsr_sylamux_sk


    subroutine zcsr_syamux_sk(a,x,ax)
    type(zcsr_matrix),intent(in)::a
    complex(q),intent(in)::x(*)
    complex(q),intent(inout)::ax(*)
      
    integer i,k,j
      
    do i=1,a%nrow
        do k=a%i(i),a%i(i+1)-1
            j=a%j(k)
            ax(i)=ax(i)+a%a(k)*x(j)
            ! lower triangular part
            if(j==i)cycle
            ax(j)=ax(j)+conjg(a%a(k))*x(i)
        enddo
    enddo
    return
      
    end subroutine zcsr_syamux_sk
      

    subroutine zcsr_amux_sk(a,x,ax)
    type(zcsr_matrix),intent(in)::a
    complex(q),intent(in)::x(*)
    complex(q),intent(inout)::ax(*)
      
    integer i,k,j
      
    do i=1,a%nrow
        do k=a%i(i),a%i(i+1)-1
            j=a%j(k)
            ax(i)=ax(i)+a%a(k)*x(j)
        enddo
    enddo
    return
      
    end subroutine zcsr_amux_sk


    subroutine alloc_sp_matrix(a,nnz,nrow,ncol)
    integer nnz,nrow,ncol
    class(sp_matrix)::a
     
    allocate(a%j(nnz)); a%j=0
    a%nrow=nrow; a%ncol=ncol
    select type (a)
    type is (dcsr_matrix)
        allocate(a%a(nnz), a%i(nrow+1)); a%a=0; a%i=0
    type is (zcsr_matrix)
        allocate(a%a(nnz), a%i(nrow+1)); a%a=0; a%i=0
    type is (dcoo_matrix)
        allocate(a%a(nnz), a%i(nnz)); a%a=0; a%i=0; a%nnz=nnz
    class default
        stop "undefined behavior!"
    end select

    return
      
    end subroutine alloc_sp_matrix


    subroutine dealloc_sp_matrix(a)
    class(sp_matrix)::a

    deallocate(a%j)
    a%nrow=0; a%ncol=0
    select type (a)
    type is (dcsr_matrix)
        deallocate(a%a, a%i)
    type is (zcsr_matrix)
        deallocate(a%a, a%i)
    type is (dcoo_matrix)
        deallocate(a%a, a%i); a%nnz=0
    class default
        stop "undefined behavior!"
    end select
    return

    end subroutine dealloc_sp_matrix

      
    subroutine zs_dcoomux(zs,dcoo,x,y)
    type(dcoo_matrix),intent(in)::dcoo
    complex(q),intent(in)::zs,x(*)
    complex(q),intent(inout)::y(*)
    
    integer inz,i

    do inz=1,dcoo%nnz
        i=dcoo%i(inz)
        y(i)=y(i)+dcoo%a(inz)*x(dcoo%j(inz))*zs
    enddo
    return

    end subroutine zs_dcoomux


    subroutine ds_coomux(zs,dcoo,x,y)
    type(dcoo_matrix),intent(in)::dcoo
    real(q),intent(in)::zs,x(*)
    real(q),intent(inout)::y(*)
    
    integer inz,i

    do inz=1,dcoo%nnz
        i=dcoo%i(inz)
        y(i)=y(i)+dcoo%a(inz)*x(dcoo%j(inz))*zs
    enddo
    return

    end subroutine ds_coomux


    subroutine zs_dcoohmux(zs,dcoo,x,y)
    type(dcoo_matrix),intent(in)::dcoo
    complex(q),intent(in)::zs,x(*)
    complex(q),intent(inout)::y(*)
    
    integer inz,j

    do inz=1,dcoo%nnz
        j=dcoo%j(inz)
        y(j)=y(j)+dcoo%a(inz)*x(dcoo%i(inz))*zs
    enddo
    return

    end subroutine zs_dcoohmux


    subroutine ds_dcoohmux(zs,dcoo,x,y)
    type(dcoo_matrix),intent(in)::dcoo
    real(q),intent(in)::zs,x(*)
    real(q),intent(inout)::y(*)
    
    integer inz,j

    do inz=1,dcoo%nnz
        j=dcoo%j(inz)
        y(j)=y(j)+dcoo%a(inz)*x(dcoo%i(inz))*zs
    enddo
    return

    end subroutine ds_dcoohmux


    subroutine zvh_dcoo_zv(dcoo,x,y)
    type(dcoo_matrix),intent(in)::dcoo
    complex(q),intent(in)::x(*)
    complex(q),intent(inout)::y

    integer inz,i,j

    do inz=1,dcoo%nnz
        i=dcoo%i(inz)
        j=dcoo%j(inz)
        y=y+conjg(x(i))*dcoo%a(inz)*x(j)
    enddo
    return

    end subroutine zvh_dcoo_zv


    subroutine dvh_coo_v(dcoo,x,y)
    type(dcoo_matrix),intent(in)::dcoo
    real(q),intent(in)::x(*)
    real(q),intent(inout)::y

    integer inz,i,j

    do inz=1,dcoo%nnz
        i=dcoo%i(inz)
        j=dcoo%j(inz)
        y=y+x(i)*dcoo%a(inz)*x(j)
    enddo
    return

    end subroutine dvh_coo_v


    subroutine zv2h_sdcoo_zv1(zs,dcoo,x1,x2,y)
    type(dcoo_matrix),intent(in)::dcoo
    complex(q),intent(in)::zs,x1(*),x2(*)
    complex(q),intent(inout)::y

    integer inz,i,j

    do inz=1,dcoo%nnz
        i=dcoo%i(inz)
        j=dcoo%j(inz)
        y=y+zs*conjg(x2(i))*dcoo%a(inz)*x1(j)
    enddo
    return

    end subroutine zv2h_sdcoo_zv1


    subroutine dv2h_scoo_v1(zs,dcoo,x1,x2,y)
    type(dcoo_matrix),intent(in)::dcoo
    real(q),intent(in)::zs,x1(*),x2(*)
    real(q),intent(inout)::y

    integer inz,i,j

    do inz=1,dcoo%nnz
        i=dcoo%i(inz)
        j=dcoo%j(inz)
        y=y+zs*x2(i)*dcoo%a(inz)*x1(j)
    enddo
    return

    end subroutine dv2h_scoo_v1


    subroutine zv2h_sdcooh_zv1(zs,dcoo,x1,x2,y)
    type(dcoo_matrix),intent(in)::dcoo
    complex(q),intent(in)::zs,x1(*),x2(*)
    complex(q),intent(inout)::y

    integer inz,i,j

    do inz=1,dcoo%nnz
        i=dcoo%i(inz)
        j=dcoo%j(inz)
        y=y+zs*conjg(x2(j))*dcoo%a(inz)*x1(i)
    enddo
    return

    end subroutine zv2h_sdcooh_zv1


    subroutine dv2h_scooh_v1(zs,dcoo,x1,x2,y)
    type(dcoo_matrix),intent(in)::dcoo
    real(q),intent(in)::zs,x1(*),x2(*)
    real(q),intent(inout)::y

    integer inz,i,j

    do inz=1,dcoo%nnz
        i=dcoo%i(inz)
        j=dcoo%j(inz)
        y=y+zs*x2(j)*dcoo%a(inz)*x1(i)
    enddo
    return

    end subroutine dv2h_scooh_v1


    subroutine zvh_dcsr_zv(dcsr,x,y)
    type(dcsr_matrix),intent(in)::dcsr
    complex(q),intent(in)::x(*)
    complex(q),intent(inout)::y

    integer inz,i,j,k

    do i=1,dcsr%nrow; do k=dcsr%i(i),dcsr%i(i+1)-1
        j=dcsr%j(k)
        y=y+conjg(x(i))*dcsr%a(k)*x(j)
    enddo; enddo
    return

    end subroutine zvh_dcsr_zv


    subroutine dvh_csr_v(dcsr,x,y)
    type(dcsr_matrix),intent(in)::dcsr
    real(q),intent(in)::x(*)
    real(q),intent(inout)::y

    integer inz,i,j,k

    do i=1,dcsr%nrow; do k=dcsr%i(i),dcsr%i(i+1)-1
        j=dcsr%j(k)
        y=y+x(i)*dcsr%a(k)*x(j)
    enddo; enddo
    return

    end subroutine dvh_csr_v


    subroutine zvh_csr_v(zcsr,x,y)
    type(zcsr_matrix),intent(in)::zcsr
    complex(q),intent(in)::x(*)
    complex(q),intent(inout)::y

    integer inz,i,j,k

    do i=1,zcsr%nrow; do k=zcsr%i(i),zcsr%i(i+1)-1
        j=zcsr%j(k)
        y=y+conjg(x(i))*zcsr%a(k)*x(j)
    enddo; enddo
    return

    end subroutine zvh_csr_v


    subroutine zvh_sydcsr_zv(dcsr,x,y)
    type(dcsr_matrix),intent(in)::dcsr
    complex(q),intent(in)::x(*)
    complex(q),intent(inout)::y

    integer inz,i,j,k

    do i=1,dcsr%nrow; do k=dcsr%i(i),dcsr%i(i+1)-1
        j=dcsr%j(k)
        y=y+conjg(x(i))*dcsr%a(k)*x(j)
        if(i==j)cycle
        y=y+conjg(x(j))*dcsr%a(k)*x(i)
    enddo; enddo
    return

    end subroutine zvh_sydcsr_zv


    subroutine dvh_sycsr_v(dcsr,x,y)
    type(dcsr_matrix),intent(in)::dcsr
    real(q),intent(in)::x(*)
    real(q),intent(inout)::y

    integer inz,i,j,k

    do i=1,dcsr%nrow; do k=dcsr%i(i),dcsr%i(i+1)-1
        j=dcsr%j(k)
        y=y+x(i)*dcsr%a(k)*x(j)
        if(i==j)cycle
        y=y+x(j)*dcsr%a(k)*x(i)
    enddo; enddo
    return

    end subroutine dvh_sycsr_v


end module sparse
