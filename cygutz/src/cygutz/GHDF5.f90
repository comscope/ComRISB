!******************************************************************************
! Author: Yong-Xin Yao <cygutz@gmail.com>
! Function: hdf5 object handling parallel independent and collective io.
!******************************************************************************

module ghdf5
    use gprec
    use hdf5
    implicit none
    include 'mpif.h'
    private

    type,public::hdf5_ob
        integer :: nmax
        !< file identifier
        integer(hid_t),private,allocatable :: f_id(:)
        !< preset 7 ids to be used.
        integer(hid_t),private,allocatable :: dset_id(:)
        !< to save a group id.
        integer(hid_t),private,allocatable :: group_id(:)
        !< predefined types
        integer(hid_t),public :: dcomplex_id, string_id
        logical,private :: done_init=.false.

        contains
        procedure::init=>gh5_init
        procedure::end=>gh5_end
        procedure::gcreate=>gh5_gcreate
        procedure::dcreate=>gh5_dcreate
        procedure::dopen=>gh5_dopen
        procedure::fopen=>gh5_fopen
        procedure::fclose=>gh5_fclose
        procedure::gopen=>gh5_gopen
        procedure::gclose=>gh5_gclose
        procedure::dclose=>gh5_dclose
        procedure::exists=>gh5_exists
        procedure,private::awrite_iscalar  =>gh5_awrite_iscalar
        procedure,private::awrite_dscalar  =>gh5_awrite_dscalar
        procedure,private::awrite_zscalar  =>gh5_awrite_zscalar
        procedure,private::awrite_sscalar  =>gh5_awrite_sscalar
        procedure,private::awrite_1d_iarray=>gh5_awrite_1d_iarray
        procedure,private::awrite_2d_iarray=>gh5_awrite_2d_iarray
        procedure,private::awrite_3d_iarray=>gh5_awrite_3d_iarray
        procedure,private::awrite_1d_darray=>gh5_awrite_1d_darray
        procedure,private::awrite_2d_darray=>gh5_awrite_2d_darray
        procedure,private::awrite_3d_darray=>gh5_awrite_3d_darray
        procedure,private::awrite_1d_zarray=>gh5_awrite_1d_zarray
        procedure,private::awrite_2d_zarray=>gh5_awrite_2d_zarray
        procedure,private::awrite_3d_zarray=>gh5_awrite_3d_zarray
        procedure,private::awrite_1d_sarray=>gh5_awrite_1d_sarray
        procedure,private::awrite_2d_sarray=>gh5_awrite_2d_sarray
        procedure,private::awrite_3d_sarray=>gh5_awrite_3d_sarray
        procedure,private::dwrite_1d_iarray     =>gh5_dwrite_1d_iarray
        procedure,private::dwrite_2d_iarray     =>gh5_dwrite_2d_iarray
        procedure,private::dwrite_3d_iarray     =>gh5_dwrite_3d_iarray
        procedure,private::dwrite_4d_iarray     =>gh5_dwrite_4d_iarray
        procedure,private::dwrite_1d_darray     =>gh5_dwrite_1d_darray
        procedure,private::dwrite_2d_darray     =>gh5_dwrite_2d_darray
        procedure,private::dwrite_3d_darray     =>gh5_dwrite_3d_darray
        procedure,private::dwrite_4d_darray     =>gh5_dwrite_4d_darray
        procedure,private::dwrite_1d_zarray     =>gh5_dwrite_1d_zarray
        procedure,private::dwrite_2d_zarray     =>gh5_dwrite_2d_zarray
        procedure,private::dwrite_3d_zarray     =>gh5_dwrite_3d_zarray
        procedure,private::dwrite_4d_zarray     =>gh5_dwrite_4d_zarray
        procedure,private::dwrite_1d_sarray     =>gh5_dwrite_1d_sarray
        procedure,private::dwrite_2d_sarray     =>gh5_dwrite_2d_sarray
        procedure,private::dwrite_3d_sarray     =>gh5_dwrite_3d_sarray
        procedure,private::dwrite_4d_sarray     =>gh5_dwrite_4d_sarray
        procedure,private::dwrite_1d_iarray_slab=>gh5_dwrite_1d_iarray_slab
        procedure,private::dwrite_2d_iarray_slab=>gh5_dwrite_2d_iarray_slab
        procedure,private::dwrite_3d_iarray_slab=>gh5_dwrite_3d_iarray_slab
        procedure,private::dwrite_4d_iarray_slab=>gh5_dwrite_4d_iarray_slab
        procedure,private::dwrite_1d_darray_slab=>gh5_dwrite_1d_darray_slab
        procedure,private::dwrite_2d_darray_slab=>gh5_dwrite_2d_darray_slab
        procedure,private::dwrite_3d_darray_slab=>gh5_dwrite_3d_darray_slab
        procedure,private::dwrite_4d_darray_slab=>gh5_dwrite_4d_darray_slab
        procedure,private::dwrite_5d_darray_slab=>gh5_dwrite_5d_darray_slab
        procedure,private::dwrite_1d_zarray_slab=>gh5_dwrite_1d_zarray_slab
        procedure,private::dwrite_2d_zarray_slab=>gh5_dwrite_2d_zarray_slab
        procedure,private::dwrite_3d_zarray_slab=>gh5_dwrite_3d_zarray_slab
        procedure,private::dwrite_4d_zarray_slab=>gh5_dwrite_4d_zarray_slab
        procedure,private::dwrite_5d_zarray_slab=>gh5_dwrite_5d_zarray_slab
        procedure,private::dwrite_1d_sarray_slab=>gh5_dwrite_1d_sarray_slab
        procedure,private::dwrite_2d_sarray_slab=>gh5_dwrite_2d_sarray_slab
        procedure,private::dwrite_3d_sarray_slab=>gh5_dwrite_3d_sarray_slab
        procedure,private::dwrite_4d_sarray_slab=>gh5_dwrite_4d_sarray_slab
        procedure,private::aread_i  =>gh5_aread_i
        procedure,private::aread_d  =>gh5_aread_d
        procedure,private::aread_z  =>gh5_aread_z
        procedure,private::aread_s  =>gh5_aread_s
        procedure,private::dread_1d_iarray     =>gh5_dread_1d_iarray
        procedure,private::dread_2d_iarray     =>gh5_dread_2d_iarray
        procedure,private::dread_3d_iarray     =>gh5_dread_3d_iarray
        procedure,private::dread_4d_iarray     =>gh5_dread_4d_iarray
        procedure,private::dread_1d_darray     =>gh5_dread_1d_darray
        procedure,private::dread_2d_darray     =>gh5_dread_2d_darray
        procedure,private::dread_3d_darray     =>gh5_dread_3d_darray
        procedure,private::dread_4d_darray     =>gh5_dread_4d_darray
        procedure,private::dread_1d_zarray     =>gh5_dread_1d_zarray
        procedure,private::dread_2d_zarray     =>gh5_dread_2d_zarray
        procedure,private::dread_3d_zarray     =>gh5_dread_3d_zarray
        procedure,private::dread_4d_zarray     =>gh5_dread_4d_zarray
        procedure,private::dread_1d_sarray     =>gh5_dread_1d_sarray
        procedure,private::dread_2d_sarray     =>gh5_dread_2d_sarray
        procedure,private::dread_3d_sarray     =>gh5_dread_3d_sarray
        procedure,private::dread_4d_sarray     =>gh5_dread_4d_sarray
        procedure,private::dread_1d_iarray_slab=>gh5_dread_1d_iarray_slab
        procedure,private::dread_2d_iarray_slab=>gh5_dread_2d_iarray_slab
        procedure,private::dread_3d_iarray_slab=>gh5_dread_3d_iarray_slab
        procedure,private::dread_4d_iarray_slab=>gh5_dread_4d_iarray_slab
        procedure,private::dread_1d_darray_slab=>gh5_dread_1d_darray_slab
        procedure,private::dread_2d_darray_slab=>gh5_dread_2d_darray_slab
        procedure,private::dread_3d_darray_slab=>gh5_dread_3d_darray_slab
        procedure,private::dread_4d_darray_slab=>gh5_dread_4d_darray_slab
        procedure,private::dread_5d_darray_slab=>gh5_dread_5d_darray_slab
        procedure,private::dread_1d_zarray_slab=>gh5_dread_1d_zarray_slab
        procedure,private::dread_2d_zarray_slab=>gh5_dread_2d_zarray_slab
        procedure,private::dread_3d_zarray_slab=>gh5_dread_3d_zarray_slab
        procedure,private::dread_4d_zarray_slab=>gh5_dread_4d_zarray_slab
        procedure,private::dread_5d_zarray_slab=>gh5_dread_5d_zarray_slab
        procedure,private::dread_1d_sarray_slab=>gh5_dread_1d_sarray_slab
        procedure,private::dread_2d_sarray_slab=>gh5_dread_2d_sarray_slab
        procedure,private::dread_3d_sarray_slab=>gh5_dread_3d_sarray_slab
        procedure,private::dread_4d_sarray_slab=>gh5_dread_4d_sarray_slab
        generic,public:: awrite=> awrite_iscalar, &
                & awrite_dscalar, &
                & awrite_zscalar, &
                & awrite_sscalar, &
                & awrite_1d_iarray, &
                & awrite_2d_iarray, &
                & awrite_3d_iarray, &
                & awrite_1d_darray, &
                & awrite_2d_darray, &
                & awrite_3d_darray, &
                & awrite_1d_zarray, &
                & awrite_2d_zarray, &
                & awrite_3d_zarray, &
                & awrite_1d_sarray, &
                & awrite_2d_sarray, &
                & awrite_3d_sarray
        generic,public:: dwrite=> dwrite_1d_iarray, &
                & dwrite_2d_iarray, &
                & dwrite_3d_iarray, &
                & dwrite_4d_iarray, &
                & dwrite_1d_darray, &
                & dwrite_2d_darray, &
                & dwrite_3d_darray, &
                & dwrite_4d_darray, &
                & dwrite_1d_zarray, &
                & dwrite_2d_zarray, &
                & dwrite_3d_zarray, &
                & dwrite_4d_zarray, &
                & dwrite_1d_sarray, &
                & dwrite_2d_sarray, &
                & dwrite_3d_sarray, &
                & dwrite_4d_sarray, &
                & dwrite_1d_iarray_slab, &
                & dwrite_2d_iarray_slab, &
                & dwrite_3d_iarray_slab, &
                & dwrite_4d_iarray_slab, &
                & dwrite_1d_darray_slab, &
                & dwrite_2d_darray_slab, &
                & dwrite_3d_darray_slab, &
                & dwrite_4d_darray_slab, &
                & dwrite_5d_darray_slab, &
                & dwrite_1d_zarray_slab, &
                & dwrite_2d_zarray_slab, &
                & dwrite_3d_zarray_slab, &
                & dwrite_4d_zarray_slab, &
                & dwrite_5d_zarray_slab, &
                & dwrite_1d_sarray_slab, &
                & dwrite_2d_sarray_slab, &
                & dwrite_3d_sarray_slab, &
                & dwrite_4d_sarray_slab
        generic,public::read=>aread_i, &
                & aread_d, &
                & aread_z, &
                & aread_s, &
                & dread_1d_iarray, &
                & dread_2d_iarray, &
                & dread_3d_iarray, &
                & dread_4d_iarray, &
                & dread_1d_darray, &
                & dread_2d_darray, &
                & dread_3d_darray, &
                & dread_4d_darray, &
                & dread_1d_zarray, &
                & dread_2d_zarray, &
                & dread_3d_zarray, &
                & dread_4d_zarray, &
                & dread_1d_sarray, &
                & dread_2d_sarray, &
                & dread_3d_sarray, &
                & dread_4d_sarray, &
                & dread_1d_iarray_slab, &
                & dread_2d_iarray_slab, &
                & dread_3d_iarray_slab, &
                & dread_4d_iarray_slab, &
                & dread_1d_darray_slab, &
                & dread_2d_darray_slab, &
                & dread_3d_darray_slab, &
                & dread_4d_darray_slab, &
                & dread_5d_darray_slab, &
                & dread_1d_zarray_slab, &
                & dread_2d_zarray_slab, &
                & dread_3d_zarray_slab, &
                & dread_4d_zarray_slab, &
                & dread_5d_zarray_slab, &
                & dread_1d_sarray_slab, &
                & dread_2d_sarray_slab, &
                & dread_3d_sarray_slab, &
                & dread_4d_sarray_slab
    end type hdf5_ob

    interface
        subroutine abort() bind(C, name="abort")
        end subroutine
    end interface

    contains
    
    
    subroutine gh5_init(this,nmax)
    class(hdf5_ob) :: this
    integer,optional::nmax

    integer :: ierr

    if(present(nmax))then
        this%nmax=nmax
    else
        this%nmax=7
    endif
      
    allocate(this%f_id(this%nmax), this%dset_id(this%nmax), &
            &this%group_id(this%nmax))
    this%f_id=-1; this%dset_id=-1; this%group_id=-1
    !< initialize fortran interface.
    call h5open_f(ierr)
    call stop_err(ierr, "h5open_f")

    !< create a new file using default properties.
    call set_dcomplex_id(this)
    call h5tcopy_f(h5t_fortran_s1, this%string_id, ierr)
    call stop_err(ierr, "h5tcopy_f")
    this%done_init=.true.
    return
      
    end subroutine gh5_init
      

    subroutine gh5_fopen(this, gh5_file, i, wmode, serialio)
    class(hdf5_ob) :: this
    character(*),intent(in) :: gh5_file,wmode
    integer,intent(in) :: i
    logical,optional,intent(in) :: serialio

    integer :: ierr
    !< property list identifier
    integer(hid_t) :: plist_id
    integer :: mode
    logical :: serial_io=.false.
      
    if(.not.this%done_init)then
        call gh5_init(this)
    endif

    if(this%f_id(i)/=-1)then
        write(0,'("error : hdf5 file handle ", i0, " is in use!")') i
        call abort()
    endif

    if(wmode=="rw".or.wmode=="wr")then
        mode=h5f_acc_rdwr_f
    elseif(wmode=="r")then
        mode=h5f_acc_rdonly_f
    elseif(wmode=="w")then
        mode=h5f_acc_trunc_f
    else
        write(0,'("error: illegal mode for gh5_open!")')
        call abort()
    endif

    if(present(serialio))then
        serial_io=serialio
    endif

    if(serial_io)then
        if(wmode=='w')then
            call h5fcreate_f(gh5_file, mode, this%f_id(i), ierr)
            call stop_err(ierr, "h5fcreate_f")
        else
            call h5fopen_f(gh5_file, mode, this%f_id(i), ierr)
            call stop_err(ierr, "h5fopen_f")
        endif
    else
        !< Setup file access property list with parallel I/O access.
        call h5pcreate_f(h5p_file_access_f, plist_id, ierr)
        call stop_err(ierr, "h5pcreate_f")
        call h5pset_fapl_mpio_f(plist_id, mpi_comm_world, mpi_info_null, &
                &ierr)
        call stop_err(ierr, "h5pset_fapl_mpio_f")
        if(wmode=='w')then
            call h5fcreate_f(gh5_file, mode, this%f_id(i), ierr, &
                    &access_prp = plist_id)
            call stop_err(ierr, "h5fcreate_f")
        else
            call h5fopen_f(gh5_file, mode, this%f_id(i), ierr, &
                    &access_prp = plist_id)
            call stop_err(ierr, "h5fopen_f")
        endif
        !< close the property list to frees resources.
        call h5pclose_f(plist_id, ierr)
        call stop_err(ierr, "h5pclose_f")
    endif
    return
      
    end subroutine gh5_fopen
    

    subroutine gh5_fclose(this, i)
    class(hdf5_ob) :: this
    integer,intent(in) :: i
   
    integer :: ierr

    !< Close the file.
    call h5fclose_f(this%f_id(i), ierr)
    call stop_err(ierr, "h5fclose_f")
    this%f_id(i)=-1
    return
      
    end subroutine gh5_fclose
    

    subroutine gh5_end(this)
    class(hdf5_ob) :: this

    integer :: ierr

    if(.not.this%done_init)return
    call h5tclose_f(this%dcomplex_id, ierr)
    call stop_err(ierr, "h5tclose_f")
    call h5tclose_f(this%string_id, ierr)
    call stop_err(ierr, "h5tclose_f")
    !< close fortran interface.
    call h5close_f(ierr)
    call stop_err(ierr, "h5close_f")
    deallocate(this%f_id, this%dset_id, this%group_id)
    return
      
    end subroutine gh5_end
    

    subroutine set_dcomplex_id(this)
    class(hdf5_ob) :: this

    integer(hsize_t) h_size, h_offset
    integer :: ierr

    h_size = q*2; h_offset = 0
    call h5tcreate_f(h5t_compound_f, h_size, this%dcomplex_id, ierr)
    call stop_err(ierr, "h5tcreate_f")
    call h5tinsert_f(this%dcomplex_id, "r", h_offset, &
            &h5kind_to_type(q, h5_real_kind), ierr)
    call stop_err(ierr, "h5tinsert_f")
    h_offset = h_offset + q
    call h5tinsert_f(this%dcomplex_id, "i", h_offset, &
            &h5kind_to_type(q, h5_real_kind), ierr)
    call stop_err(ierr, "h5tinsert_f")
    return
      
    end subroutine set_dcomplex_id
    

    subroutine gh5set_string_id_of_len(this, ns)
    class(hdf5_ob) :: this
    integer,intent(in) :: ns

    integer(hsize_t) h_size
    integer :: ierr

    h_size = ns
    call h5tset_size_f(this%string_id, h_size, ierr)
    call stop_err(ierr, "h5tset_size_f")
    return

    end subroutine gh5set_string_id_of_len


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! write scalar attribute
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine gh5_awrite_iscalar(this, i, a_name, i_g)
    class(hdf5_ob) :: this
    integer, intent(in), target :: i
    !< attribute name
    character(len=*), intent(in) :: a_name
    integer, intent(in) :: i_g
   
    type(c_ptr) :: gh5_ptr
 
    gh5_ptr = c_loc(i)
    call gh5_awrite_scalar(this, gh5_ptr, h5t_native_integer, a_name, &
            &i_g)
    return
      
    end subroutine gh5_awrite_iscalar


    subroutine gh5_awrite_dscalar(this, d, a_name, i)
    class(hdf5_ob) :: this
    real(q), intent(in), target :: d
    character(len=*), intent(in) :: a_name
    integer, intent(in) :: i

    type(c_ptr) :: gh5_ptr
    
    gh5_ptr = c_loc(d)
    call gh5_awrite_scalar(this, gh5_ptr, h5t_native_double, a_name, &
            &i)
    return
      
    end subroutine gh5_awrite_dscalar


    subroutine gh5_awrite_zscalar(this, z, a_name, i)
    class(hdf5_ob) :: this
    complex(q), intent(in), target :: z
    character(len=*), intent(in) :: a_name
    integer, intent(in) :: i

    type(c_ptr) :: gh5_ptr

    gh5_ptr = c_loc(z)
    call gh5_awrite_scalar(this, gh5_ptr, this%dcomplex_id, a_name, &
            &i)
    return

    end subroutine gh5_awrite_zscalar


    subroutine gh5_awrite_sscalar(this, s, ns, a_name, i)
    class(hdf5_ob) :: this
    integer, intent(in) :: ns, i
    character(ns), intent(in), target :: s
    character(len=*), intent(in) :: a_name

    type(c_ptr) :: gh5_ptr

    call gh5set_string_id_of_len(this, ns)
    gh5_ptr = c_loc(s(1:1))
    call gh5_awrite_scalar(this, gh5_ptr, this%string_id, a_name, &
            &i)
    return

    end subroutine gh5_awrite_sscalar


    subroutine gh5_awrite_scalar(this, s_ptr, h5t_id, a_name, i)
    class(hdf5_ob) :: this
    type(c_ptr), intent(in) :: s_ptr
    !< attribute name
    character(len=*), intent(in) :: a_name
    !< h5t_id defines the data type
    integer(hid_t), intent(in) :: h5t_id
    integer, intent(in) :: i
   
    integer(hid_t) :: attr_id, dspace_id
    integer :: ierr
    logical :: lexist

    call h5aexists_f(this%group_id(i), a_name, lexist, ierr)
    call stop_err(ierr, "h5aexists_f")
    if(lexist)then
        call h5aopen_name_f(this%group_id(i), a_name, attr_id, ierr) 
        call stop_err(ierr, "h5aopen_name_f")
    else
        !< Create scalar attribute
        call h5screate_f(h5s_scalar_f, dspace_id, ierr)
        call stop_err(ierr, "h5screate_f")
        call h5acreate_f(this%group_id(i), a_name, h5t_id, dspace_id, &
                &attr_id, ierr)
        call stop_err(ierr, "h5acreate_f")
    endif
    call h5awrite_f(attr_id, h5t_id, s_ptr, ierr)
    call stop_err(ierr, "h5awrite_f")
    ! Close the attribute.
    call h5aclose_f(attr_id, ierr)
    call stop_err(ierr, "h5aclose_f")
    if(.not.lexist)then
        ! Terminate access to the data space.
        call h5sclose_f(dspace_id, ierr)
        call stop_err(ierr, "h5sclose_f")
    endif
    return
      
    end subroutine gh5_awrite_scalar
    

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! write array attribute
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
    subroutine gh5_awrite_1d_iarray(this, iarray, n1, a_name, i)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, i
    integer, target, intent(in) :: iarray(n1)
    character(len=*), intent(in) :: a_name
    
    integer(hsize_t) :: h_dims(1)
    type(c_ptr) :: gh5_ptr
 
    h_dims = n1
    gh5_ptr = c_loc(iarray)
    call gh5_awrite_array(this, gh5_ptr, 1, h_dims, h5t_native_integer, &
            &a_name, i)
    return
      
    end subroutine gh5_awrite_1d_iarray


    subroutine gh5_awrite_2d_iarray(this, iarray, n1, n2, a_name, &
            &i)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, n2, i
    integer, target, intent(in) :: iarray(n1, n2)
    character(len=*), intent(in) :: a_name

    integer(hsize_t) :: h_dims(2)
    type(c_ptr) :: gh5_ptr
    
    h_dims = (/n1, n2/)
    gh5_ptr = c_loc(iarray)
    call gh5_awrite_array(this, gh5_ptr, 2, h_dims, h5t_native_integer, &
            &a_name, i)
    return
      
    end subroutine gh5_awrite_2d_iarray


    subroutine gh5_awrite_3d_iarray(this, iarray, n1, n2, n3, a_name, &
            &i)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, n2, n3, i
    integer, target, intent(in) :: iarray(n1, n2, n3)
    character(len=*), intent(in) :: a_name
   
    integer(hsize_t) :: h_dims(3)
    type(c_ptr) :: gh5_ptr

    h_dims = (/n1, n2, n3/)
    gh5_ptr = c_loc(iarray)
    call gh5_awrite_array(this, gh5_ptr, 3, h_dims, h5t_native_integer, &
            &a_name, i)
    return
      
    end subroutine gh5_awrite_3d_iarray


    subroutine gh5_awrite_1d_darray(this, darray, n1, a_name, i)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, i
    real(q), target, intent(in) :: darray(n1)
    character(len=*), intent(in) :: a_name
   
    integer(hsize_t) :: h_dims(1)
    type(c_ptr) :: gh5_ptr

    h_dims = n1
    gh5_ptr = c_loc(darray)
    call gh5_awrite_array(this, gh5_ptr, 1, h_dims, h5t_native_double, &
            &a_name, i)
    return
      
    end subroutine gh5_awrite_1d_darray
    

    subroutine gh5_awrite_2d_darray(this, darray, n1, n2, a_name, &
            &i)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, n2, i
    real(q), target, intent(in) :: darray(n1, n2)
    character(len=*), intent(in) :: a_name
    
    integer(hsize_t) :: h_dims(2)
    type(c_ptr) :: gh5_ptr

    h_dims = (/n1, n2/)
    gh5_ptr = c_loc(darray)
    call gh5_awrite_array(this, gh5_ptr, 2, h_dims, h5t_native_double, &
            &a_name, i)
    return
      
    end subroutine gh5_awrite_2d_darray


    subroutine gh5_awrite_3d_darray(this, darray, n1, n2, n3, a_name, &
            &i)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, n2, n3, i
    real(q), target, intent(in) :: darray(n1, n2, n3)
    character(len=*), intent(in) :: a_name
    
    integer(hsize_t) :: h_dims(3)
    type(c_ptr) :: gh5_ptr

    h_dims = (/n1, n2, n3/)
    gh5_ptr = c_loc(darray)
    call gh5_awrite_array(this, gh5_ptr, 3, h_dims, h5t_native_double, &
            &a_name, i)
    return
      
    end subroutine gh5_awrite_3d_darray


    subroutine gh5_awrite_1d_zarray(this, zarray, n1, a_name, i)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, i
    complex(q), target, intent(in) :: zarray(n1)
    character(len=*), intent(in) :: a_name
      
    integer(hsize_t) :: h_dims(1)
    type(c_ptr) :: gh5_ptr

    h_dims = n1
    gh5_ptr = c_loc(zarray)
    call gh5_awrite_array(this, gh5_ptr, 1, h_dims, this%dcomplex_id, &
            &a_name, i)
    return
      
    end subroutine gh5_awrite_1d_zarray


    subroutine gh5_awrite_2d_zarray(this, zarray, n1, n2, a_name, &
            &i)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, n2, i
    complex(q), target, intent(in) :: zarray(n1, n2)
    character(len=*), intent(in) :: a_name
      
    integer(hsize_t) :: h_dims(2)
    type(c_ptr) :: gh5_ptr

    h_dims = (/n1, n2/)
    gh5_ptr = c_loc(zarray)
    call gh5_awrite_array(this, gh5_ptr, 2, h_dims, this%dcomplex_id, &
            &a_name, i)
    return
      
    end subroutine gh5_awrite_2d_zarray


    subroutine gh5_awrite_3d_zarray(this, zarray, n1, n2, n3, a_name, &
            &i)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, n2, n3, i
    complex(q), target, intent(in) :: zarray(n1, n2, n3)
    character(len=*), intent(in) :: a_name
    
    integer(hsize_t) :: h_dims(3)
    type(c_ptr) :: gh5_ptr

    h_dims = (/n1, n2, n3/)
    gh5_ptr = c_loc(zarray)
    call gh5_awrite_array(this, gh5_ptr, 3, h_dims, this%dcomplex_id, &
            &a_name, i)
    return
      
    end subroutine gh5_awrite_3d_zarray


    subroutine gh5_awrite_1d_sarray(this, sarray, ns, n1, a_name, &
            &i)
    class(hdf5_ob) :: this
    integer, intent(in) :: ns, n1, i
    character(ns), target, intent(in) :: sarray(n1)
    character(len=*), intent(in) :: a_name
    
    integer(hsize_t) :: h_dims(1)
    type(c_ptr) :: gh5_ptr

    h_dims = n1
    call gh5set_string_id_of_len(this, ns)
    gh5_ptr = c_loc(sarray(1)(1:1))
    call gh5_awrite_array(this, gh5_ptr, 1, h_dims, this%string_id, &
            &a_name, i)
    return
      
    end subroutine gh5_awrite_1d_sarray


    subroutine gh5_awrite_2d_sarray(this, sarray, ns, n1, n2, a_name, &
            &i)
    class(hdf5_ob) :: this
    integer, intent(in) :: ns, n1, n2, i
    character(ns), target, intent(in) :: sarray(n1, n2)
    character(len=*), intent(in) :: a_name
   
    integer(hsize_t) :: h_dims(2)
    type(c_ptr) :: gh5_ptr

    h_dims = (/n1, n2/)
    call gh5set_string_id_of_len(this, ns)
    gh5_ptr = c_loc(sarray(1,1)(1:1))
    call gh5_awrite_array(this, gh5_ptr, 2, h_dims, this%string_id, &
            &a_name, i)
    return
      
    end subroutine gh5_awrite_2d_sarray


    subroutine gh5_awrite_3d_sarray(this, sarray, ns, n1, n2, n3, a_name, &
            &i)
    class(hdf5_ob) :: this
    integer, intent(in) :: ns, n1, n2, n3, i
    character(ns), target, intent(in) :: sarray(n1, n2, n3)
    character(len=*), intent(in) :: a_name
   
    integer(hsize_t) :: h_dims(3)
    type(c_ptr) :: gh5_ptr

    h_dims = (/n1, n2, n3/)
    call gh5set_string_id_of_len(this, ns)
    gh5_ptr = c_loc(sarray(1,1,1)(1:1))
    call gh5_awrite_array(this, gh5_ptr, 3, h_dims, this%string_id, &
            &a_name, i)
    return
      
    end subroutine gh5_awrite_3d_sarray


    subroutine gh5_awrite_array(this, f_ptr, n, h_dims, h5t_id, &
            &a_name, i)
    class(hdf5_ob) :: this
    type(c_ptr), intent(in) :: f_ptr
    integer, intent(in) :: n, i
    integer(hsize_t), intent(in) :: h_dims(n)
    character(len=*), intent(in) :: a_name
    integer(hid_t), intent(in) :: h5t_id
     
    integer(hid_t) :: attr_id, dspace_id
    integer :: ierr
    logical :: lexist

    call h5aexists_f(this%group_id(i), a_name, lexist, ierr)
    call stop_err(ierr, "h5aexists_f")
    if(lexist)then
        call h5aopen_name_f(this%group_id(i), a_name, attr_id, ierr)
        call stop_err(ierr, "h5aopen_name_f")
    else
        ! create the dataspace.
        call h5screate_simple_f(n, h_dims, dspace_id, ierr)
        call stop_err(ierr, "h5screate_simple_f")
        ! create dataset attribute.
        call h5acreate_f(this%group_id(i), a_name, h5t_id, dspace_id, &
                &attr_id, ierr)
        call stop_err(ierr, "h5acreate_f")
    endif
    ! write the attribute data.
    call h5awrite_f(attr_id, h5t_id, f_ptr, ierr)
    call stop_err(ierr, "h5awrite_f")
    ! close the attribute.
    call h5aclose_f(attr_id, ierr)
    call stop_err(ierr, "h5aclose_f")
    if(.not.lexist)then
       ! terminate access to the data space.
        call h5sclose_f(dspace_id, ierr)
        call stop_err(ierr, "h5sclose_f")
    endif
    return
      
    end subroutine gh5_awrite_array


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! read scalar attribute
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine gh5_aread_i(this, ibuf, a_name, i)
    class(hdf5_ob) :: this
    integer, intent(out), target :: ibuf
    character(len=*), intent(in) :: a_name
    integer, intent(in) :: i

    type(c_ptr) :: gh5_ptr
    
    gh5_ptr = c_loc(ibuf)
    call gh5_aread(this, gh5_ptr, h5t_native_integer, a_name, i)
    return
      
    end subroutine gh5_aread_i


    subroutine gh5_aread_d(this, dbuf, a_name, i)
    class(hdf5_ob) :: this
    real(q), intent(out), target :: dbuf
    character(len=*), intent(in) :: a_name
    integer, intent(in) :: i

    type(c_ptr) :: gh5_ptr
    
    gh5_ptr = c_loc(dbuf)
    call gh5_aread(this, gh5_ptr, h5t_native_double, a_name, i)
    return
      
    end subroutine gh5_aread_d


    subroutine gh5_aread_z(this, zbuf, a_name, i)
    class(hdf5_ob) :: this
    complex(q), intent(out), target :: zbuf
    character(len=*), intent(in) :: a_name
    integer, intent(in) :: i

    type(c_ptr) :: gh5_ptr
    
    gh5_ptr = c_loc(zbuf)
    call gh5_aread(this, gh5_ptr, this%dcomplex_id, a_name, i)
    return
      
    end subroutine gh5_aread_z
 

    subroutine gh5_aread_s(this, sbuf, ns, a_name, i)
    class(hdf5_ob) :: this
    integer, intent(in) :: ns
    character(ns), intent(out), target :: sbuf
    character(len=*), intent(in) :: a_name
    integer, intent(in) :: i

    type(c_ptr) :: gh5_ptr

    call gh5set_string_id_of_len(this, ns)
    gh5_ptr = c_loc(sbuf(1:1))
    call gh5_aread(this, gh5_ptr, this%string_id, a_name, i)
    return

    end subroutine gh5_aread_s


    subroutine gh5_aread(this, f_ptr, h5t_id, a_name, i)
    class(hdf5_ob) :: this
    type(c_ptr), intent(inout) :: f_ptr
    character(len=*), intent(in) :: a_name
    integer(hid_t), intent(in) :: h5t_id
    integer, intent(in) :: i
    
    integer(hid_t) :: attr_id
    integer :: ierr

    call h5aopen_f(this%group_id(i), a_name, attr_id, ierr)
    call stop_err(ierr, "h5aopen_f")
    call h5aread_f(attr_id, h5t_id, f_ptr, ierr)
    call stop_err(ierr, "h5aread_f")
    ! Close the attribute.
    call h5aclose_f(attr_id, ierr)
    call stop_err(ierr, "h5aclose_f")
    return
      
    end subroutine gh5_aread


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! write array dataset
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine gh5_dwrite_1d_iarray(this, iarray, n1, path, i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, i
    integer, target, intent(in) :: iarray(n1)
    character(len=*), intent(in) :: path
    logical, intent(in), optional :: serialio

    integer(hsize_t) :: h_dims(1)
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.
    
    h_dims = n1
    gh5_ptr = c_loc(iarray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dwrite_array(this, gh5_ptr, 1, h_dims, h5t_native_integer, &
            &i, serial_io, path=path)
    return
      
    end subroutine gh5_dwrite_1d_iarray


    subroutine gh5_dwrite_2d_iarray(this, iarray, n1, n2, path, i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, n2, i
    integer, target, intent(in) :: iarray(n1, n2)
    character(len=*), intent(in) :: path
    logical, intent(in), optional :: serialio

    integer(hsize_t) :: h_dims(2)
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.
    
    h_dims = (/n1, n2/)
    gh5_ptr = c_loc(iarray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dwrite_array(this, gh5_ptr, 2, h_dims, h5t_native_integer, &
            &i, serial_io, path=path)
    return
      
    end subroutine gh5_dwrite_2d_iarray


    subroutine gh5_dwrite_3d_iarray(this, iarray, n1, n2, n3, path, &
            &i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, n2, n3, i
    integer, target, intent(in) :: iarray(n1, n2, n3)
    character(len=*), intent(in) :: path
    logical, intent(in), optional :: serialio

    integer(hsize_t) :: h_dims(3)
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.
    
    h_dims = (/n1, n2, n3/)
    gh5_ptr = c_loc(iarray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dwrite_array(this, gh5_ptr, 3, h_dims, h5t_native_integer, &
            &i, serial_io, path=path)
    return
      
    end subroutine gh5_dwrite_3d_iarray


    subroutine gh5_dwrite_4d_iarray(this, iarray, n1, n2, n3, n4, path, &
            &i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, n2, n3, n4, i
    integer, target, intent(in) :: iarray(n1, n2, n3, n4)
    character(len=*), intent(in) :: path
    logical, intent(in), optional :: serialio

    integer(hsize_t) :: h_dims(4)
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.
    
    h_dims = (/n1, n2, n3, n4/)
    gh5_ptr = c_loc(iarray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dwrite_array(this, gh5_ptr, 4, h_dims, h5t_native_integer, &
            &i, serial_io, path=path)
    return
      
    end subroutine gh5_dwrite_4d_iarray


    subroutine gh5_dwrite_1d_darray(this, darray, n1, path, i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, i
    real(q), target, intent(in) :: darray(n1)
    character(len=*), intent(in) :: path
    logical, intent(in), optional :: serialio

    integer(hsize_t) :: h_dims(1)
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_dims = n1
    gh5_ptr = c_loc(darray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dwrite_array(this, gh5_ptr, 1, h_dims, h5t_native_double, &
            &i, serial_io, path=path)
    return

    end subroutine gh5_dwrite_1d_darray


    subroutine gh5_dwrite_2d_darray(this, darray, n1, n2, path, i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, n2, i
    real(q), target, intent(in) :: darray(n1, n2)
    character(len=*), intent(in) :: path
    logical, intent(in), optional :: serialio

    integer(hsize_t) :: h_dims(2)
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_dims = (/n1, n2/)
    gh5_ptr = c_loc(darray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dwrite_array(this, gh5_ptr, 2, h_dims, h5t_native_double, &
            &i, serial_io, path=path)
    return

    end subroutine gh5_dwrite_2d_darray


    subroutine gh5_dwrite_3d_darray(this, darray, n1, n2, n3, path, &
            &i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, n2, n3, i
    real(q), target, intent(in) :: darray(n1, n2, n3)
    character(len=*), intent(in) :: path
    logical, intent(in), optional :: serialio

    integer(hsize_t) :: h_dims(3)
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_dims = (/n1, n2, n3/)
    gh5_ptr = c_loc(darray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dwrite_array(this, gh5_ptr, 3, h_dims, h5t_native_double, &
            &i, serial_io, path=path)
    return

    end subroutine gh5_dwrite_3d_darray


    subroutine gh5_dwrite_4d_darray(this, darray, n1, n2, n3, n4, path, &
            &i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, n2, n3, n4, i
    real(q), target, intent(in) :: darray(n1, n2, n3, n4)
    character(len=*), intent(in) :: path
    logical, intent(in), optional :: serialio

    integer(hsize_t) :: h_dims(4)
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_dims = (/n1, n2, n3, n4/)
    gh5_ptr = c_loc(darray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dwrite_array(this, gh5_ptr, 4, h_dims, h5t_native_double, &
            &i, serial_io, path=path)
    return

    end subroutine gh5_dwrite_4d_darray


    subroutine gh5_dwrite_1d_zarray(this, zarray, n1, path, i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, i
    complex(q), target, intent(in) :: zarray(n1)
    character(len=*), intent(in) :: path
    logical, intent(in), optional :: serialio

    integer(hsize_t) :: h_dims(1)
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_dims = n1
    gh5_ptr = c_loc(zarray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dwrite_array(this, gh5_ptr, 1, h_dims, this%dcomplex_id, &
            &i, serial_io, path=path)
    return

    end subroutine gh5_dwrite_1d_zarray


    subroutine gh5_dwrite_2d_zarray(this, zarray, n1, n2, path, i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, n2, i
    complex(q), target, intent(in) :: zarray(n1, n2)
    character(len=*), intent(in) :: path
    logical, intent(in), optional :: serialio

    integer(hsize_t) :: h_dims(2)
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_dims = (/n1, n2/)
    gh5_ptr = c_loc(zarray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dwrite_array(this, gh5_ptr, 2, h_dims, this%dcomplex_id, &
            &i, serial_io, path=path)
    return

    end subroutine gh5_dwrite_2d_zarray


    subroutine gh5_dwrite_3d_zarray(this, zarray, n1, n2, n3, path, &
            &i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, n2, n3, i
    complex(q), target, intent(in) :: zarray(n1, n2, n3)
    character(len=*), intent(in) :: path
    logical, intent(in), optional :: serialio

    integer(hsize_t) :: h_dims(3)
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_dims = (/n1, n2, n3/)
    gh5_ptr = c_loc(zarray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dwrite_array(this, gh5_ptr, 3, h_dims, this%dcomplex_id, &
            &i, serial_io, path=path)
    return

    end subroutine gh5_dwrite_3d_zarray


    subroutine gh5_dwrite_4d_zarray(this, zarray, n1, n2, n3, n4, path, &
            &i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, n2, n3, n4, i
    complex(q), target, intent(in) :: zarray(n1, n2, n3, n4)
    character(len=*), intent(in) :: path
    logical, intent(in), optional :: serialio

    integer(hsize_t) :: h_dims(4)
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_dims = (/n1, n2, n3, n4/)
    gh5_ptr = c_loc(zarray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dwrite_array(this, gh5_ptr, 4, h_dims, this%dcomplex_id, &
            &i, serial_io, path=path)
    return

    end subroutine gh5_dwrite_4d_zarray


    subroutine gh5_dwrite_1d_sarray(this, sarray, ns, n1, path, i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: ns, n1, i
    character(ns), target, intent(in) :: sarray(n1)
    character(len=*), intent(in) :: path
    logical, intent(in), optional :: serialio

    integer(hsize_t) :: h_dims(1)
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_dims = n1
    call gh5set_string_id_of_len(this, ns)
    gh5_ptr = c_loc(sarray(1)(1:1))
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dwrite_array(this, gh5_ptr, 1, h_dims, this%string_id, &
            &i, serial_io, path=path)
    return

    end subroutine gh5_dwrite_1d_sarray


    subroutine gh5_dwrite_2d_sarray(this, sarray, ns, n1, n2, path, &
            &i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: ns, n1, n2, i
    character(ns), target, intent(in) :: sarray(n1, n2)
    character(len=*), intent(in) :: path
    logical, intent(in), optional :: serialio

    integer(hsize_t) :: h_dims(2)
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_dims = (/n1, n2/)
    call gh5set_string_id_of_len(this, ns)
    gh5_ptr = c_loc(sarray(1,1)(1:1))
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dwrite_array(this, gh5_ptr, 2, h_dims, this%string_id, &
            &i, serial_io, path=path)
    return

    end subroutine gh5_dwrite_2d_sarray


    subroutine gh5_dwrite_3d_sarray(this, sarray, ns, n1, n2, n3, path, &
            &i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: ns, n1, n2, n3, i
    character(ns), target, intent(in) :: sarray(n1, n2, n3)
    character(len=*), intent(in) :: path
    logical, intent(in), optional :: serialio

    integer(hsize_t) :: h_dims(3)
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_dims = (/n1, n2, n3/)
    call gh5set_string_id_of_len(this, ns)
    gh5_ptr = c_loc(sarray(1,1,1)(1:1))
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dwrite_array(this, gh5_ptr, 3, h_dims, this%string_id, &
            &i, serial_io, path=path)
    return

    end subroutine gh5_dwrite_3d_sarray


    subroutine gh5_dwrite_4d_sarray(this, sarray, ns, n1, n2, n3, n4, path, &
            &i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: ns, n1, n2, n3, n4, i
    character(ns), target, intent(in) :: sarray(n1, n2, n3, n4)
    character(len=*), intent(in) :: path
    logical, intent(in), optional :: serialio

    integer(hsize_t) :: h_dims(4)
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_dims = (/n1, n2, n3, n4/)
    call gh5set_string_id_of_len(this, ns)
    gh5_ptr = c_loc(sarray(1,1,1,1)(1:1))
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dwrite_array(this, gh5_ptr, 4, h_dims, this%string_id, &
            &i, serial_io, path=path)
    return

    end subroutine gh5_dwrite_4d_sarray


    subroutine gh5_dwrite_1d_iarray_slab(this, iarray, n1, &
            &o1, i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, o1, i
    integer, target, intent(in) :: iarray(n1)
    logical, intent(in), optional :: serialio

    integer(hsize_t), dimension(1) :: h_dims, h_offsets
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_offsets = (/o1/)
    h_dims = (/n1/)
    gh5_ptr = c_loc(iarray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dwrite_array(this, gh5_ptr, 1, h_dims, h5t_native_integer, &
            &i, serial_io, h_offsets=h_offsets)
    return

    end subroutine gh5_dwrite_1d_iarray_slab

    
    subroutine gh5_dwrite_2d_iarray_slab(this, iarray, n1, n2, &
            &o1, o2, i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, n2, o1, o2, i
    integer, target, intent(in) :: iarray(n1, n2)
    logical, intent(in), optional :: serialio

    integer(hsize_t), dimension(2) :: h_dims, h_offsets
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_offsets = (/o1, o2/)
    h_dims = (/n1, n2/)
    gh5_ptr = c_loc(iarray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dwrite_array(this, gh5_ptr, 2, h_dims, h5t_native_integer, &
            &i, serial_io, h_offsets=h_offsets)
    return

    end subroutine gh5_dwrite_2d_iarray_slab

   
    subroutine gh5_dwrite_3d_iarray_slab(this, iarray, n1, n2, n3, &
            &o1, o2, o3, i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, n2, n3, o1, o2, o3, i
    integer, target, intent(in) :: iarray(n1, n2, n3)
    logical, intent(in), optional :: serialio

    integer(hsize_t), dimension(3) :: h_dims, h_offsets
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_offsets = (/o1, o2, o3/)
    h_dims = (/n1, n2, n3/)
    gh5_ptr = c_loc(iarray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dwrite_array(this, gh5_ptr, 3, h_dims, h5t_native_integer, &
            &i, serial_io, h_offsets=h_offsets)
    return

    end subroutine gh5_dwrite_3d_iarray_slab
 

    subroutine gh5_dwrite_4d_iarray_slab(this, iarray, n1, n2, n3, n4, &
            &o1, o2, o3, o4, i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, n2, n3, n4, o1, o2, o3, o4, i
    integer, target, intent(in) :: iarray(n1, n2, n3, n4)
    logical, intent(in), optional :: serialio

    integer(hsize_t), dimension(4) :: h_dims, h_offsets
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_offsets = (/o1, o2, o3, o4/)
    h_dims = (/n1, n2, n3, n4/)
    gh5_ptr = c_loc(iarray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dwrite_array(this, gh5_ptr, 4, h_dims, h5t_native_integer, &
            &i, serial_io, h_offsets=h_offsets)
    return

    end subroutine gh5_dwrite_4d_iarray_slab


    subroutine gh5_dwrite_1d_darray_slab(this, darray, n1, &
            &o1, i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, o1, i
    real(q), target, intent(in) :: darray(n1)
    logical, intent(in), optional :: serialio

    integer(hsize_t), dimension(1) :: h_dims, h_offsets
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_offsets = (/o1/)
    h_dims = (/n1/)
    gh5_ptr = c_loc(darray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dwrite_array(this, gh5_ptr, 1, h_dims, h5t_native_double, &
            &i, serial_io, h_offsets=h_offsets)
    return

    end subroutine gh5_dwrite_1d_darray_slab

    
    subroutine gh5_dwrite_2d_darray_slab(this, darray, n1, n2, &
            &o1, o2, i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, n2, o1, o2, i
    real(q), target, intent(in) :: darray(n1, n2)
    logical, intent(in), optional :: serialio

    integer(hsize_t), dimension(2) :: h_dims, h_offsets
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_offsets = (/o1, o2/)
    h_dims = (/n1, n2/)
    gh5_ptr = c_loc(darray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dwrite_array(this, gh5_ptr, 2, h_dims, h5t_native_double, &
            &i, serial_io, h_offsets=h_offsets)
    return

    end subroutine gh5_dwrite_2d_darray_slab

   
    subroutine gh5_dwrite_3d_darray_slab(this, darray, n1, n2, n3, &
            &o1, o2, o3, i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, n2, n3, o1, o2, o3, i
    real(q), target, intent(in) :: darray(n1, n2, n3)
    logical, intent(in), optional :: serialio

    integer(hsize_t), dimension(3) :: h_dims, h_offsets
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_offsets = (/o1, o2, o3/)
    h_dims = (/n1, n2, n3/)
    gh5_ptr = c_loc(darray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dwrite_array(this, gh5_ptr, 3, h_dims, h5t_native_double, &
            &i, serial_io, h_offsets=h_offsets)
    return

    end subroutine gh5_dwrite_3d_darray_slab
 

    subroutine gh5_dwrite_4d_darray_slab(this, darray, n1, n2, n3, n4, &
            &o1, o2, o3, o4, i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, n2, n3, n4, o1, o2, o3, o4, i
    real(q), target, intent(in) :: darray(n1, n2, n3, n4)
    logical, intent(in), optional :: serialio

    integer(hsize_t), dimension(4) :: h_dims, h_offsets
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_offsets = (/o1, o2, o3, o4/)
    h_dims = (/n1, n2, n3, n4/)
    gh5_ptr = c_loc(darray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dwrite_array(this, gh5_ptr, 4, h_dims, h5t_native_double, &
            &i, serial_io, h_offsets=h_offsets)
    return

    end subroutine gh5_dwrite_4d_darray_slab
 

    subroutine gh5_dwrite_5d_darray_slab(this, darray, n1, n2, n3, n4, n5, &
            &o1, o2, o3, o4, o5, i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, n2, n3, n4, n5, o1, o2, o3, o4, o5, i
    real(q), target, intent(in) :: darray(n1, n2, n3, n4, n5)
    logical, intent(in), optional :: serialio

    integer(hsize_t), dimension(5) :: h_dims, h_offsets
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_offsets = (/o1, o2, o3, o4, o5/)
    h_dims = (/n1, n2, n3, n4, n5/)
    gh5_ptr = c_loc(darray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dwrite_array(this, gh5_ptr, 5, h_dims, h5t_native_double, &
            &i, serial_io, h_offsets=h_offsets)
    return

    end subroutine gh5_dwrite_5d_darray_slab


    subroutine gh5_dwrite_1d_zarray_slab(this, zarray, n1, &
            &o1, i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, o1, i
    complex(q), target, intent(in) :: zarray(n1)
    logical, intent(in), optional :: serialio

    integer(hsize_t), dimension(1) :: h_dims, h_offsets
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_offsets = (/o1/)
    h_dims = (/n1/)
    gh5_ptr = c_loc(zarray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dwrite_array(this, gh5_ptr, 1, h_dims, this%dcomplex_id, &
            &i, serial_io, h_offsets=h_offsets)
    return

    end subroutine gh5_dwrite_1d_zarray_slab

    
    subroutine gh5_dwrite_2d_zarray_slab(this, zarray, n1, n2, &
            &o1, o2, i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, n2, o1, o2, i
    complex(q), target, intent(in) :: zarray(n1, n2)
    logical, intent(in), optional :: serialio

    integer(hsize_t), dimension(2) :: h_dims, h_offsets
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_offsets = (/o1, o2/)
    h_dims = (/n1, n2/)
    gh5_ptr = c_loc(zarray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dwrite_array(this, gh5_ptr, 2, h_dims, this%dcomplex_id, &
            &i, serial_io, h_offsets=h_offsets)
    return

    end subroutine gh5_dwrite_2d_zarray_slab

   
    subroutine gh5_dwrite_3d_zarray_slab(this, zarray, n1, n2, n3, &
            &o1, o2, o3, i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, n2, n3, o1, o2, o3, i
    complex(q), target, intent(in) :: zarray(n1, n2, n3)
    logical, intent(in), optional :: serialio

    integer(hsize_t), dimension(3) :: h_dims, h_offsets
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_offsets = (/o1, o2, o3/)
    h_dims = (/n1, n2, n3/)
    gh5_ptr = c_loc(zarray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dwrite_array(this, gh5_ptr, 3, h_dims, this%dcomplex_id, &
            &i, serial_io, h_offsets=h_offsets)
    return

    end subroutine gh5_dwrite_3d_zarray_slab
 

    subroutine gh5_dwrite_4d_zarray_slab(this, zarray, n1, n2, n3, n4, &
            &o1, o2, o3, o4, i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, n2, n3, n4, o1, o2, o3, o4, i
    complex(q), target, intent(in) :: zarray(n1, n2, n3, n4)
    logical, intent(in), optional :: serialio

    integer(hsize_t), dimension(4) :: h_dims, h_offsets
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_offsets = (/o1, o2, o3, o4/)
    h_dims = (/n1, n2, n3, n4/)
    gh5_ptr = c_loc(zarray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dwrite_array(this, gh5_ptr, 4, h_dims, this%dcomplex_id, &
            &i, serial_io, h_offsets=h_offsets)
    return

    end subroutine gh5_dwrite_4d_zarray_slab
 

    subroutine gh5_dwrite_5d_zarray_slab(this, zarray, n1, n2, n3, n4, n5, &
            &o1, o2, o3, o4, o5, i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, n2, n3, n4, n5, o1, o2, o3, o4, o5, i
    complex(q), target, intent(in) :: zarray(n1, n2, n3, n4, n5)
    logical, intent(in), optional :: serialio

    integer(hsize_t), dimension(5) :: h_dims, h_offsets
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_offsets = (/o1, o2, o3, o4, o5/)
    h_dims = (/n1, n2, n3, n4, n5/)
    gh5_ptr = c_loc(zarray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dwrite_array(this, gh5_ptr, 5, h_dims, this%dcomplex_id, &
            &i, serial_io, h_offsets=h_offsets)
    return

    end subroutine gh5_dwrite_5d_zarray_slab


    subroutine gh5_dwrite_1d_sarray_slab(this, sarray, ns, n1, &
            &o1, i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: ns, n1, o1, i
    character(ns), target, intent(in) :: sarray(n1)
    logical, intent(in), optional :: serialio

    integer(hsize_t), dimension(1) :: h_dims, h_offsets
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_offsets = (/o1/)
    h_dims = (/n1/)
    call gh5set_string_id_of_len(this, ns)
    gh5_ptr = c_loc(sarray(1)(1:1))
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dwrite_array(this, gh5_ptr, 1, h_dims, this%string_id, &
            &i, serial_io, h_offsets=h_offsets)
    return

    end subroutine gh5_dwrite_1d_sarray_slab


    subroutine gh5_dwrite_2d_sarray_slab(this, sarray, ns, n1, n2, &
            &o1, o2, i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: ns, n1, n2, o1, o2, i
    character(ns), target, intent(in) :: sarray(n1, n2)
    logical, intent(in), optional :: serialio

    integer(hsize_t), dimension(2) :: h_dims, h_offsets
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_offsets = (/o1, o2/)
    h_dims = (/n1, n2/)
    call gh5set_string_id_of_len(this, ns)
    gh5_ptr = c_loc(sarray(1,1)(1:1))
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dwrite_array(this, gh5_ptr, 2, h_dims, this%string_id, &
            &i, serial_io, h_offsets=h_offsets)
    return

    end subroutine gh5_dwrite_2d_sarray_slab


    subroutine gh5_dwrite_3d_sarray_slab(this, sarray, ns, n1, n2, n3, &
            &o1, o2, o3, i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: ns, n1, n2, n3, o1, o2, o3, i
    character(ns), target, intent(in) :: sarray(n1, n2, n3)
    logical, intent(in), optional :: serialio

    integer(hsize_t), dimension(3) :: h_dims, h_offsets
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_offsets = (/o1, o2, o3/)
    h_dims = (/n1, n2, n3/)
    call gh5set_string_id_of_len(this, ns)
    gh5_ptr = c_loc(sarray(1,1,1)(1:1))
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dwrite_array(this, gh5_ptr, 3, h_dims, this%string_id, &
            &i, serial_io, h_offsets=h_offsets)
    return

    end subroutine gh5_dwrite_3d_sarray_slab
   

    subroutine gh5_dwrite_4d_sarray_slab(this, sarray, ns, n1, n2, n3, n4, &
            &o1, o2, o3, o4, i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: ns, n1, n2, n3, n4, o1, o2, o3, o4, i
    character(ns), target, intent(in) :: sarray(n1, n2, n3, n4)
    logical, intent(in), optional :: serialio

    integer(hsize_t), dimension(4) :: h_dims, h_offsets
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_offsets = (/o1, o2, o3, o4/)
    h_dims = (/n1, n2, n3, n4/)
    call gh5set_string_id_of_len(this, ns)
    gh5_ptr = c_loc(sarray(1,1,1,1)(1:1))
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dwrite_array(this, gh5_ptr, 4, h_dims, this%string_id, &
            &i, serial_io, h_offsets=h_offsets)
    return

    end subroutine gh5_dwrite_4d_sarray_slab

   
    subroutine gh5_dwrite_array(this, f_ptr, n, h_dims, h5t_id, i, &
            &serialio, path, h_offsets)
    class(hdf5_ob) :: this
    type(c_ptr), intent(in) :: f_ptr
    integer, intent(in) :: n, i
    integer(hsize_t), intent(in) :: h_dims(n)
    character(len=*), intent(in),optional :: path
    integer(hid_t), intent(in) :: h5t_id
    integer(hsize_t), intent(in), optional :: h_offsets(n)
    logical, intent(in) :: serialio
     
    integer(hid_t) :: fspace_id, dspace_id, dsetid, plist_id
    logical :: lexist
    integer :: ierr

    ! create the dataspace.
    call h5screate_simple_f(n, h_dims, dspace_id, ierr)
    call stop_err(ierr, "h5screate_simple_f")
    if(serialio)then
        plist_id=h5p_default_f
    else
        ! Create property list for collective dataset write
        call h5pcreate_f(h5p_dataset_xfer_f, plist_id, ierr)
        call stop_err(ierr, "h5pcreate_f")
        call h5pset_dxpl_mpio_f(plist_id, h5fd_mpio_collective_f, ierr)
        call stop_err(ierr, "h5pset_dxpl_mpio_f")
    endif
    if(.not.present(path))then
        call h5dget_space_f(this%dset_id(i), fspace_id, ierr)
        call stop_err(ierr, "h5dget_space_f")
        call h5sselect_hyperslab_f(fspace_id, h5s_select_set_f, &
                &h_offsets, h_dims, ierr)
        call stop_err(ierr, "h5sselect_hyperslab_f")
        ! Write the dataset collectively.
        call h5dwrite_f(this%dset_id(i), h5t_id, f_ptr, &
                &ierr, file_space_id=fspace_id, &
                &mem_space_id=dspace_id, xfer_prp=plist_id)
        call stop_err(ierr, "h5dwrite_f")
        ! terminate access to the data space.
        call h5sclose_f(fspace_id, ierr)
        call stop_err(ierr, "h5sclose_f")
    else
        call h5lexists_f(this%f_id(i), path, lexist, ierr)
        call stop_err(ierr, "h5lexists_f")
        if(lexist)then
            call h5dopen_f(this%f_id(i), path, dsetid, ierr)
            call stop_err(ierr, "h5dopen_f")
        else
            ! create the dataset with default properties.
            call h5dcreate_f(this%f_id(i), path, h5t_id, dspace_id, &
                    &dsetid, ierr)
            call stop_err(ierr, "h5dcreate_f")
        endif
        ! Write the dataset collectively.
        call h5dwrite_f(dsetid, h5t_id, f_ptr, ierr, xfer_prp=plist_id)
        call stop_err(ierr, "h5dwrite_f")
        ! end access to the dataset and release resources used by it.
        call h5dclose_f(dsetid, ierr)
        call stop_err(ierr, "h5dclose_f")
    endif
    ! terminate access to the data space.
    call h5sclose_f(dspace_id, ierr)
    call stop_err(ierr, "h5sclose_f")
    if(.not.serialio)then
        !< close the property list to frees resources.
        call h5pclose_f(plist_id, ierr)
        call stop_err(ierr, "h5pclose_f")
    endif
    return
      
    end subroutine gh5_dwrite_array


    subroutine gh5_dread_1d_iarray(this, iarray, n1, path, i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, i
    integer, target, intent(out) :: iarray(n1)
    character(len=*), intent(in) :: path
    logical, intent(in), optional :: serialio

    integer(hsize_t), dimension(1) :: h_dims
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_dims = n1
    gh5_ptr = c_loc(iarray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dread_array(this, gh5_ptr, 1, h_dims, h5t_native_integer, &
            &i, serial_io, path=path)
    return

    end subroutine gh5_dread_1d_iarray


    subroutine gh5_dread_2d_iarray(this, iarray, n1, n2, path, i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, n2, i
    integer, target, intent(out) :: iarray(n1, n2)
    character(len=*), intent(in) :: path
    logical, intent(in), optional :: serialio

    integer(hsize_t), dimension(2) :: h_dims
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_dims = (/n1, n2/)
    gh5_ptr = c_loc(iarray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dread_array(this, gh5_ptr, 2, h_dims, h5t_native_integer, &
            &i, serial_io, path=path)
    return

    end subroutine gh5_dread_2d_iarray


    subroutine gh5_dread_3d_iarray(this, iarray, n1, n2, n3, path, &
            &i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, n2, n3, i
    integer, target, intent(out) :: iarray(n1, n2, n3)
    character(len=*), intent(in) :: path
    logical, intent(in), optional :: serialio

    integer(hsize_t), dimension(3) :: h_dims
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_dims = (/n1, n2, n3/)
    gh5_ptr = c_loc(iarray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dread_array(this, gh5_ptr, 3, h_dims, h5t_native_integer, &
            &i, serial_io, path=path)
    return

    end subroutine gh5_dread_3d_iarray


    subroutine gh5_dread_4d_iarray(this, iarray, n1, n2, n3, n4, path, &
            &i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, n2, n3, n4, i
    integer, target, intent(out) :: iarray(n1, n2, n3, n4)
    character(len=*), intent(in) :: path
    logical, intent(in), optional :: serialio

    integer(hsize_t), dimension(4) :: h_dims
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_dims = (/n1, n2, n3, n4/)
    gh5_ptr = c_loc(iarray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dread_array(this, gh5_ptr, 4, h_dims, h5t_native_integer, &
            &i, serial_io, path=path)
    return

    end subroutine gh5_dread_4d_iarray


    subroutine gh5_dread_1d_darray(this, darray, n1, path, i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, i
    real(q), target, intent(out) :: darray(n1)
    character(len=*), intent(in) :: path
    logical, intent(in), optional :: serialio

    integer(hsize_t), dimension(1) :: h_dims
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_dims = n1
    gh5_ptr = c_loc(darray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dread_array(this, gh5_ptr, 1, h_dims, h5t_native_double, &
            &i, serial_io, path=path)
    return

    end subroutine gh5_dread_1d_darray


    subroutine gh5_dread_2d_darray(this, darray, n1, n2, path, &
            &i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, n2, i
    real(q), target, intent(out) :: darray(n1, n2)
    character(len=*), intent(in) :: path
    logical, intent(in), optional :: serialio

    integer(hsize_t), dimension(2) :: h_dims
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_dims = (/n1, n2/)
    gh5_ptr = c_loc(darray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dread_array(this, gh5_ptr, 2, h_dims, h5t_native_double, &
            &i, serial_io, path=path)
    return

    end subroutine gh5_dread_2d_darray


    subroutine gh5_dread_3d_darray(this, darray, n1, n2, n3, path, &
            &i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, n2, n3, i
    real(q), target, intent(out) :: darray(n1, n2, n3)
    character(len=*), intent(in) :: path
    logical, intent(in), optional :: serialio

    integer(hsize_t), dimension(3) :: h_dims
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_dims = (/n1, n2, n3/)
    gh5_ptr = c_loc(darray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dread_array(this, gh5_ptr, 3, h_dims, h5t_native_double, &
            &i, serial_io, path=path)
    return

    end subroutine gh5_dread_3d_darray


    subroutine gh5_dread_4d_darray(this, darray, n1, n2, n3, n4, path, &
            &i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, n2, n3, n4, i
    real(q), target, intent(out) :: darray(n1, n2, n3, n4)
    character(len=*), intent(in) :: path
    logical, intent(in), optional :: serialio

    integer(hsize_t), dimension(4) :: h_dims
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_dims = (/n1, n2, n3, n4/)
    gh5_ptr = c_loc(darray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dread_array(this, gh5_ptr, 4, h_dims, h5t_native_double, &
            &i, serial_io, path=path)
    return

    end subroutine gh5_dread_4d_darray


    subroutine gh5_dread_1d_zarray(this, zarray, n1, path, i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, i
    complex(q), target, intent(out) :: zarray(n1)
    character(len=*), intent(in) :: path
    logical, intent(in), optional :: serialio

    integer(hsize_t), dimension(1) :: h_dims
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_dims = n1
    gh5_ptr = c_loc(zarray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dread_array(this, gh5_ptr, 1, h_dims, this%dcomplex_id, &
            &i, serial_io, path=path)
    return

    end subroutine gh5_dread_1d_zarray


    subroutine gh5_dread_2d_zarray(this, zarray, n1, n2, path, i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, n2, i
    complex(q), target, intent(out) :: zarray(n1, n2)
    character(len=*), intent(in) :: path
    logical, intent(in), optional :: serialio

    integer(hsize_t), dimension(2) :: h_dims
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_dims = (/n1, n2/)
    gh5_ptr = c_loc(zarray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dread_array(this, gh5_ptr, 2, h_dims, this%dcomplex_id, &
            &i, serial_io, path=path)
    return

    end subroutine gh5_dread_2d_zarray


    subroutine gh5_dread_3d_zarray(this, zarray, n1, n2, n3, path, i, &
            &serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, n2, n3, i
    complex(q), target, intent(out) :: zarray(n1, n2, n3)
    character(len=*), intent(in) :: path
    logical, intent(in), optional :: serialio

    integer(hsize_t), dimension(3) :: h_dims
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_dims = (/n1, n2, n3/)
    gh5_ptr = c_loc(zarray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dread_array(this, gh5_ptr, 3, h_dims, this%dcomplex_id, &
            &i, serial_io, path=path)
    return

    end subroutine gh5_dread_3d_zarray


    subroutine gh5_dread_4d_zarray(this, zarray, n1, n2, n3, n4, path, &
            &i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, n2, n3, n4, i
    complex(q), target, intent(out) :: zarray(n1, n2, n3, n4)
    character(len=*), intent(in) :: path
    logical, intent(in), optional :: serialio

    integer(hsize_t), dimension(4) :: h_dims
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_dims = (/n1, n2, n3, n4/)
    gh5_ptr = c_loc(zarray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dread_array(this, gh5_ptr, 4, h_dims, this%dcomplex_id, &
            &i, serial_io, path=path)
    return

    end subroutine gh5_dread_4d_zarray


    subroutine gh5_dread_1d_sarray(this, sarray, ns, n1, path, i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: ns, n1, i
    character(ns), target, intent(out) :: sarray(n1)
    character(len=*), intent(in) :: path
    logical, intent(in), optional :: serialio

    integer(hsize_t) :: h_dims(1)
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_dims = n1
    call gh5set_string_id_of_len(this, ns)
    gh5_ptr = c_loc(sarray(1)(1:1))
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dread_array(this, gh5_ptr, 1, h_dims, this%string_id, &
            &i, serial_io, path=path)
    return

    end subroutine gh5_dread_1d_sarray


    subroutine gh5_dread_2d_sarray(this, sarray, ns, n1, n2, path, i, &
            &serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: ns, n1, n2, i
    character(ns), target, intent(out) :: sarray(n1, n2)
    character(len=*), intent(in) :: path
    logical, intent(in), optional :: serialio

    integer(hsize_t) :: h_dims(2)
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_dims = (/n1, n2/)
    call gh5set_string_id_of_len(this, ns)
    gh5_ptr = c_loc(sarray(1,1)(1:1))
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dread_array(this, gh5_ptr, 2, h_dims, this%string_id, &
            &i, serial_io, path=path)
    return

    end subroutine gh5_dread_2d_sarray


    subroutine gh5_dread_3d_sarray(this, sarray, ns, n1, n2, n3, path, i, &
            &serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: ns, n1, n2, n3, i
    character(ns), target, intent(out) :: sarray(n1, n2, n3)
    character(len=*), intent(in) :: path
    logical, intent(in), optional :: serialio

    integer(hsize_t) :: h_dims(3)
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_dims = (/n1, n2, n3/)
    call gh5set_string_id_of_len(this, ns)
    gh5_ptr = c_loc(sarray(1,1,1)(1:1))
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dread_array(this, gh5_ptr, 3, h_dims, this%string_id, &
            &i, serial_io, path=path)
    return

    end subroutine gh5_dread_3d_sarray


    subroutine gh5_dread_4d_sarray(this, sarray, ns, n1, n2, n3, n4, path, &
            &i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: ns, n1, n2, n3, n4, i
    character(ns), target, intent(out) :: sarray(n1, n2, n3, n4)
    character(len=*), intent(in) :: path
    logical, intent(in), optional :: serialio

    integer(hsize_t) :: h_dims(4)
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_dims = (/n1, n2, n3, n4/)
    call gh5set_string_id_of_len(this, ns)
    gh5_ptr = c_loc(sarray(1,1,1,1)(1:1))
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dread_array(this, gh5_ptr, 4, h_dims, this%string_id, &
            &i, serial_io, path=path)
    return

    end subroutine gh5_dread_4d_sarray


    subroutine gh5_dread_1d_iarray_slab(this, iarray, n1, &
            &o1, i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, o1, i
    integer, target, intent(out) :: iarray(n1)
    logical, intent(in), optional :: serialio

    integer(hsize_t), dimension(1) :: h_dims, h_offsets
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_offsets = (/o1/)
    h_dims = (/n1/)
    gh5_ptr = c_loc(iarray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dread_array(this, gh5_ptr, 1, h_dims, h5t_native_integer, &
            &i, serial_io, h_offsets=h_offsets)
    return

    end subroutine gh5_dread_1d_iarray_slab


    subroutine gh5_dread_2d_iarray_slab(this, iarray, n1, n2, &
            &o1, o2, i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, n2, o1, o2, i
    integer, target, intent(out) :: iarray(n1, n2)
    logical, intent(in), optional :: serialio

    integer(hsize_t), dimension(2) :: h_dims, h_offsets
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_offsets = (/o1, o2/)
    h_dims = (/n1, n2/)
    gh5_ptr = c_loc(iarray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dread_array(this, gh5_ptr, 2, h_dims, h5t_native_integer, &
            &i, serial_io, h_offsets=h_offsets)
    return

    end subroutine gh5_dread_2d_iarray_slab


    subroutine gh5_dread_3d_iarray_slab(this, iarray, n1, n2, n3, &
            &o1, o2, o3, i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, n2, n3, o1, o2, o3, i
    integer, target, intent(out) :: iarray(n1, n2, n3)
    logical, intent(in), optional :: serialio

    integer(hsize_t), dimension(3) :: h_dims, h_offsets
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_offsets = (/o1, o2, o3/)
    h_dims = (/n1, n2, n3/)
    gh5_ptr = c_loc(iarray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dread_array(this, gh5_ptr, 3, h_dims, h5t_native_integer, &
            &i, serial_io, h_offsets=h_offsets)
    return

    end subroutine gh5_dread_3d_iarray_slab

    
    subroutine gh5_dread_4d_iarray_slab(this, iarray, n1, n2, n3, n4, &
            &o1, o2, o3, o4, i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, n2, n3, n4, o1, o2, o3, o4, i
    integer, target, intent(out) :: iarray(n1, n2, n3, n4)
    logical, intent(in), optional :: serialio

    integer(hsize_t), dimension(4) :: h_dims, h_offsets
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_offsets = (/o1, o2, o3, o4/)
    h_dims = (/n1, n2, n3, n4/)
    gh5_ptr = c_loc(iarray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dread_array(this, gh5_ptr, 4, h_dims, h5t_native_integer, &
            &i, serial_io, h_offsets=h_offsets)
    return

    end subroutine gh5_dread_4d_iarray_slab


    subroutine gh5_dread_1d_darray_slab(this, darray, n1, &
            &o1, i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, o1, i
    real(q), target, intent(out) :: darray(n1)
    logical, intent(in), optional :: serialio

    integer(hsize_t), dimension(1) :: h_dims, h_offsets
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_offsets = (/o1/)
    h_dims = (/n1/)
    gh5_ptr = c_loc(darray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dread_array(this, gh5_ptr, 1, h_dims, h5t_native_double, &
            &i, serial_io, h_offsets=h_offsets)
    return

    end subroutine gh5_dread_1d_darray_slab


    subroutine gh5_dread_2d_darray_slab(this, darray, n1, n2, &
            &o1, o2, i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, n2, o1, o2, i
    real(q), target, intent(out) :: darray(n1, n2)
    logical, intent(in), optional :: serialio

    integer(hsize_t), dimension(2) :: h_dims, h_offsets
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_offsets = (/o1, o2/)
    h_dims = (/n1, n2/)
    gh5_ptr = c_loc(darray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dread_array(this, gh5_ptr, 2, h_dims, h5t_native_double, &
            &i, serial_io, h_offsets=h_offsets)
    return

    end subroutine gh5_dread_2d_darray_slab


    subroutine gh5_dread_3d_darray_slab(this, darray, n1, n2, n3, &
            &o1, o2, o3, i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, n2, n3, o1, o2, o3, i
    real(q), target, intent(out) :: darray(n1, n2, n3)
    logical, intent(in), optional :: serialio

    integer(hsize_t), dimension(3) :: h_dims, h_offsets
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_offsets = (/o1, o2, o3/)
    h_dims = (/n1, n2, n3/)
    gh5_ptr = c_loc(darray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dread_array(this, gh5_ptr, 3, h_dims, h5t_native_double, &
            &i, serial_io, h_offsets=h_offsets)
    return

    end subroutine gh5_dread_3d_darray_slab

    
    subroutine gh5_dread_4d_darray_slab(this, darray, n1, n2, n3, n4, &
            &o1, o2, o3, o4, i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, n2, n3, n4, o1, o2, o3, o4, i
    real(q), target, intent(out) :: darray(n1, n2, n3, n4)
    logical, intent(in), optional :: serialio

    integer(hsize_t), dimension(4) :: h_dims, h_offsets
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_offsets = (/o1, o2, o3, o4/)
    h_dims = (/n1, n2, n3, n4/)
    gh5_ptr = c_loc(darray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dread_array(this, gh5_ptr, 4, h_dims, h5t_native_double, &
            &i, serial_io, h_offsets=h_offsets)
    return

    end subroutine gh5_dread_4d_darray_slab

    
    subroutine gh5_dread_5d_darray_slab(this, darray, n1, n2, n3, n4, n5, &
            &o1, o2, o3, o4, o5, i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, n2, n3, n4, n5, o1, o2, o3, o4, o5, i
    real(q), target, intent(out) :: darray(n1, n2, n3, n4, n5)
    logical, intent(in), optional :: serialio

    integer(hsize_t), dimension(5) :: h_dims, h_offsets
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_offsets = (/o1, o2, o3, o4, o5/)
    h_dims = (/n1, n2, n3, n4, n5/)
    gh5_ptr = c_loc(darray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dread_array(this, gh5_ptr, 5, h_dims, h5t_native_double, &
            &i, serial_io, h_offsets=h_offsets)
    return

    end subroutine gh5_dread_5d_darray_slab


    subroutine gh5_dread_1d_zarray_slab(this, zarray, n1, &
            &o1, i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, o1, i
    complex(q), target, intent(out) :: zarray(n1)
    logical, intent(in), optional :: serialio

    integer(hsize_t), dimension(1) :: h_dims, h_offsets
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_offsets = (/o1/)
    h_dims = (/n1/)
    gh5_ptr = c_loc(zarray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dread_array(this, gh5_ptr, 1, h_dims, this%dcomplex_id, &
            &i, serial_io, h_offsets=h_offsets)
    return

    end subroutine gh5_dread_1d_zarray_slab


    subroutine gh5_dread_2d_zarray_slab(this, zarray, n1, n2, &
            &o1, o2, i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, n2, o1, o2, i
    complex(q), target, intent(out) :: zarray(n1, n2)
    logical, intent(in), optional :: serialio

    integer(hsize_t), dimension(2) :: h_dims, h_offsets
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_offsets = (/o1, o2/)
    h_dims = (/n1, n2/)
    gh5_ptr = c_loc(zarray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dread_array(this, gh5_ptr, 2, h_dims, this%dcomplex_id, &
            &i, serial_io, h_offsets=h_offsets)
    return

    end subroutine gh5_dread_2d_zarray_slab


    subroutine gh5_dread_3d_zarray_slab(this, zarray, n1, n2, n3, &
            &o1, o2, o3, i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, n2, n3, o1, o2, o3, i
    complex(q), target, intent(out) :: zarray(n1, n2, n3)
    logical, intent(in), optional :: serialio

    integer(hsize_t), dimension(3) :: h_dims, h_offsets
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_offsets = (/o1, o2, o3/)
    h_dims = (/n1, n2, n3/)
    gh5_ptr = c_loc(zarray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dread_array(this, gh5_ptr, 3, h_dims, this%dcomplex_id, &
            &i, serial_io, h_offsets=h_offsets)
    return

    end subroutine gh5_dread_3d_zarray_slab

    
    subroutine gh5_dread_4d_zarray_slab(this, zarray, n1, n2, n3, n4, &
            &o1, o2, o3, o4, i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, n2, n3, n4, o1, o2, o3, o4, i
    complex(q), target, intent(out) :: zarray(n1, n2, n3, n4)
    logical, intent(in), optional :: serialio

    integer(hsize_t), dimension(4) :: h_dims, h_offsets
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_offsets = (/o1, o2, o3, o4/)
    h_dims = (/n1, n2, n3, n4/)
    gh5_ptr = c_loc(zarray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dread_array(this, gh5_ptr, 4, h_dims, this%dcomplex_id, &
            &i, serial_io, h_offsets=h_offsets)
    return

    end subroutine gh5_dread_4d_zarray_slab


    subroutine gh5_dread_5d_zarray_slab(this, zarray, n1, n2, n3, n4, n5, &
            &o1, o2, o3, o4, o5, i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: n1, n2, n3, n4, n5, o1, o2, o3, o4, o5, i
    complex(q), target, intent(out) :: zarray(n1, n2, n3, n4, n5)
    logical, intent(in), optional :: serialio

    integer(hsize_t), dimension(5) :: h_dims, h_offsets
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_offsets = (/o1, o2, o3, o4, o5/)
    h_dims = (/n1, n2, n3, n4, n5/)
    gh5_ptr = c_loc(zarray)
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dread_array(this, gh5_ptr, 5, h_dims, this%dcomplex_id, &
            &i, serial_io, h_offsets=h_offsets)
    return

    end subroutine gh5_dread_5d_zarray_slab


    subroutine gh5_dread_1d_sarray_slab(this, sarray, ns, n1, &
            &o1, i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: ns, n1, o1, i
    character(ns), target, intent(out) :: sarray(n1)
    logical, intent(in), optional :: serialio

    integer(hsize_t), dimension(1) :: h_dims, h_offsets
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_offsets = (/o1/)
    h_dims = (/n1/)
    call gh5set_string_id_of_len(this, ns)
    gh5_ptr = c_loc(sarray(1)(1:1))
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dread_array(this, gh5_ptr, 1, h_dims, this%string_id, &
            &i, serial_io, h_offsets=h_offsets)
    return

    end subroutine gh5_dread_1d_sarray_slab


    subroutine gh5_dread_2d_sarray_slab(this, sarray, ns, n1, n2, &
            &o1, o2, i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: ns, n1, n2, o1, o2, i
    character(ns), target, intent(out) :: sarray(n1, n2)
    logical, intent(in), optional :: serialio

    integer(hsize_t), dimension(2) :: h_dims, h_offsets
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_offsets = (/o1, o2/)
    h_dims = (/n1, n2/)
    call gh5set_string_id_of_len(this, ns)
    gh5_ptr = c_loc(sarray(1,1)(1:1))
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dread_array(this, gh5_ptr, 2, h_dims, this%string_id, &
            &i, serial_io, h_offsets=h_offsets)
    return

    end subroutine gh5_dread_2d_sarray_slab


    subroutine gh5_dread_3d_sarray_slab(this, sarray, ns, n1, n2, n3, &
            &o1, o2, o3, i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: ns, n1, n2, n3, o1, o2, o3, i
    character(ns), target, intent(out) :: sarray(n1, n2, n3)
    logical, intent(in), optional :: serialio

    integer(hsize_t), dimension(3) :: h_dims, h_offsets
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_offsets = (/o1, o2, o3/)
    h_dims = (/n1, n2, n3/)
    call gh5set_string_id_of_len(this, ns)
    gh5_ptr = c_loc(sarray(1,1,1)(1:1))
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dread_array(this, gh5_ptr, 3, h_dims, this%string_id, &
            &i, serial_io, h_offsets=h_offsets)
    return

    end subroutine gh5_dread_3d_sarray_slab

    
    subroutine gh5_dread_4d_sarray_slab(this, sarray, ns, n1, n2, n3, n4, &
            &o1, o2, o3, o4, i, serialio)
    class(hdf5_ob) :: this
    integer, intent(in) :: ns, n1, n2, n3, n4, o1, o2, o3, o4, i
    character(ns), target, intent(out) :: sarray(n1, n2, n3, n4)
    logical, intent(in), optional :: serialio

    integer(hsize_t), dimension(4) :: h_dims, h_offsets
    type(c_ptr) :: gh5_ptr
    logical :: serial_io=.false.

    h_offsets = (/o1, o2, o3, o4/)
    h_dims = (/n1, n2, n3, n4/)
    call gh5set_string_id_of_len(this, ns)
    gh5_ptr = c_loc(sarray(1,1,1,1)(1:1))
    if(present(serialio))then
        serial_io=serialio
    endif
    call gh5_dread_array(this, gh5_ptr, 4, h_dims, this%string_id, &
            &i, serial_io, h_offsets=h_offsets)
    return

    end subroutine gh5_dread_4d_sarray_slab


    subroutine gh5_dread_array(this, f_ptr, n, h_dims, h5t_id, i, &
            &serialio, path, h_offsets)
    class(hdf5_ob) :: this
    type(c_ptr), intent(inout) :: f_ptr
    integer, intent(in) :: n, i ! index for either dset_id or f_id
    integer(hsize_t), intent(in) :: h_dims(n)
    character(len=*), intent(in), optional :: path
    integer(hid_t), intent(in) :: h5t_id
    integer(hsize_t), intent(in), optional :: h_offsets(n)
    logical, intent(in) :: serialio

    integer(hid_t) :: fspace_id, dspace_id, dsetid, plist_id
    integer :: ierr

    if(serialio)then
        plist_id=h5p_default_f
    else
        ! Create property list for collective dataset write
        call h5pcreate_f(h5p_dataset_xfer_f, plist_id, ierr)
        call stop_err(ierr, "h5pcreate_f")
        call h5pset_dxpl_mpio_f(plist_id, h5fd_mpio_collective_f, ierr)
        call stop_err(ierr, "h5pset_dxpl_mpio_f")
    endif
    if(.not.present(path))then
        ! create the dataspace.
        call h5screate_simple_f(n, h_dims, dspace_id, ierr)
        call stop_err(ierr, "h5screate_simple_f")
        ! Define and select the hyperslab to use for reading.
        call h5dget_space_f(this%dset_id(i), fspace_id, ierr)
        call stop_err(ierr, "h5dget_space_f")
        call h5sselect_hyperslab_f(fspace_id, h5s_select_set_f, &
                &h_offsets, h_dims, ierr)
        call stop_err(ierr, "h5sselect_hyperslab_f")
        call h5dread_f(this%dset_id(i), h5t_id, f_ptr, &
                &ierr, file_space_id=fspace_id, mem_space_id=dspace_id, &
                &xfer_prp=plist_id)
        call stop_err(ierr, "h5dread_f")
        ! terminate access to the data space.
        call h5sclose_f(dspace_id, ierr)
        call stop_err(ierr, "h5sclose_f")
        call h5sclose_f(fspace_id, ierr)
        call stop_err(ierr, "h5sclose_f")
    else
        ! open path using the default propertie
        call h5dopen_f(this%f_id(i), path, dsetid, ierr)
        call stop_err(ierr, "h5dopen_f")
        ! read the data using the default properties.
        call h5dread_f(dsetid, h5t_id, f_ptr, &
                &ierr, xfer_prp=plist_id)
        call stop_err(ierr, "h5dread_f")
        ! end access to the dataset and release resources used by it.
        call h5dclose_f(dsetid, ierr)
        call stop_err(ierr, "h5dclose_f")
    endif
    if(.not.serialio)then
        !< close the property list to frees resources.
        call h5pclose_f(plist_id, ierr)
        call stop_err(ierr, "h5pclose_f")
    endif
    return
      
    end subroutine gh5_dread_array


    subroutine gh5_dcreate(this, n, dims, path, h5t_id, i_f, i_d)
    class(hdf5_ob) :: this
    integer, intent(in) :: n
    integer, intent(in) :: dims(n)
    character(len=*), intent(in) :: path
    integer(hid_t), intent(in) :: h5t_id
    integer, intent(in) :: i_f, i_d
   
    integer(hsize_t) :: h_dims(n)
    integer(hid_t) :: dspace_id
    logical :: lexist
    integer :: ierr

    call h5lexists_f(this%f_id(i_f), path, lexist, ierr)
    call stop_err(ierr, "h5lexists_f")
    if(lexist)then
        call h5dopen_f(this%f_id(i_f), path, this%dset_id(i_d), ierr)
        call stop_err(ierr, "h5dopen_f")
    else
        ! create the dataset
        h_dims = dims
        ! create the dataspace.
        call h5screate_simple_f(n, h_dims, dspace_id, ierr)
        call stop_err(ierr, "h5screate_simple_f")
        ! create the dataset with default properties.
        call h5dcreate_f(this%f_id(i_f), path, h5t_id, dspace_id, &
                &this%dset_id(i_d), ierr)
        call stop_err(ierr, "h5dcreate_f")
        ! close dspace
        call h5sclose_f(dspace_id, ierr)
        call stop_err(ierr, "h5sclose_f")
    endif
    return

    end subroutine gh5_dcreate


    subroutine gh5_dopen(this, path, i_f, i_d)
    class(hdf5_ob) :: this
    character(len=*), intent(in) :: path
    integer, intent(in) :: i_d, i_f

    integer :: ierr

    call h5dopen_f(this%f_id(i_f), path, this%dset_id(i_d), ierr)
    call stop_err(ierr, "h5dopen_f")
    return

    end subroutine gh5_dopen


    subroutine gh5_dclose(this, i)
    class(hdf5_ob) :: this
    integer, intent(in) :: i

    integer :: ierr

    call h5dclose_f(this%dset_id(i), ierr)
    call stop_err(ierr, "h5dclose_f")
    this%dset_id(i)=-1
    return

    end subroutine gh5_dclose


    subroutine gh5_gcreate(this, path, i_f, i_g)
    class(hdf5_ob) :: this
    character(len=*), intent(in) :: path
    integer, intent(in) :: i_f
    integer, intent(in), optional :: i_g

    integer(hid_t) :: gid
    integer :: ierr
      
    call h5gcreate_f(this%f_id(i_f), path, gid, ierr)
    call stop_err(ierr, "h5gcreate_f")
    if(present(i_g))then
        this%group_id(i_g)=gid
    else
        call h5gclose_f(gid, ierr)
        call stop_err(ierr, "h5gclose_f")
    endif
    return
      
    end subroutine gh5_gcreate

   
    subroutine gh5_gopen(this, path, i_f, i_g)
    class(hdf5_ob) :: this
    character(len=*), intent(in) :: path
    integer, intent(in) :: i_f, i_g

    integer :: ierr

    call h5gopen_f(this%f_id(i_f), path, this%group_id(i_g), ierr)
    call stop_err(ierr, "h5gopen_f")
    return

    end subroutine gh5_gopen


    subroutine gh5_gclose(this, i)
    class(hdf5_ob) :: this
    integer, intent(in) :: i

    integer :: ierr

    call h5gclose_f(this%group_id(i), ierr)
    call stop_err(ierr, "h5gclose_f")
    this%group_id(i)=-1
    return

    end subroutine gh5_gclose


    subroutine gh5_exists(this, i, path, lexist)
    class(hdf5_ob) :: this
    integer, intent(in) :: i
    character(len=*), intent(in) :: path
    logical, intent(out) :: lexist

    integer :: ierr

    if(path(1:1)/='/')then ! attribute
        call h5aexists_f(this%group_id(i), path, lexist, ierr)
        call stop_err(ierr, "h5aexists_f")
    else ! data set
        call h5lexists_f(this%f_id(i), path, lexist, ierr)
        call stop_err(ierr, "h5lexists_f")
    endif
    return

    end subroutine gh5_exists


    subroutine stop_err(ierr, place)
    character(len=*), intent(in) :: place
    integer, intent(in) :: ierr

    if(ierr/=0)then
        write(0,'("ierr = ", i0, " at ", a, "!")')ierr, place
        call abort()
    endif

    end subroutine stop_err

      
end module ghdf5
