module ghdf5_sm
    use gprec
    use ghdf5
    use sparse
    use hdf5
    implicit none
    private

    type,extends(hdf5_ob),public::hdf5_sm_ob
        contains
        procedure::write_sp_matrix=>gh5_write_sp_matrix
    end type

    contains

    subroutine gh5_write_sp_matrix(this, spmat, path, f_id, serialio)
    class(hdf5_sm_ob)::this
    class(sp_matrix), target, intent(in) :: spmat
    character(len=*), intent(in) :: path
    integer, intent(in) :: f_id
    logical, optional, intent(in) :: serialio

    logical :: serial_io=.false.
    integer nnz

    if(present(serialio))then
        serial_io=serialio
    endif
    call this%gcreate(path, f_id, 1)
    call this%awrite(spmat%nrow, "nrow", 1)
    call this%awrite(spmat%ncol, "ncol", 1)
    call this%awrite(1, "base", 1)
    select type(spmat)
    type is(dcsr_matrix)
        call this%dwrite(spmat%i, spmat%nrow + 1, path//"/indptr", f_id, &
                &serialio=serial_io)
        nnz = spmat%i(spmat%nrow + 1) -1
        call this%dwrite(spmat%j, nnz, path//"/indices", f_id, &
                &serialio=serial_io)
        call this%dwrite(spmat%a, nnz, path//"/data", f_id, &
                &serialio=serial_io)
    type is(zcsr_matrix)
        call this%dwrite(spmat%i, spmat%nrow + 1, path//"/indptr", f_id, &
                &serialio=serial_io)
        nnz = spmat%i(spmat%nrow + 1) -1
        call this%dwrite(spmat%j, nnz, path//"/indices", f_id, &
                &serialio=serial_io)
        call this%dwrite(spmat%a, nnz, path//"/data", f_id, &
                &serialio=serial_io)
    type is(dcoo_matrix)
        call this%awrite(spmat%nnz, "nnz", 1)
        call this%dwrite(spmat%i, spmat%nnz, path//"/i", f_id, &
                &serialio=serial_io)
        call this%dwrite(spmat%j, spmat%nnz, path//"/j", f_id, &
                &serialio=serial_io)
        call this%dwrite(spmat%a, spmat%nnz, path//"/data", f_id, &
                &serialio=serial_io)
    end select
    call this%gclose(1)
    return

    end subroutine gh5_write_sp_matrix


end module ghdf5_sm
