!******************************************************************************
! Author: Yong-Xin Yao <cygutz@gmail.com>
! Function: calculate quasi-particle (qp) band structure
!           and qp local quantities.
!******************************************************************************
program gbands
    use ghdf5
    use bandstru
    use bandsder
    use localstder
    use gmpi
    use gtime
    use dcstd
    use gkernel
    implicit none
   
    type(time_ob) :: time
    type(mpi_ob) :: mpi
    type(hdf5_ob) :: gh5
    type(localstder_ob) :: loc
    type(bandsder_ob) :: bnd
    type(dcstd_ob) :: dc
    type(gkernel_ob) :: gk

    call time%start(1)
    call mpi%init()
    call gh5%init()

    if(mpi%io>0)then
        open(mpi%io,file='Gutz.log',status='unknown',Access='append')
    endif

    call loc%init(gh5,mpi%io)
    call bnd%init(gh5,mpi,loc)
    ! read hk0 with the symmetry-adapted basis
    call bnd%read_hk0(gh5)
    if(bnd%mode_hk==0)then
        write(0,'("error: need to run CyGInit first!")')
        stop
    endif
    call dc%init(loc,gh5,mpi%io)
    call gk%set_params(gh5)
    call gk%load_x(gh5)
    call gk%calc_band(bnd,loc,dc,mpi,time)
    call loc%dump_hembed_list(gh5,gk%iembeddiag)
    if(gk%iwrite>0)then
        call bnd%write_bands(loc,gh5)
    endif
    if(gk%ideri>0)then
        call bnd%calc_derivatives(loc,mpi,mpi%io)
        call loc%h5save_pdlc(gh5)
    endif
    call time%finish(1)
    call time%print_usage('CyGBand',1,mpi%io)
    if(mpi%io>0)then
        close(mpi%io)
    endif
    call gh5%end()
    call mpi%finalize()
      
end program gbands
