!******************************************************************************
! Author: Yong-Xin Yao <cygutz@gmail.com>
! Function: calculate vector function error of Gutzwiller
!           nonlinear equations.
!******************************************************************************
program gerr
    use ghdf5
    use bandsder
    use localstder
    use gmpi
    use gtime
    use gkerder
    implicit none
   
    type(time_ob) :: time
    type(mpi_ob) :: mpi
    type(hdf5_ob) :: gh5
    type(localstder_ob) :: loc
    type(bandsder_ob) :: bnd
    type(gkerder_ob) :: gk

    call time%start(1)
    call mpi%init()
    call gh5%init()

    if(mpi%io>0)then
        open(mpi%io,file='Gutz.log',status='unknown',Access='append')
    endif

    call loc%init(gh5,mpi%io)
    call gk%set_params(gh5)
    call gk%load_x(gh5)
    call loc%load_hembed_list(gh5,mpi%io,gk%iembeddiag)
    call loc%nks_pp(mpi%io)
    call loc%ncvar_pp(mpi%io)
    if(gk%iembeddiag/=10)then
        ! not hartree-fock calculation.
        call loc%calc_isimix(mpi%io,mode=1)
        call loc%r01_pp(mpi%io)
    endif
    call gk%save_verr(gh5,loc,mpi%io)

    if(gk%ideri>0)then
        call loc%setup2_pnrplr(gh5,mpi%io)
        call gk%save_jacobian(gh5,loc,mpi%io)
    endif
    call time%finish(1)
    call time%print_usage('CyGErr',1,mpi%io)
    if(mpi%io>0)then
        close(mpi%io)
    endif
    call gh5%end()
    call mpi%finalize()
      
end program gerr
