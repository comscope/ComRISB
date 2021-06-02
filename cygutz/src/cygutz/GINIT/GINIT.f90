!******************************************************************************
! Author: Yong-Xin Yao <cygutz@gmail.com>
! Function: initial processing, including necessary transformation 
!           of the bare Hamiltonian and append it 
!           to the "GBareH.h5" files, 
!           and generating initial solution vectors.
!******************************************************************************
program ginit
    use ghdf5
    use bandstru
    use localstore
    use dcstd
    use gmpi
    use gkernel
    use gtime
    use gutil, only: rm_file
    implicit none
   
    type(time_ob) :: time
    type(mpi_ob) :: mpi
    type(hdf5_ob) :: gh5
    type(localstore_ob) :: loc
    type(bandstru_ob) :: bnd
    type(dcstd_ob) :: dc
    type(gkernel_ob) :: gk

    ! total execution time
    call time%start(1)
    call mpi%init()
    call gh5%init()

    ! log file
    if(mpi%io>0)then
        open(mpi%io,file='Gutz.log',status='replace')
        !! delete HEmbed.h5 file if existing.
        call rm_file("HEmbed.h5")
    endif

    call loc%init(gh5,mpi%io)
    call bnd%init(gh5,mpi,loc)
    ! read hk0 with the default basis (before symmetry rotations).
    call bnd%read_hk0(gh5)
    ! rotate hk0 to symmetry-adapted basis if mode_hk=0.
    if(bnd%mode_hk==0)then
        call bnd%rotate_hk0(loc)
    else
        ! add back h1e for bare hamiltonian calculation.
        call bnd%rm_h1e_from_hk0(loc,1)
    endif

    ! correlated block energy window.
    call bnd%corr_ebwidth(mpi%io)

    ! check the bare band dispersions.
    call bnd%calc_all(loc,time,mpi)
    call bnd%calc_fermi(mpi%io)
    call bnd%calc_nks(loc,mpi)
    call loc%nks_pp(mpi%io)

    ! double counting
    call dc%init(loc,gh5,mpi%io)

    call gk%set_params(gh5)
    ! initial {R, \lambda} or n (hf) vector
    call loc%init_rln(gh5,mpi%io,gk%r_default,gk%iembeddiag)
    call gk%init_x(loc)
    call gk%dump_x(gh5)

    if(bnd%mode_hk==0)then
        ! single out the local one-body part.
        call bnd%rm_h1e_from_hk0(loc,-1)
        ! update the revised nonlocal bare Hamiltonian
        call bnd%write_hk0_sab(gh5)
    endif
    ! check of symmetrization error
    call loc%symm_check(mpi%io)

    call time%finish(1)
    call time%print_usage('CyGInit',1,mpi%io)
    if(mpi%io>0)then
        close(mpi%io)
    endif
    call gh5%end()
    call mpi%finalize()
      
end program ginit
