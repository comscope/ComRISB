!******************************************************************************
! Author: Yong-Xin Yao <cygutz@gmail.com>
! Function: calculate once more quasi-particle (qp) band structure,
!           total energy of the model, and optionally the renormalized
!           electron density.
!******************************************************************************
program gfinal
    use ghdf5
    use bandstru
    use localstore
    use gmpi
    use gtime
    use dcstd
    use gkernel
    implicit none
   
    type(time_ob) :: time
    type(mpi_ob) :: mpi
    type(hdf5_ob) :: gh5
    type(localstore_ob) :: loc
    type(bandstru_ob) :: bnd
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
    else
        ! add back h1e for bare hamiltonian calculation.
        call bnd%rm_h1e_from_hk0(loc,1)
    endif
    ! check the bare band dispersions.
    call bnd%calc_all(loc,time,mpi)
    call bnd%calc_fermi(mpi%io)
    call bnd%calc_bz_integration_correction()
    call bnd%calc_ebnd(mpi, bare_mode=1)
    ! single out the local one-body part.
    call bnd%rm_h1e_from_hk0(loc,-1)

    call dc%init(loc,gh5,mpi%io)
    call gk%set_params(gh5)
    call gk%load_x(gh5)
    call gk%calc_band(bnd,loc,dc,mpi,time)
    call loc%load_hembed_list(gh5,mpi%io,gk%iembeddiag)
    call loc%ncvar_pp(mpi%io)
    call loc%ncphy_pp(mpi%io)

    if(gk%iembeddiag/=10)then
        ! not hartree-fock calculation.
        call loc%calc_isimix(mpi%io,mode=1)
        call loc%r01_pp(mpi%io)
    endif
    call gk%calc_total_energy(bnd,loc,dc,mpi)

    if(gk%iupdaterho>0)then
        call bnd%calc_rnrl(loc)
        call loc%map_bnd_matrix(loc%nc_phy,bnd%nc_phy,.false.)
        if(gk%iupdaterho==1)then
            call bnd%calc_kswt(loc,mpi,time)
            call bnd%write_kswt(loc,gh5)
        endif
    endif

    call dc%chk_err(loc,mpi%io)
    call bnd%calc_num_electrons(mpi%io)
    call gk%save_results(bnd,loc,dc,gh5)
    call loc%symm_check(mpi%io)

    call time%finish(1)
    call time%print_usage('CyGFinal',1,mpi%io)
    if(mpi%io>0)then
        close(mpi%io)
    endif
    call gh5%end()
    call mpi%finalize()
      
end program gfinal
