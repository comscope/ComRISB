module bandsder
    use gprec
    use bandstru
    use gmpi
    use localstder
    use gconstant, only: zi
    use gutil, only: anmxbmm,calc_loewner,  &
            &fermi_fun2_1derivative,fermi_fun2, &
            &gauss_fun2,gauss_fun2_1derivative, &
            &fermi_fun_1derivative,gauss_fun_1derivative,  &
            &set_mu_delta,annxb
    implicit none
    private

    type,extends(bandstru_ob),public::bandsder_ob
        real(q),pointer :: fw_der(:,:,:)
        real(q) :: sum_fw_der=0
        real(q),allocatable :: pmupr(:)  ! dmu/dr expanded in terms of hs.
        real(q),allocatable :: pmupl(:)  ! dmu/dlambda expanded in terms of hs.
        contains
        procedure::calc_derivatives

    end type

    contains


    subroutine calc_derivatives(this,loc,mpi,io)
    class(bandsder_ob) :: this
    class(localstder_ob) :: loc
    class(mpi_ob) :: mpi
    integer,intent(in)::io

    call loc%init_der()
    call set_mu_delta(this%ef,this%kpt%delta)
    call set_fermi_wt_derivative(this,mpi,io)
    call set_pmupl(this,loc,mpi)
    call set_pmupr(this,loc,mpi)
    call set_pdm(this,loc,mpi)
    call loc%calc_pdm_pp(io)
    call loc%calc_pdplr()
    call loc%calc_plcplr()
    call loc%calc_pdlc_pp(io)
    return

    end subroutine calc_derivatives


    subroutine set_pdm(this,loc,mpi)
    class(bandsder_ob) :: this
    class(localstder_ob) :: loc
    class(mpi_ob) :: mpi

    integer isp,ispp,ikpl,nbase,ib,ikp,isym,nmaxins,nbbase
    integer i,ii,ih,ibase,jbase,dim_hs,na2,na22,na2base,iibase
    real(q) wtk
    complex(q),pointer::vk(:,:),hkr1(:,:)
    real(q),pointer::ek(:)
    real(q)::loemat(this%nmaxin*this%riso,this%nmaxin*this%riso)
    complex(q),target::hkr(this%nmaxin*this%riso,this%nmaxin*this%riso)
    complex(q)::hsp(this%nmaxin*this%riso,this%nmaxin*this%riso),&  ! H's
            &hspl(this%nmaxin*this%riso,this%nmaxin*this%riso)

    nmaxins=this%nmaxin*this%riso
    hkr=0
    if(loc%iso==1)then
        allocate(vk(nmaxins,nmaxins),ek(nmaxins), &
                &hkr1(this%nmaxin,this%nmaxin))
    else
        hkr1=>hkr
    endif
    do ikpl=1,this%kpt%diml
        ikp=this%kpt%base+ikpl
        wtk=this%kpt%wt(ikp)/this%sym%nop
        do isym=1,this%sym%nop
            if(loc%iso==2)then
                nbbase=this%ne(2,ikp,1)-1
                vk=>this%vk(:,:,isym,ikpl,1)
                ek=>this%ek(nbbase+1:nbbase+nmaxins,ikp,1)
            else
                vk=0
                ek=0
            endif
            nbase=0
            do isp=1,this%nspin
                ispp=min(isp,this%nspin_in)
                nbbase=this%ne(2,ikp,ispp)-1
                ! set up hk0 r^\dagger
                hkr1=this%hk0(:,:,isym,ikpl,ispp)
                call anmxbmm('c',hkr1(:,1:loc%nasotot), &
                        &this%r(:,:,isp),this%nmaxin,loc%nasotot) ! hkr^+
                if(loc%iso==1)then
                    hkr(isp::2,isp::2)=hkr1
                    ! padding by zeros should have no effect.
                    ! basis in spin-faster order
                    ! bands still in orbital faster order
                    ek(nbase+1:nbase+this%nmaxin)=this%ek(nbbase+1: &
                            &nbbase+this%nmaxin,ikp,isp)
                    vk(isp::2,nbase+1:nbase+this%nmaxin)= &
                            &this%vk(:,:,isym,ikpl,isp)
                    nbase=nbase+this%nmaxin
                    if(this%nspin==1)then
                        hkr(2::2,2::2)=hkr1
                        ek(nbase+1:)=this%ek(nbbase+1:nbbase+this%nmaxin, &
                                &ikp,isp)
                        vk(2::2,nbase+1:)=this%vk(:,:,isym,ikpl,isp)
                    endif
                endif
            enddo
            ! get loewner matrix
            select case(this%kpt%ismear)
            case(-1)
                call calc_loewner(ek,loemat,nmaxins,fermi_fun2, &
                        &fermi_fun2_1derivative)
            case(0)
                call calc_loewner(ek,loemat,nmaxins,gauss_fun2, &
                        &gauss_fun2_1derivative)
            case default
                stop ' error in set_pdm: &
                        &this%kpt%ismear not supported!'
            end select
            ! loop over hs
            ibase=0
            jbase=0
            na2base=0
            do i=1,loc%num_imp
                na2=loc%na2_list(i)
                if(loc%imap_list(i)/=i)then
                    na2base=na2base+na2
                    cycle
                endif
                na22=na2**2
                dim_hs=loc%hm_l%dim_hs(i)
                do ih=1,dim_hs
                    hsp=0
                    do ii=i,loc%num_imp
                        ! possible repeated blocks
                        if(loc%imap_list(ii)/=i)exit
                        iibase=na2base+(ii-i)*na2
                        hsp(iibase+1:iibase+na2,iibase+1:iibase+na2)= &
                                loc%co(ii)%hs_l(:,:,ih)
                    enddo
                    ! derivative wrt lambda
                    hspl=hsp
                    do ii=1,nmaxins
                        hspl(ii,ii)=hspl(ii,ii)-this%pmupl(ibase+1)
                    enddo
                    call loc%dfah(vk,loemat,hspl,loc%pdmpl(:,ibase+1), &
                            &nmaxins,.true.,wtk)
                    ! p d0 / p lambda
                    ! hspl is now p f(a) / p h
                    call annxb('n','n',hkr,hspl,nmaxins,2)
                    ! part of p d0 / p l
                    call loc%big_to_bkmat(hspl,loc%pd0pl(:,ibase+1), &
                            &nmaxins,.true.,wtk)
                    if(loc%hm_r%dim_hs(i)==dim_hs)then
                        ! non-mott solution
                        call calc_pdmpr(this,loc,hsp,hkr,vk,loemat,nmaxins, &
                                &ibase*loc%r_factor,wtk)
                        call calc_pd0pr2(this,loc,hsp,vk,nmaxins, &
                                &ibase*loc%r_factor,ikp,ikpl,isym,wtk)
                    endif
                    ibase=ibase+1
                enddo
                dim_hs=loc%hm_r%dim_hs(i)
                if(loc%hm_l%dim_hs(i)/=dim_hs)then
                    ! special Mott solution case.
                    do ih=1,dim_hs
                        hsp=0
                        do ii=i,loc%num_imp
                            ! possible repeated blocks
                            if(loc%imap_list(ii)/=i)exit
                            iibase=na2base+(ii-i)*na2
                            hsp(iibase+1:iibase+na2,iibase+1:iibase+na2)= &
                                    &loc%co(ii)%hs_r(:,:,ih)
                        enddo
                        call calc_pdmpr(this,loc,hsp,hkr,vk,loemat,nmaxins, &
                                &jbase,wtk)
                        call calc_pd0pr2(this,loc,hsp,vk,nmaxins,jbase, &
                                &ikp,ikpl,isym,wtk)
                        jbase=jbase+loc%r_factor
                    enddo
                endif
                na2base=na2base+na2
            enddo
        enddo
    enddo
    call mpi%sum_all(loc%pd0pl(:,1),loc%na2112*loc%hm_l%dimhst)
    call mpi%sum_all(loc%pd0pr(:,1),loc%na2112*loc%hm_r%dimhst*loc%r_factor)
    call mpi%sum_all(loc%pdmpl(:,1),loc%na2112*loc%hm_l%dimhst)
    call mpi%sum_all(loc%pdmpr(:,1),loc%na2112*loc%hm_r%dimhst*loc%r_factor)
   return

    end subroutine set_pdm


    subroutine calc_pdmpr(this,loc,h,hkr,vk,loemat,n,ibase,wt)
    class(bandsder_ob) :: this
    class(localstder_ob) :: loc
    integer,intent(in)::n,ibase
    complex(q),intent(in)::hkr(n,n),vk(n,n),h(n,n)
    real(q),intent(in)::loemat(n,n),wt

    integer i
    complex(q) hp(n,n),hhkr(n,n)

    ! derivative wrt r
    ! hs \epsilon_k r^\dagger
    hhkr=h
    call annxb('n','n',hhkr,hkr,n,1)
    hp=hhkr+conjg(transpose(hhkr))
    do i=1,n
        hp(i,i)=hp(i,i)-this%pmupr(ibase+1)
    enddo
    call loc%dfah(vk,loemat,hp,loc%pdmpr(:,ibase+1),n,.true.,wt)
    ! now hp is p dm / p r
    call annxb('n','n',hkr,hp,n,2)
    ! part of p d0 / p r
    call loc%big_to_bkmat(hp,loc%pd0pr(:,ibase+1),n,.true.,wt)
    if(loc%r_factor==2)then  ! complex-r version, imaginary part
        hp=(hhkr-conjg(transpose(hhkr)))*zi
        do i=1,n
            hp(i,i)=hp(i,i)-this%pmupr(ibase+2)
        enddo
        call loc%dfah(vk,loemat,hp,loc%pdmpr(:,ibase+2),n, &
                &.true.,wt)
        call annxb('n','n',hkr,hp,n,2)
        call loc%big_to_bkmat(hp,loc%pd0pr(:,ibase+2),n, &
                &.true.,wt)
    endif
    return

    end subroutine calc_pdmpr


    subroutine calc_pd0pr2(this,loc,h,vk,n,ibase,ikp,ikpl,isym,wt)
    class(bandsder_ob) :: this
    class(localstder_ob) :: loc
    integer,intent(in)::n,ibase,ikp,ikpl,isym
    complex(q),intent(in)::h(n,n),vk(n,n)
    real(q),intent(in)::wt

    integer i,isp,ispp,ibbase,ib,nemin,nemax,nbands,rso
    complex(q) fh(n,n),hk(n,n),hkh(n,n)
 
    fh=0
    ibbase=0
    hk=0
    rso=2/loc%iso
    ! f(H): fermi weight in orbital basis
    do i=1,rso
        isp=min(i,this%nspin)
        ispp=min(isp,this%nspin_in)
        nemin=this%ne(2,ikp,ispp)-1
        nemax=this%ne(3,ikp,ispp)
        nbands=nemax-nemin
        do ib=1,nbands
            fh(:,ibbase+ib)=vk(:,ibbase+ib)*this%ferwe(nemin+ib,ikp,isp) &
                    &/this%rspo/this%kpt%wt(ikp)
        enddo
        ! hk0
        hk(i::rso,i::rso)=this%hk0(:,:,isym,ikpl,ispp)
        ibbase=ibbase+this%nmaxin
    enddo
    call annxb('n','c',fh,vk,n,1)
    ! real part
    hkh=hk
    ! hk0 * h
    call annxb('n','n',hkh,h,n,1)
    ! hk0 * h * f(H)
    call annxb('n','n',hkh,fh,n,1)
    call loc%big_to_bkmat(hkh,loc%pd0pr(:,ibase+1),n, &
            &.true.,wt)
    if(loc%r_factor==2)then  ! complex-r version, imaginary part
        hkh=h*zi
        call annxb('n','n',hk,hkh,n,2)
        ! hk0 * ih * f(H)
        call annxb('n','n',hkh,fh,n,1)
        call loc%big_to_bkmat(hk,loc%pd0pr(:,ibase+2),n, &
                &.true.,wt)
    endif
    return

    end subroutine calc_pd0pr2


    subroutine set_pmupl(this,loc,mpi)
    ! derivative of chemical potential mu wrt ls of lambda coef:
    ! p mu/ p ls = Tr[f'(A)hs] / Tr[f'(A)]
    ! we will use subtoutines of calculating nks for the evaluation,
    ! because they follow the same math.
    ! should be called before calculating nks, otherwise got overwritten.
    class(bandsder_ob) :: this
    class(localstder_ob) :: loc
    class(mpi_ob) :: mpi

    complex(q)::nks_buf(loc%na2112),nks_coef_buf(loc%hm_l%dimhst)

    allocate(this%pmupl(loc%hm_l%dimhst))
    this%pmupl=0
    if(this%kpt%ensemble/=0)then
        ! grand canonical ensemble.
        return
    endif
    ! save the nks data.
    nks_buf=loc%nks
    nks_coef_buf=loc%nks_coef
    call this%calc_nks(loc,mpi,this%fw_der)
    call loc%nks_pp(mpi%io,obj="nkp")
    this%pmupl=loc%nks_coef/this%sum_fw_der
    ! adding degeneracy factors.
    call add_pmupl_degenaracies(this,loc)
    ! recover nks data.
    loc%nks=nks_buf
    loc%nks_coef=nks_coef_buf
    return

    end subroutine set_pmupl


    subroutine add_pmupl_degenaracies(this,loc)
    class(bandsder_ob) :: this
    class(localstder_ob) :: loc

    integer i,imap,hsbase,deg,dim_hs

    hsbase=0
    do i=1,loc%num_imp
        imap=loc%imap_list(i)
        if(imap/=i)then
            cycle
        endif
        deg=count(loc%imap_list==i)
        dim_hs=loc%hm_l%dim_hs(i)
        this%pmupl(hsbase+1:hsbase+dim_hs)= &
                &this%pmupl(hsbase+1:hsbase+dim_hs)*deg
        hsbase=hsbase+dim_hs
    enddo
    return

    end subroutine add_pmupl_degenaracies


    subroutine set_pmupr(this,loc,mpi)
    ! derivative of chemical potential mu wrt ls of r coef:
    ! p mu/ p rs_real = Tr[f'(A)(hs hk0 r^dagger+h.c.)] / Tr[f'(A)]
    ! p mu/ p rs_imag = Tr[f'(A)i(hs hk0 r^dagger-h.c.)] / Tr[f'(A)]
    class(bandsder_ob) :: this
    class(localstder_ob) :: loc
    class(mpi_ob) :: mpi

    complex(q),target::h0rv(loc%na2tot,this%nmaxin*this%riso)
    complex(q),pointer::h0rv1(:,:)
    ! merge both spin channels, to be consistent with hs convention.
    complex(q) vk(loc%na2tot,this%nmaxin*this%riso)
    integer isp,ispp,ikpl,nbase,ib,ikp,isym,nmaxins
    real(q) wtk
    complex(q) zes
    real(q),parameter::tol=1.e-6_q
    integer i,ih,ibase,hsbase,dim_hs,na2,na22,na2base
    complex(q),pointer::p_hs(:,:),p_vl(:,:),p_vr(:,:)
    complex(q),target::zbuf1(loc%na2tot*this%nmaxin*this%riso), &
            &zbuf2(loc%na2tot*this%nmaxin*this%riso), &
            &zbuf3(loc%na2max*loc%na2max)
    complex(q),external::zdotc

    allocate(this%pmupr(loc%hm_r%dimhst*loc%r_factor))
    this%pmupr=0
    if(this%kpt%ensemble/=0)then
        ! grand canonical ensemble.
        return
    endif
    if(loc%iso==2)then
        h0rv1=>h0rv
    else
        allocate(h0rv1(loc%nasotot,this%nmaxin))
    endif
    nmaxins=this%nmaxin*this%riso
    wtk=1._q/this%sym%nop/this%rspo
    do ikpl=1,this%kpt%diml
        ikp=this%kpt%base+ikpl
        do isym=1,this%sym%nop
            h0rv=0
            do isp=1,this%nspin
                ispp=min(isp,this%nspin_in)
                ! hk0
                h0rv1=this%hk0(1:loc%nasotot,:,isym,ikpl,ispp)
                ! hk0*r^\dagger
                call anmxbmm('c',h0rv1(:,1:loc%nasotot),this%r(:,:,isp), &
                        &loc%nasotot,loc%nasotot)
                ! hk0*r^\dagger*v
                call anmxbmm('n',h0rv1,this%vk(:,:,isym,ikpl,isp),&
                        &loc%nasotot,this%nmaxin)
                nbase=this%ne(2,ikp,ispp)-1
                do ib=1,this%nmaxin
                    h0rv1(:,ib)=h0rv1(:,ib)*this%fw_der(nbase+ib,ikp,isp)*wtk
                enddo
                if(loc%iso==1)then
                    h0rv(isp::2,:this%nmaxin)=h0rv1
                    vk(isp::2,:this%nmaxin)=this%vk(:loc%nasotot,:, &
                            &isym,ikpl,isp)
                    if(this%nspin==1)then
                        h0rv(2::2,this%nmaxin+1:)=h0rv1
                        vk(2::2,this%nmaxin+1:)=this%vk(:loc%nasotot,:, &
                                &isym,ikpl,isp)
                    endif
                else
                    vk=this%vk(:loc%nasotot,:,isym,ikpl,isp)
                endif
            enddo  ! isp

            hsbase=0
            ibase=0
            na2base=0
            do i=1,loc%num_imp
                na2=loc%na2_list(i)
                na22=na2**2
                dim_hs=loc%hm_r%dim_hs(i)
                if(dim_hs<=0)then
                    na2base=na2base+na2
                    cycle
                endif
                p_vl(1:na2,1:nmaxins)=>zbuf1(:na2*nmaxins)
                p_vl=vk(na2base+1:na2base+na2,:)
                p_vr(1:na2,1:nmaxins)=>zbuf2(:na2*nmaxins)
                p_hs(1:na2,1:na2)=>zbuf3(1:na22)
                do ih=1,dim_hs
                    p_vr=h0rv(na2base+1:na2base+na2,:)
                    zbuf3(1:na22)=loc%hm_r%hs(hsbase+1:hsbase+na22)
                    call annxb('n',p_hs,p_vr,na2,nmaxins)
                    zes=zdotc(na2*nmaxins,p_vl(1,1),1,p_vr(1,1),1)
                    this%pmupr(ibase+(ih-1)*loc%r_factor+1)= &
                            &this%pmupr(ibase+(ih-1)*loc%r_factor+1)+ &
                            &2*real(zes)
                    if(loc%r_factor==2)then
                        this%pmupr(ibase+ih*2)=this%pmupr(ibase+ih*2)- &
                                &2*aimag(zes)
                    endif
                    hsbase=hsbase+na22
                enddo
                if(i<loc%num_imp)then
                    ! a new inequivalent site to be expected 
                    if(loc%imap_list(i+1)==i+1)then
                        ibase=ibase+dim_hs*loc%r_factor
                    endif
                endif
                na2base=na2base+na2
            enddo
        enddo
    enddo
    call mpi%sum_all(this%pmupr,loc%hm_r%dimhst*loc%r_factor)
    this%pmupr=this%pmupr/this%sum_fw_der
    return

    end subroutine set_pmupr


    subroutine set_fermi_wt_derivative(this,mpi,io)
    ! calculate the array of fermi weight derivatives
    class(bandsder_ob) :: this
    class(mpi_ob) :: mpi
    integer,intent(in)::io

    integer isp,ispp,ikp,ib
    real(q) dt

    ! same format as ferwe
    allocate(this%fw_der(this%nmax,this%kpt%dim,this%nspin))
    this%fw_der=0
    do isp=1,this%nspin
        ispp=min(isp,this%nspin_in)
        do ikp=1,this%kpt%dim
            do ib=1,this%ne(1,ikp,ispp)
                dt=(this%ek(ib,ikp,isp)-this%ef)/this%kpt%delta
                select case(this%kpt%ismear)
                case(-1)
                    this%fw_der(ib,ikp,isp)=fermi_fun_1derivative(dt)
                case(0)
                    this%fw_der(ib,ikp,isp)=gauss_fun_1derivative(dt)
                case default
                    stop ' error in set_fermi_wt_derivative: &
                            &this%kpt%ismear not supported!'
                end select
                this%fw_der(ib,ikp,isp)=this%fw_der(ib,ikp,isp) &
                        &*this%rspo*this%kpt%wt(ikp)
            enddo
        enddo
    enddo
    this%sum_fw_der=sum(this%fw_der)
    if(io>0)then
        write(io,'(" (-)sum of fermi weight derivatives", f10.4)') &
                &-this%sum_fw_der
    endif
    return

    end subroutine set_fermi_wt_derivative



end module bandsder
