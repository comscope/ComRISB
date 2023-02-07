#ifdef constraint_s2
subroutine set_s2_spci()
    use gprec
    use gspci
    use gutil
    implicit none
    integer irun,ival,istate,i1,i2,jstate,nnz,nhalf, &
            &itmp,mask1,ibs1,ibs2,ia1,ia2,inz0, &
            &itmp1,itmp2,iup,idn
    real(q) sz1,sz2

    nhalf=dmem%norb/2
    dmem%s2op%nrow=dmem%nstates
    dmem%s2op%ncol=dmem%nstates
    allocate(dmem%s2op%i(dmem%nstates+1)); dmem%s2op%i(1)=1
    mask1=ishft(1,dmem%norb)-1
    do irun=1,2
        if(irun==2)then
            allocate(dmem%s2op%j(nnz),dmem%s2op%a(nnz))
        endif
        nnz=0
        do ival=dmem%nval_bot,dmem%nval_top
            do i1=dmem%idx(ival),dmem%idx(ival+1)-1
                if(irun==2)then
                    sz1=dmem%bs_sz(i1)
                endif
                do i2=dmem%idx_l(dmem%norb-ival), &
                        &dmem%idx_l(dmem%norb-ival+1)-1
                    istate=dmem%i_phi(i1)%i(i2)
                    if(istate==0)then
                        cycle
                    endif
                    if(irun==2)then
                        sz2=dmem%bs_sz_l(i2)
                    endif
                    nnz=nnz+1 ! part sz^2 + sz
                    if(irun==2)then
                        inz0=nnz
                        dmem%s2op%j(nnz)=istate
                        dmem%s2op%a(nnz)=(sz1+sz2)**2+sz1+sz2
                    endif
                    itmp=ior(dmem%bs(i1),ishft(dmem%bs_l(i2),dmem%norb))
                    do ia1=1,dmem%norb
                        itmp1=itmp
                        ! s_ia1^+
                        ! spin-dn orbital
                        idn=2*ia1-2+dmem%idn
                        if(btest(itmp1,idn))then
                            itmp1=ibclr(itmp1,idn)
                        else
                            cycle
                        endif
                        ! spin-up orbital
                        iup=2*ia1-2+dmem%iup
                        if(btest(itmp1,iup))then
                            cycle
                        else
                            itmp1=ibset(itmp1,iup)
                        endif
                        do ia2=1,dmem%norb
                            ! s_ia2^- s_ia1^+
                            itmp2=itmp1
                            if(ia2==ia1)then
                                if(irun==2)then
                                    dmem%s2op%a(inz0)=dmem%s2op%a(inz0)+1._q
                                endif
                                cycle
                            ! Only keep lower trianglar part
                            elseif(ia1<=nhalf.and.ia2>ia1.and.ia2<=nhalf)then
                                cycle
                            elseif(ia1>nhalf.and.ia2<=nhalf)then
                                cycle
                            elseif(ia1>nhalf.and.ia2>ia1.and.ia2>nhalf)then
                                cycle
                            else
                                iup=2*ia2-2+dmem%iup
                                if(btest(itmp2,iup))then
                                    itmp2=ibclr(itmp2,iup)
                                else
                                    cycle
                                endif
                                idn=2*ia2-2+dmem%idn
                                if(btest(itmp2,idn))then
                                    cycle
                                else
                                    itmp2=ibset(itmp2,idn)
                                endif
                            endif
                            nnz=nnz+1
                            if(irun==1)cycle
                            ibs1=dmem%ibs(iand(itmp2,mask1))
                            ibs2=dmem%ibs_l(ishft(itmp2,-dmem%norb))
                            jstate=dmem%i_phi(ibs1)%i(ibs2)
                            dmem%s2op%j(nnz)=jstate
                            dmem%s2op%a(nnz)=1._q
                        enddo
                    enddo
                    dmem%s2op%i(istate+1)=nnz+1
                enddo
            enddo
        enddo
    enddo
    return

end subroutine set_s2_spci


subroutine chk_eval_s2(lstop)
    use gprec
    use gspci
    use sparse
    implicit none
    logical,intent(in)::lstop

#ifdef real_version
    real(q) eval
#else
    complex(q) eval
#endif

    eval=0
    call vh_sycsr_v(dmem%s2op,dmem%v,eval)
    write(0,*) "<S^2> = ", eval
    if(lstop)then
        if(abs(eval)>1.d-8)then
            write(0,'(" Warning: <S^2> is not zero!")')
        endif
    endif
    return

end subroutine chk_eval_s2


subroutine av1_gspci_s2(v1,v2)
    use gprec
    use gspci
    use sparse
    implicit none
#ifdef real_version
    real(q),  &
#else
    complex(q),  &
#endif
            &intent(in)::v1(*)
#ifdef real_version
    real(q),  &
#else
    complex(q),  &
#endif
            &intent(inout)::v2(*)

    ! s2op
    if(abs(dmem%lambda_j2)>1.d-12)then
        call csr_syamux(dmem%lambda_j2,dmem%s2op,v1,v2)
    endif
    return

end subroutine av1_gspci_s2


subroutine add_hdns_spci_s2(a,n)
    use gprec
    use gspci
    implicit none
    integer,intent(in)::n
#ifdef real_version
    real(q),  &
#else
    complex(q),  &
#endif
            &intent(inout)::a(n,n)

    integer i,inz,j

    ! s2, here we search for s=0 solution.
    if(abs(dmem%lambda_j2)>1.d-12)then
        do i=1,dmem%s2op%nrow
            do inz=dmem%s2op%i(i),dmem%s2op%i(i+1)-1
                j=dmem%s2op%j(inz)
                a(i,j)=a(i,j)+dmem%s2op%a(inz)*dmem%lambda_j2
                if(i==j)cycle
                a(j,i)=a(j,i)+dmem%s2op%a(inz)*dmem%lambda_j2
            enddo
        enddo
    endif
    return

end subroutine add_hdns_spci_s2
#endif constraint_s2


#ifdef nc_csr
! In physical subspace.
subroutine setup_npcoo_op1()
    use gprec
    use gspci
    use sparse
    use gutil
    implicit none
    integer nstates,ia1,ia2,isum,nbase,ival,itmp1,isgn,i1,i2,ibs

    allocate(dmem%ncabop(dmem%norb,dmem%norb))

    nstates=dmem%idx(dmem%nval_top)-dmem%idx(dmem%nval_bot)
    do ia1=1,dmem%norb; do ia2=1,ia1
        call alloc_sp_matrix(dmem%npcoo(ia1,ia2),nstates,nstates,nstates)
        nbase=dmem%idx(dmem%nval_bot)-1
        isum=0
        ! < state | c_ia1^\dagger c_ia2
        do ival=dmem%nval_bot,dmem%nval_top
            do i1=dmem%idx(ival),dmem%idx(ival+1)-1
                itmp1=dmem%bs(i1)
                isgn=1 
                call act_state(itmp1,ia1-1,.false.,isgn)
                if(isgn==0)cycle
                call act_state(itmp1,ia2-1,.true.,isgn)
                if(isgn==0)cycle
                ibs=dmem%ibs(itmp1)
                isum=isum+1
                dmem%npcoo(ia1,ia2)%i(isum)=i1-nbase
                dmem%npcoo(ia1,ia2)%j(isum)=ibs-nbase
                dmem%npcoo(ia1,ia2)%a(isum)=real(isgn,q)
            enddo
        enddo
        dmem%npcoo(ia1,ia2)%nnz=isum
    enddo; enddo
    return

end subroutine setup_ncab_op
#endif nc_csr


subroutine modify_m_struct(svec,lvec)
    use gprec
    use gspci
    implicit none
    complex(q),intent(in)::svec(dmem%norb,dmem%norb,3), &
            &lvec(dmem%norb,dmem%norb,3)

    integer i,j

    do i=1,dmem%norb; do j=1,dmem%norb
        if(dmem%m_struct(i,j)>0)cycle
        if(sum(abs(svec(i,j,:)))>1.d-12.or. &
                &sum(abs(svec(j,i,:)))>1.d-12.or. &
                &sum(abs(lvec(i,j,:)))>1.d-12.or. &
                &sum(abs(lvec(j,i,:)))>1.d-12)then
            dmem%m_struct(i,j)=101
        endif
    enddo; enddo
    return

end subroutine modify_m_struct


subroutine set_full_fock_states_l_spci()
    use gprec
    use gspci
    implicit none
    integer i

    if(dmem%norb_mott>0)then
        call set_dmem_mott()
    else
        dmem%bs_l=>dmem%bs; dmem%ibs_l=>dmem%ibs
        dmem%idx_l=>dmem%idx; dmem%bs_sz_l=>dmem%bs_sz
    endif
    call calc_dimphi_sz()
    return

end subroutine set_full_fock_states_l_spci


subroutine calc_dimphi_sz()
    use gprec
    use gspci
    use gutil
    implicit none

    integer ival,i,j,isum
    real(q) sz

    allocate(dmem%i_phi(dmem%idx(dmem%nval_bot): &
            &dmem%idx(dmem%nval_top+1)-1))
    isum=0
    do ival=dmem%nval_bot,dmem%nval_top
        do i=dmem%idx(ival),dmem%idx(ival+1)-1
            allocate(dmem%i_phi(i)%i(dmem%idx_l(dmem%norb-ival): &
                    &dmem%idx_l(dmem%norb-ival+1)-1))
            dmem%i_phi(i)%i=0
            sz=dmem%bs_sz(i)
            do j=dmem%idx_l(dmem%norb-ival),dmem%idx_l(dmem%norb-ival+1)-1
                if(dmem%mode==201)then
                    ! diagonal fock solver
                    if(dmem%bs(i)+dmem%bs_l(j)+1==2**dmem%norb)then
                        isum=isum+1
                        dmem%i_phi(i)%i(j)=isum
                        dmem%i_phi(i)%imin=isum
                        dmem%i_phi(i)%imax=isum
                        exit
                    endif
                else
                    if(abs(sz+dmem%bs_sz_l(j))<1.d-6)then
                        isum=isum+1
                        dmem%i_phi(i)%i(j)=isum
                        if(dmem%i_phi(i)%imin==0)then
                            dmem%i_phi(i)%imin=isum
                        endif
                        dmem%i_phi(i)%imax=isum
                    endif
                endif
            enddo
        enddo
    enddo
    dmem%nstates=isum
    write(0,'(" dim_phi = ", i0)')dmem%nstates
    return


end subroutine calc_dimphi_sz


subroutine set_dmem_mott()
    use gprec
    use gspci
    use gutil
    implicit none
    integer i,n,nbs,isum,sume

    nbs=ishft(1,dmem%norb)
    allocate(dmem%idx_l(0:dmem%norb+1),dmem%bs_l(nbs),dmem%ibs_l(0:nbs-1), &
            &dmem%bs_sz_l(nbs))
    dmem%ibs_l=0
    dmem%idx_l(0)=1; isum=1
    do n=0,dmem%norb
        if(n<=dmem%norb-dmem%nelect_mott)then
            do i=dmem%idx(n),dmem%idx(n+1)-1
                call sum_empty_slots_fs(dmem%bs(i),dmem%norb_mott, &
                        &dmem%iorb_mott,sume)
                if(sume==dmem%nelect_mott)then
                    dmem%bs_l(isum)=dmem%bs(i)
                    dmem%bs_sz_l(isum)=dmem%bs_sz(i)
                    dmem%ibs_l(dmem%bs(i))=isum
                    isum=isum+1
                endif
            enddo
        endif
        dmem%idx_l(n+1)=isum
    enddo
    return

end subroutine set_dmem_mott


subroutine solve_hembed_spci_drive()
    use gprec
    use gspci
    implicit none
    integer method
    real(q) w(2)

    if(dmem%norb2<=12)then
        method=0
    else
        method=1
    endif

    if(.not.associated(dmem%v))then
        allocate(dmem%v(dmem%nstates))
        w(1)=0
    else
        w(1)=1
    endif
    call solve_hembed_spci(w(1:2),method)
    if(method==0)then
        write(0,'(" two lowest eigen-values: ", 2f14.6)')w(1:2)
    endif
    dmem%etot=w(1)
    return

end subroutine solve_hembed_spci_drive


subroutine solve_hembed_spci(w,method)
    use gprec
    use gspci
    use gutil
    implicit none
    real(q),intent(inout)::w(2)
    integer method
    external::av_gspci

    if(method==0)then
        call lapack_diag_spci(dmem%v,dmem%nstates,w)
    else

#ifdef debug_mode
        call zprimme_diag_chk(dmem%nstates,10,av_gspci)
#endif
        call primme_diag(dmem%v,dmem%nstates,w,av_gspci)
    endif
    return

end subroutine solve_hembed_spci


subroutine solve_hembed_spci_chk()
    use gspci
    use gutil
    implicit none

    integer nchk
    external::av_gspci

    nchk=min(10,dmem%nstates-1)
    call zprimme_diag_chk(dmem%nstates,nchk,av_gspci)

end subroutine solve_hembed_spci_chk


subroutine lapack_diag_spci(v,n,w)
    use gprec
    use gspci
    use gutil
    implicit none
    integer,intent(in) :: n
#ifdef real_version
    real(q),  &
#else
    complex(q),  &
#endif
            &intent(out) :: v(n)

    real(q),intent(out) :: w(2)

    integer i
    real(q) w_(n)
#ifdef real_version
    real(q) z(n,n),v1(n)
#else
    complex(q) z(n,n),v1(n)
#endif

    z=0
    call setup_hdns_spci(z,n)
    call hermev('v','l',z,w_,n)
    v=z(:,1)
    w=w_(:2)
    return

end subroutine lapack_diag_spci


subroutine set_mncoo_spci()
    use gprec
    use gspci
    implicit none
    integer i,j

    allocate(dmem%mcoo(dmem%norb,dmem%norb), &
            &dmem%nvcoo(dmem%norb,dmem%norb), &
            &dmem%npcoo(dmem%norb,dmem%norb))
    do i=1,dmem%norb
        do j=1,dmem%norb

            if(dmem%m_struct(i,j)<=0)cycle
            ! c_j^\dagger f_i, note act to left
            dmem%mcoo(j,i)%nrow=dmem%nstates
            dmem%mcoo(j,i)%ncol=dmem%nstates
            call set_dcoo_ij_spci(dmem%mcoo(j,i),j,i+dmem%norb,.false.,.true.)

            if(i<j)cycle ! hermitian

            ! f_j f_i^\dagger, note act to left
            dmem%nvcoo(j,i)%nrow=dmem%nstates
            dmem%nvcoo(j,i)%ncol=dmem%nstates
            call set_dcoo_ij_spci(dmem%nvcoo(j,i),j+dmem%norb,i+dmem%norb, &
                    &.true.,.false.)

            ! c_i^\dagger c_j  
            dmem%npcoo(i,j)%nrow=dmem%nstates
            dmem%npcoo(i,j)%ncol=dmem%nstates
            call set_dcoo_ij_spci(dmem%npcoo(i,j),i,j,.false.,.true.)
        enddo
    enddo
    return

end subroutine set_mncoo_spci


subroutine set_dcoo_ij_spci(dcoo,icf,jcf,iact,jact)
    use gprec
    use gspci
    use sparse
    use gutil
    implicit none
    type(dcoo_matrix),intent(inout)::dcoo
    logical,intent(in)::iact,jact
    integer,intent(in)::icf,jcf

    integer irun,ival,istate,i1,i2,jstate,nnz,dval, &
            &itmp,isgn,mask1,ibs1,ibs2,nfs,nfs_l

    dval=0
    if(icf<=dmem%norb)then
        if(iact)then
            dval=dval+1
        else
            dval=dval-1
        endif
    endif
    if(jcf<=dmem%norb)then
        if(jact)then
            dval=dval+1
        else
            dval=dval-1
        endif
    endif
    mask1=ishft(1,dmem%norb)-1
    do irun=1,2
        if(irun==2)then
            dcoo%nnz=nnz
            allocate(dcoo%i(nnz),dcoo%j(nnz),dcoo%a(nnz))
        endif
        nnz=0
        do ival=dmem%nval_bot,dmem%nval_top
            if(dval==1)then
                if(ival==dmem%nval_top)then
                    exit
                endif
            elseif(dval==-1.and.ival==dmem%nval_bot)then
                cycle
            elseif(abs(dval)>1)then
                write(0,'(" Error in set_dcoo_ij_spci!")')
                stop 2
            endif
            nfs_l=dmem%idx_l(dmem%norb-(ival+dval)+1)- &
                    &dmem%idx_l(dmem%norb-(ival+dval))
            if(nfs_l<=0)then
                cycle
            endif
            do i1=dmem%idx(ival),dmem%idx(ival+1)-1
                do i2=dmem%idx_l(dmem%norb-ival), &
                        &dmem%idx_l(dmem%norb-ival+1)-1
                    istate=dmem%i_phi(i1)%i(i2)
                    if(istate==0)then
                        cycle
                    endif
                    itmp=ior(dmem%bs(i1),ishft(dmem%bs_l(i2),dmem%norb))
                    isgn=1
                    call act_state(itmp,icf-1,iact,isgn)
                    if(isgn==0)cycle
                    call act_state(itmp,jcf-1,jact,isgn)
                    if(isgn==0)cycle
                    ibs2=dmem%ibs_l(ishft(itmp,-dmem%norb))
                    if(ibs2<=0)cycle
                    nnz=nnz+1
                    if(irun==1)cycle
                    ibs1=dmem%ibs(iand(itmp,mask1))
                    jstate=dmem%i_phi(ibs1)%i(ibs2)
#ifdef DEBUG
                    if(jstate<=0)then
                        write(0,'(" Error in set_dcoo_ij_spci with &
                                &jstate = ", i0)')jstate
                        stop 4
                    endif
#endif
                    dcoo%i(nnz)=istate
                    dcoo%j(nnz)=jstate
                    dcoo%a(nnz)=isgn
                enddo
            enddo
        enddo
    enddo
    return

end subroutine set_dcoo_ij_spci


subroutine calc_ucsr_spci()
    use gprec, only:dp=>q
    use gspci
    use gutil
    implicit none
#ifdef real_version
    real(dp)  &
#else
    complex(dp)  &
#endif
            &z_row(maxval(dmem%idx(dmem%nval_bot+1:dmem%nval_top+1) &
            &-dmem%idx(dmem%nval_bot:dmem%nval_top)))
    integer nstates,nnz,irun,ival,istate,i,jstate,j,p,q,r,s, &
            &itmp(4),isgn(4)

    nstates=dmem%idx(dmem%nval_top+1)-dmem%idx(dmem%nval_bot)
    dmem%ucsr%nrow=nstates
    dmem%ucsr%ncol=nstates
    allocate(dmem%ucsr%i(nstates+1)); dmem%ucsr%i(1)=1

    do irun=1,2
        if(irun==2)then
            allocate(dmem%ucsr%j(nnz),dmem%ucsr%a(nnz))
        endif
        nnz=0; istate=0
        do ival=dmem%nval_bot,dmem%nval_top
            do i=dmem%idx(ival),dmem%idx(ival+1)-1
                istate=istate+1
                z_row=0
                do p=1,dmem%norb
                    ! <v|p^\dagger
                    isgn(1)=1
                    itmp(1)=dmem%bs(i)
                    call act_state(itmp(1),p-1,.false.,isgn(1))
                    if(isgn(1)==0)cycle

                    ! one-body part
                    do q=1,dmem%norb
                        if(abs(dmem%h1e(p,q))<1.d-10)cycle
                        ! <v|p^\dagger q
                        isgn(2)=isgn(1)
                        itmp(2)=itmp(1)
                        call act_state(itmp(2),q-1,.true.,isgn(2))
                        if(isgn(2)==0)cycle
                        jstate=dmem%ibs(itmp(2))
                        if(jstate>i)cycle
                        jstate=jstate-dmem%idx(ival)+1
                        z_row(jstate)=z_row(jstate)+isgn(2)*dmem%h1e(p,q)
                    enddo

                    ! two-body
                    do q=1,p ! take care of factor 1/2
                        ! <v|p^\dagger q^\dagger
                        isgn(2)=isgn(1)
                        itmp(2)=itmp(1)
                        call act_state(itmp(2),q-1,.false.,isgn(2))
                        if(isgn(2)==0)cycle
                        do r=1,dmem%norb
                            ! <v|p^\dagger q^\dagger r 
                            isgn(3)=isgn(2)
                            itmp(3)=itmp(2)
                            call act_state(itmp(3),r-1,.true.,isgn(3))
                            if(isgn(3)==0)cycle
                            do s=1,dmem%norb
                                if(abs(dmem%v2e(p,s,q,r))<1.d-10)cycle
                                ! <v|p^\dagger q^\dagger r s
                                isgn(4)=isgn(3)
                                itmp(4)=itmp(3)
                                call act_state(itmp(4),s-1,.true.,isgn(4))
                                if(isgn(4)==0)cycle
                                jstate=dmem%ibs(itmp(4))
                                if(jstate>i)cycle
                                jstate=jstate-dmem%idx(ival)+1
                                z_row(jstate)=z_row(jstate)+isgn(4)* &
                                        &dmem%v2e(p,s,q,r)
                            enddo
                        enddo
                    enddo
                enddo
                if(irun==1)then
                    nnz=nnz+count(abs(z_row(1:i-dmem%idx(ival)+1))>1.d-10)
                    dmem%ucsr%i(istate+1)=nnz+1
                    cycle
                else
                    do j=dmem%idx(ival),i
                        jstate=j-dmem%idx(ival)+1
                        if(abs(z_row(jstate))<=1.d-10)cycle
                        nnz=nnz+1
                        dmem%ucsr%j(nnz)=j-dmem%idx(dmem%nval_bot)+1
                        dmem%ucsr%a(nnz)=z_row(jstate)
                    enddo
                endif
            enddo
        enddo
    enddo
    return

end subroutine calc_ucsr_spci


subroutine setup_hdns_spci_dlh(a,n)
    use gprec
    use gspci
    implicit none
    integer,intent(in)::n
#ifdef real_version
    real(q),intent(out)::a(n,n)
#else
    complex(q),intent(out)::a(n,n)
#endif

    integer i,j,inz,r,s,nfs,nfs_l,i1,j1,i1_,i2,ival

    a=0

    ! D_ij
    do i=1,dmem%norb
        do j=1,dmem%norb
            if(abs(dmem%daalpha(i,j))<1.d-10)cycle
            do inz=1,dmem%mcoo(j,i)%nnz
                r=dmem%mcoo(j,i)%i(inz)
                s=dmem%mcoo(j,i)%j(inz)
                a(r,s)=a(r,s)+dmem%mcoo(j,i)%a(inz)*dmem%daalpha(i,j)
                a(s,r)=a(s,r)+dmem%mcoo(j,i)%a(inz)*  &
#ifndef real_version
                        &conjg  &
#endif
                        &(dmem%daalpha(i,j))
            enddo
        enddo
    enddo

    ! lambda_c (la2)
    do i=1,dmem%norb
        do j=1,i
            if(abs(dmem%lambdac(i,j))<1.d-10)cycle
            do inz=1,dmem%nvcoo(j,i)%nnz
                r=dmem%nvcoo(j,i)%i(inz)
                s=dmem%nvcoo(j,i)%j(inz)
                a(r,s)=a(r,s)+dmem%nvcoo(j,i)%a(inz)*dmem%lambdac(i,j)
                if(i==j)cycle
                a(s,r)=a(s,r)+dmem%nvcoo(j,i)%a(inz)*  &
#ifndef real_version
                        &conjg  &
#endif
                        &(dmem%lambdac(i,j))
            enddo
        enddo
    enddo

    ! h_loc (ucsr)
    do ival=dmem%nval_bot,dmem%nval_top
        nfs=dmem%idx(ival+1)-dmem%idx(ival)
        nfs_l=dmem%idx_l(dmem%norb-ival+1)-dmem%idx_l(dmem%norb-ival)
        if(nfs_l<=0)cycle
        do i1=dmem%idx(ival),dmem%idx(ival+1)-1
            i1_=i1-dmem%idx(dmem%nval_bot)+1
            do inz=dmem%ucsr%i(i1_),dmem%ucsr%i(i1_+1)-1
                j1=dmem%ucsr%j(inz)+dmem%idx(dmem%nval_bot)-1
                do i2=dmem%idx_l(dmem%norb-ival), &
                        &dmem%idx_l(dmem%norb-ival+1)-1
                    r=dmem%i_phi(i1)%i(i2)
                    s=dmem%i_phi(j1)%i(i2)
                    if(r==0.or.s==0)cycle
                    a(r,s)=a(r,s)+dmem%ucsr%a(inz)
                    if(r==s)cycle
                    a(s,r)=a(s,r)+  &
#ifndef real_version
                            &conjg  &
#endif
                            &(dmem%ucsr%a(inz))
                enddo
            enddo
        enddo
    enddo
    return

end subroutine setup_hdns_spci_dlh


subroutine av1_gspci_dlh(v1,v2)
    use gprec
    use gspci
    use sparse
    implicit none
#ifdef real_version
    real(q),  &
#else
    complex(q),  &
#endif
            &intent(in)::v1(*)
#ifdef real_version
    real(q),  &
#else
    complex(q),  &
#endif
            &intent(inout)::v2(*)

    integer i,j

    ! D_ij
    do i=1,dmem%norb
        do j=1,dmem%norb
            if(abs(dmem%daalpha(i,j))<1.d-10)cycle
            call sp_amux (dmem%daalpha(i,j),dmem%mcoo(j,i),v1,v2)
            call coohmux(  &
#ifdef real_version
                    &dmem%daalpha(i,j), &
#else
                    &conjg(dmem%daalpha(i,j)), &
#endif
                    &dmem%mcoo(j,i),v1,v2)
        enddo
    enddo

    ! lambda_c (la2)
    do i=1,dmem%norb
        do j=1,i
            if(abs(dmem%lambdac(i,j))<1.d-10)cycle
            call sp_amux (dmem%lambdac(i,j),dmem%nvcoo(j,i),v1,v2)
            if(i==j)cycle
            call coohmux(dmem%lambdac(j,i),dmem%nvcoo(j,i),v1,v2)
        enddo
    enddo

    ! h_loc (ucsr)
    call act_ucsr_spci(v1,v2)
    return

end subroutine av1_gspci_dlh


subroutine act_ucsr_spci(v1,v2)
    use gprec
    use gspci
    implicit none
#ifdef real_version
    real(q),intent(in)::v1(*)
    real(q),intent(inout)::v2(*)
#else
    complex(q),intent(in)::v1(*)
    complex(q),intent(inout)::v2(*)
#endif

    integer ival,nfs_l,i1,i1_,i2,j1,inz,r,s

    do ival=dmem%nval_bot,dmem%nval_top
        nfs_l=dmem%idx_l(dmem%norb-ival+1)-dmem%idx_l(dmem%norb-ival)
        if(nfs_l<=0)cycle
        do i1=dmem%idx(ival),dmem%idx(ival+1)-1
            i1_=i1-dmem%idx(dmem%nval_bot)+1
            do inz=dmem%ucsr%i(i1_),dmem%ucsr%i(i1_+1)-1
                j1=dmem%ucsr%j(inz)+dmem%idx(dmem%nval_bot)-1
                do i2=dmem%idx_l(dmem%norb-ival), &
                        &dmem%idx_l(dmem%norb-ival+1)-1
                    r=dmem%i_phi(i1)%i(i2)
                    if(r<=0)then
                        cycle
                    endif
                    s=dmem%i_phi(j1)%i(i2)
                    if(s<=0)then
                        cycle
                    endif
                    v2(r)=v2(r)+dmem%ucsr%a(inz)*v1(s)
                    if(r==s)cycle
                    v2(s)=v2(s)+ &
#ifndef real_version
                        &conjg &
#endif
                        &(dmem%ucsr%a(inz))*v1(r)
                enddo
            enddo
        enddo
    enddo
    return

end subroutine act_ucsr_spci


subroutine calc_dm_spci()
    use gprec
    use gspci
    use sparse
    implicit none
    integer i,j,i_,j_

    dmem%dm=0

    ! c_i^\dagger c_j
    do i=1,dmem%norb
        do j=1,i
            call vh_sp_v(dmem%npcoo(i,j),dmem%v,dmem%dm(i,j))
            if(i==j)cycle
            dmem%dm(j,i)=  &
#ifndef real_version
                    &conjg  &
#endif
                    &(dmem%dm(i,j))
        enddo
    enddo

    ! f_j f_i^\dagger = \delta_i,j - f_i^\dagger f_j
    do i=1,dmem%norb
        i_=i+dmem%norb
        do j=1,i
            j_=j+dmem%norb
            call vh_sp_v(dmem%nvcoo(j,i),dmem%v,dmem%dm(i_,j_))
            if(i_==j_)then
                dmem%dm(i_,j_)=1-dmem%dm(i_,j_)
            else
                dmem%dm(i_,j_)=-dmem%dm(i_,j_)
                dmem%dm(j_,i_)=  &
#ifndef real_version
                        &conjg  &
#endif
                        &(dmem%dm(i_,j_))
            endif
        enddo
    enddo

    ! c_j^\dagger f_i
    do i=1,dmem%norb
        i_=i+dmem%norb
        do j=1,dmem%norb
            call vh_sp_v(dmem%mcoo(j,i),dmem%v,dmem%dm(j,i_))
            dmem%dm(i_,j)=  &
#ifndef real_version
                    &conjg  &
#endif
                    &(dmem%dm(j,i_))
        enddo
    enddo
    return

end subroutine calc_dm_spci


subroutine av_gspci(v1,v2,k,primme)
    use gprec
    use gspci
    use gconstant
    implicit none
    integer,intent(in)::k
    integer(8),intent(in)::primme
#ifdef real_version
    real(q),intent(in)::v1(*)
    real(q),intent(out)::v2(*)
#else
    complex(q),intent(in)::v1(*)
    complex(q),intent(out)::v2(*)
#endif

    integer n,i

    n=dmem%nstates
    do i=1,k
        v2(n*(i-1)+1:n*i)=0
        call av1_gspci(v1(n*(i-1)+1),v2(n*(i-1)+1))
    enddo
    return

end subroutine av_gspci


#ifndef real_version
subroutine chk_eval_j2_sumway(cm_vec,label)
    use gprec
    use gspci
    use sparse
    implicit none
    character(*),intent(in)::label
    complex(q),intent(in)::cm_vec(dmem%norb,dmem%norb,3)

    integer i
    complex(q) eval

    eval=0
    do i=1,3
        call calc_eval_gspci_j2op(cm_vec(:,:,i),dmem%v,eval)
    enddo
    write(0,*) label,"(sum_way) = ", eval
    return

end subroutine chk_eval_j2_sumway


subroutine av1_gspci_j2_sumway(lambda_j2,cm_vec,v1,v2)
    use gprec
    use gspci
    use sparse
    implicit none
    real(q),intent(in)::lambda_j2
    complex(q),intent(in)::cm_vec(dmem%norb,dmem%norb,3)
    complex(q),intent(in)::v1(*)
    complex(q),intent(inout)::v2(*)

    complex(q) cm(dmem%norb,dmem%norb)
    integer i

    if(abs(lambda_j2)>1.d-12)then
        do i=1,3
            cm=cm_vec(:,:,i)*sqrt(lambda_j2)
            call av1_gspci_j2op(cm,v1,v2)
        enddo
    endif
    return

end subroutine av1_gspci_j2_sumway


! Using equation J^2 = J-J+ + Jz^2 + Jz
subroutine av1_gspci_j2_npzway(lambda_j2,cm_n,cm_z,v1,v2)
    use gprec
    use gspci
    use sparse
    implicit none
    real(q),intent(in)::lambda_j2
    complex(q),intent(in)::cm_n(dmem%norb,dmem%norb),cm_z(dmem%norb,dmem%norb)
    complex(q),intent(in)::v1(*)
    complex(q),intent(inout)::v2(*)

    if(abs(lambda_j2)>1.d-12)then
        call av1_gspci_jz2_pjz(lambda_j2,cm_z,v1,v2)
        call av1_gspci_jnjp(lambda_j2,cm_n,v1,v2)
    endif
    return

end subroutine av1_gspci_j2_npzway


! A component of angular momentum vector acting on |v1>
! |v2> += \sum_{a,b}{cm_{a,b} (nc_{a,b} + nf_{a,b})} |v1>
subroutine av1_gspci_jop(cm,v1,v2)
    use gprec
    use gspci
    use sparse
    implicit none
    complex(q),intent(in)::cm(dmem%norb,dmem%norb)
    complex(q),intent(in)::v1(dmem%nstates)
    complex(q),intent(inout)::v2(dmem%nstates)

    integer i,j

    do i=1,dmem%norb; do j=1,i
        if(abs(cm(i,j))>1.d-12)then
            call sp_amux( cm(i,j),dmem%npcoo(i,j),v1,v2)
            call sp_amux(-cm(i,j),dmem%nvcoo(j,i),v1,v2)
            if(i==j)then
                v2=v2+v1*cm(i,j)
            endif
        endif
        if(i==j)cycle
        if(abs(cm(j,i))>1.d-12)then
            call coohmux( cm(j,i),dmem%npcoo(i,j),v1,v2)
            call coohmux(-cm(j,i),dmem%nvcoo(j,i),v1,v2)
        endif
    enddo; enddo
    return

end subroutine av1_gspci_jop


! Expectation value.
! <v2| \sum_{a,b}{cm_{a,b} (nc_{a,b} + nf_{a,b})} |v1>
subroutine zv2h_gspci_jop_v1(cm,v1,v2,zes)
    use gprec
    use gspci
    use sparse
    implicit none
    complex(q),intent(in)::cm(dmem%norb,dmem%norb)
    complex(q),intent(in)::v1(dmem%nstates),v2(dmem%nstates)
    complex(q),intent(inout)::zes
    complex(q),external::zdotc

    integer i,j

    do i=1,dmem%norb; do j=1,i
        if(abs(cm(i,j))>1.d-12)then
            call vh_sp_v( cm(i,j),dmem%npcoo(i,j),v1,v2,zes)
            call vh_sp_v(-cm(i,j),dmem%nvcoo(j,i),v1,v2,zes)
            if(i==j)then
                zes=zes+zdotc(dmem%nstates,v2,1,v1,1)*cm(i,j)
            endif
        endif
        if(i==j)cycle
        if(abs(cm(j,i))>1.d-12)then
            call vh_cooh_v( cm(j,i),dmem%npcoo(i,j),v1,v2,zes)
            call vh_cooh_v(-cm(j,i),dmem%nvcoo(j,i),v1,v2,zes)
        endif
    enddo; enddo
    return

end subroutine zv2h_gspci_jop_v1


! The square of one component of angular momentum operator scting on |v1>
! |v2> += (s/l/j)_{x/y/z}^2 |v1>
subroutine av1_gspci_j2op(cm,v1,v2)
    use gprec
    use gspci
    use sparse
    implicit none
    complex(q),intent(in)::cm(dmem%norb,dmem%norb)
    complex(q),intent(in)::v1(dmem%nstates)
    complex(q),intent(inout)::v2(dmem%nstates)

    complex(q) v1p(dmem%nstates)

    v1p=0
    call av1_gspci_jop(cm,v1,v1p)
    call av1_gspci_jop(cm,v1p,v2)
    return

end subroutine av1_gspci_j2op


! |v2> += lambda_j2 * (J_z^2 + Jz) |v1>
subroutine av1_gspci_jz2_pjz(lambda_j2,cm_z,v1,v2)
    use gprec
    use gspci
    use sparse
    implicit none
    real(q),intent(in)::lambda_j2
    complex(q),intent(in)::cm_z(dmem%norb,dmem%norb)
    complex(q),intent(in)::v1(dmem%nstates)
    complex(q),intent(inout)::v2(dmem%nstates)

    complex(q) v1p(dmem%nstates)
    complex(q) cm(dmem%norb,dmem%norb)

    v1p=0
    cm=cm_z*sqrt(lambda_j2)
    call av1_gspci_jop(cm,v1,v1p)
    v2=v2+sqrt(lambda_j2)*v1p
    call av1_gspci_jop(cm,v1p,v2)
    return

end subroutine av1_gspci_jz2_pjz


! |v2> += lambda_j2 * (J-J+) |v1>
subroutine av1_gspci_jnjp(lambda_j2,cm_n,v1,v2)
    use gprec
    use gspci
    use sparse
    implicit none
    real(q),intent(in)::lambda_j2
    complex(q),intent(in)::cm_n(dmem%norb,dmem%norb)
    complex(q),intent(in)::v1(dmem%nstates)
    complex(q),intent(inout)::v2(dmem%nstates)

    complex(q) v1p(dmem%nstates)
    complex(q) cm(dmem%norb,dmem%norb)

    v1p=0
    cm=conjg(transpose(cm_n))*sqrt(lambda_j2)
    call av1_gspci_jop(cm,v1,v1p)
    cm=cm_n*sqrt(lambda_j2)
    call av1_gspci_jop(cm,v1p,v2)
    return

end subroutine av1_gspci_jnjp


! The expectation value of square of one component of angular momentum 
! operator with respect to |v1>, i.e.,
! <v1| (s/l/j)_{x/y/z}^2 |v1>
subroutine calc_eval_gspci_j2op(cm,v1,zes)
    use gprec
    use gspci
    use sparse
    implicit none
    complex(q),intent(in)::cm(dmem%norb,dmem%norb)
    complex(q),intent(in)::v1(dmem%nstates)
    complex(q),intent(inout)::zes

    complex(q) v1p(dmem%nstates)

    v1p=0
    call av1_gspci_jop(cm,v1,v1p)
    call zv2h_gspci_jop_v1(cm,v1p,v1,zes)
    return

end subroutine calc_eval_gspci_j2op
#endif

subroutine calc_save_n_blks(imp)
    use gprec
    use sparse
    use gutil
    use gspci
    implicit none
    integer,intent(in)::imp

    integer ival,i,j
    ! have to be allocatable array!
    type(dcoo_matrix),allocatable::npcoo(:,:)

    allocate(npcoo(dmem%norb,dmem%norb))
    call gh5%fopen('HEmbedNij_'//trim(int_to_str(imp))//'.h5',1,"w", &
            &serialio=.true.)
    do ival=dmem%nval_bot,dmem%nval_top
        call gh5%gcreate('/valence_block_'//trim(int_to_str(ival)),1)
        call calc_npcoo_nblk(npcoo,ival)
        do i=1,dmem%norb; do j=1,i
            if(npcoo(i,j)%nnz==0)cycle
            call gh5%write_sp_matrix(npcoo(i,j),'/valence_block_'// &
                    &trim(int_to_str(ival))//'/n_'//trim(int_to_str(i-1))// &
                    &'_'//trim(int_to_str(j-1)),1,serialio=.true.)
            call dealloc_sp_matrix(npcoo(i,j))
        enddo; enddo
    enddo
    call gh5%fclose(1)
    deallocate(npcoo)
    return

end subroutine calc_save_n_blks


subroutine calc_save_rho_cp_blks()
    use gprec
    use sparse
    use gutil
    use ghdf5
    use gspci
    implicit none
    integer ival,nstates,i,j
#ifdef real_version
    real(q),  &
#else
    complex(q),  &
#endif
            &allocatable::rho(:,:)
    ! have to be allocatable array!
    type(dcoo_matrix),allocatable::npcoo(:,:)

    allocate(npcoo(dmem%norb,dmem%norb))
    call gh5%fopen('HEmbedAnalysis_'//trim(int_to_str(dmem%imp-1))//'.h5', &
            &1,"w",serialio=.true.)
    do ival=dmem%nval_bot,dmem%nval_top
        nstates=dmem%idx(ival+1)-dmem%idx(ival)
        allocate(rho(nstates,nstates))
        call calc_reduced_rho_nblk(rho,nstates,ival)
        if(maxval(abs(rho))>1.d-16)then
            call gh5%gcreate('/valence_block_'// &
                    &trim(int_to_str(ival)),1)
            call gh5%dwrite(rho,nstates,nstates,'/valence_block_'// &
                    &trim(int_to_str(ival))//'/RHO',1,serialio=.true.)
            call calc_npcoo_nblk(npcoo,ival)
            do i=1,dmem%norb; do j=1,i
                if(npcoo(i,j)%nnz==0)cycle
                call gh5%write_sp_matrix(npcoo(i,j),'/valence_block_'// &
                        &trim(int_to_str(ival))//'/NP_'//trim(int_to_str(i))// &
                        &'_'//trim(int_to_str(j)),1,serialio=.true.)
                call dealloc_sp_matrix(npcoo(i,j))
            enddo; enddo
        endif
        deallocate(rho)
    enddo
    call gh5%fclose(1)
    deallocate(npcoo)
    return

end subroutine calc_save_rho_cp_blks


subroutine calc_save_varrho_blks(imp)
    use gprec
    use sparse
    use gutil
    use ghdf5
    use gspci
    implicit none
    integer,intent(in)::imp

    integer ival,nstates,i,j
#ifdef real_version
    real(q),  &
#else
    complex(q),  &
#endif
            &allocatable::rho(:,:)
    ! have to be allocatable array!
    type(dcoo_matrix),allocatable::npcoo(:,:)

    allocate(npcoo(dmem%norb,dmem%norb))
    call gh5%fopen('HEmbedVarRho_'//trim(int_to_str(imp))//'.h5',1,"w", &
            &serialio=.true.)
    do ival=dmem%nval_bot,dmem%nval_top
        nstates=dmem%idx(ival+1)-dmem%idx(ival)
        allocate(rho(nstates,nstates))
        call calc_variational_reduced_rho_nblk(rho,nstates,ival)
        if(maxval(abs(rho))>1.d-16)then
            call gh5%gcreate('/valence_block_'// &
                    &trim(int_to_str(dmem%norb-ival)),1)
            call gh5%dwrite(rho,nstates,nstates,'/valence_block_'// &
                    &trim(int_to_str(dmem%norb-ival))//'/VRHO',1, &
                    &serialio=.true.)
        endif
        deallocate(rho)
    enddo
    call gh5%fclose(1)
    deallocate(npcoo)
    return

end subroutine calc_save_varrho_blks


subroutine calc_save_varrho_blks_slow(imp)
    use gprec
    use sparse
    use gutil
    use ghdf5
    use gspci
    implicit none
    integer,intent(in)::imp

    integer ival,nstates,i,j
#ifdef real_version
    real(q),  &
#else
    complex(q),  &
#endif
            &allocatable::rho(:,:)
    ! have to be allocatable array!
    type(dcoo_matrix),allocatable::npcoo(:,:)

    allocate(npcoo(dmem%norb,dmem%norb))
    call gh5%fopen('HEmbedVarRho_'//trim(int_to_str(imp))//'.h5',1,"w", &
            &serialio=.true.)
    call switch_vbasis_order(dmem%v,dmem%nstates)
    do ival=dmem%nval_bot,dmem%nval_top
        nstates=dmem%idx(ival+1)-dmem%idx(ival)
        allocate(rho(nstates,nstates))
        call calc_reduced_rho_nblk(rho,nstates,ival)
        ! rotate rho from particle basis to hole basis.
        call transform_p2h(rho,nstates,ival)
        if(maxval(abs(rho))>1.d-16)then
            call gh5%gcreate('/valence_block_'// &
                    &trim(int_to_str(dmem%norb-ival)),1)
            call gh5%dwrite(rho,nstates,nstates,'/valence_block_'// &
                    &trim(int_to_str(dmem%norb-ival))//'/VRHO',1, &
                    &serialio=.true.)
        endif
        deallocate(rho)
    enddo
    call gh5%fclose(1)
    deallocate(npcoo)
    return

end subroutine calc_save_varrho_blks_slow


subroutine calc_npcoo_nblk(npcoo,ival)
    use gprec
    use gspci
    use sparse
    use gutil
    implicit none
    integer,intent(in)::ival
    type(dcoo_matrix),intent(out)::npcoo(dmem%norb,dmem%norb)

    integer nstates,ia1,ia2,i1,isum,itmp1,isgn,ibs,nbase

    nstates=dmem%idx(ival+1)-dmem%idx(ival)
    nbase=dmem%idx(ival)-1
    do ia1=1,dmem%norb; do ia2=1,ia1
        call alloc_sp_matrix(npcoo(ia1,ia2),nstates,nstates,nstates)
        isum=0
        ! < state | c_ia1^\dagger c_ia2
        do i1=dmem%idx(ival),dmem%idx(ival+1)-1
            itmp1=dmem%bs(i1)
            isgn=1
            call act_state(itmp1,ia1-1,.false.,isgn)
            if(isgn==0)cycle
            call act_state(itmp1,ia2-1,.true.,isgn)
            if(isgn==0)cycle
            ibs=dmem%ibs(itmp1)
            isum=isum+1
            npcoo(ia1,ia2)%i(isum)=i1-nbase
            npcoo(ia1,ia2)%j(isum)=ibs-nbase
            npcoo(ia1,ia2)%a(isum)=real(isgn,q)
        enddo
        if(isum==0)then
            call dealloc_sp_matrix(npcoo(ia1,ia2))
        else
            npcoo(ia1,ia2)%nnz=isum
        endif
    enddo; enddo
    return

end subroutine calc_npcoo_nblk


! Generate reduced many-body denity matrix in valence block n.
subroutine calc_reduced_rho_nblk(rho,n,ival)
    use gprec
    use gspci
    implicit none
    integer,intent(in)::n,ival
#ifdef real_version
    real(q),  &
#else
    complex(q),  &
#endif
            &intent(out)::rho(n,n)
    
    integer i,i_,j,j_,ndim,nfs_l
    real(q) sz
#ifdef real_version
    real(q),external::ddot
#else
    complex(q),external::zdotc
#endif

    rho=0
    nfs_l=dmem%idx_l(dmem%norb-ival+1)-dmem%idx_l(dmem%norb-ival)
    if(nfs_l<=0)return

    do i=dmem%idx(ival),dmem%idx(ival+1)-1
        sz=dmem%bs_sz(i)
        ndim=dmem%i_phi(i)%imax-dmem%i_phi(i)%imin+1
        if(ndim<=0)cycle
        i_=i-dmem%idx(ival)+1
        do j=dmem%idx(ival),i
            if(abs(sz-dmem%bs_sz(j))>1.d-6)then
                cycle
            endif
            j_=j-dmem%idx(ival)+1
            rho(i_,j_)=  &
#ifdef real_version
                    &ddot&
#else
                    &zdotc&
#endif
                    &(ndim,dmem%v(dmem%i_phi(j)%imin),1,&
                    &dmem%v(dmem%i_phi(i)%imin),1)
            if(i==j)cycle
            rho(j_,i_)=  &
#ifndef real_version
                    &conjg  &
#endif
                    &(rho(i_,j_))
        enddo
    enddo
    return


end subroutine calc_reduced_rho_nblk


! rotate solution vector phy(x)var to var(x)phy.
subroutine switch_vbasis_order(v,n)
    use gprec
    use gspci
    implicit none
    integer,intent(in)::n
#ifdef real_version
    real(q),  &
#else
    complex(q),  &
#endif
            &intent(inout)::v(n)

#ifdef real_version
    real(q) vp(n)
#else
    complex(q) vp(n)
#endif
    integer ival,i,j,ij,ji,nstate1,nstate2

    vp=0
    do ival=dmem%nval_bot,dmem%nval_top
        do i=dmem%idx(ival),dmem%idx(ival+1)-1
            do j=dmem%idx_l(dmem%norb-ival),dmem%idx_l(dmem%norb-ival+1)-1
                nstate1=dmem%i_phi(i)%i(j)
                ij=dmem%ibs_l(ibits(not(dmem%bs(i)),0,dmem%norb))
                ji=dmem%ibs(ibits(not(dmem%bs_l(j)),0,dmem%norb))
                nstate2=dmem%i_phi(ji)%i(ij)
                ! sanity check 
                if((nstate1>0).neqv.(nstate2>0))then
                    write(0,"('basis not consistent!')")
                    stop 3
                endif
                if(nstate1>0)then
                    vp(nstate2)=v(nstate1)
                endif
            enddo
        enddo
    enddo
    v=vp
    return
    

end subroutine switch_vbasis_order


subroutine transform_p2h(rho,n,ival)
    use gprec
    use gspci
    use gutil
    implicit none
    integer,intent(in)::n,ival
#ifdef real_version
    real(q),  &
#else
    complex(q),  &
#endif
            &intent(inout)::rho(n,n)

    integer i,j,i_,j_
#ifdef real_version
    real(q) u(n,n)
#else
    complex(q) u(n,n)
#endif

    if(n<=0)return
    u=0._q
    do i=dmem%idx(ival),dmem%idx(ival+1)-1
        i_=i-dmem%idx(ival)+1
        j=dmem%ibs_l(ibits(not(dmem%bs(i)),0,dmem%norb))
        j_=j-dmem%idx(dmem%norb-ival)+1
        u(i_,j_)=1._q
    enddo
    call uhau(rho,u,n,n)
    return

end subroutine transform_p2h


! Generate reduced variational many-body denity matrix in valence block n.
subroutine calc_variational_reduced_rho_nblk(rho,n,ival)
    use gprec
    use gspci
    implicit none
    integer,intent(in)::n,ival
#ifdef real_version
    real(q),  &
#else
    complex(q),  &
#endif
            &intent(out)::rho(n,n)
    
    integer i,ip,ip_,j,jp,jp_,ndim,nfs_l,istate,jstate
    real(q) sz

    rho=0
    nfs_l=dmem%idx_l(dmem%norb-ival+1)-dmem%idx_l(dmem%norb-ival)
    if(nfs_l<=0)return

    do i=dmem%idx(ival),dmem%idx(ival+1)-1
        sz=dmem%bs_sz(i)
        do ip=dmem%idx_l(dmem%norb-ival),dmem%idx_l(dmem%norb-ival+1)-1
            istate=dmem%i_phi(i)%i(ip)
            if(istate==0)then
                cycle
            endif
            ip_=ip-dmem%idx_l(dmem%norb-ival)+1
            do jp=dmem%idx_l(dmem%norb-ival),dmem%idx_l(dmem%norb-ival+1)-1
                jstate=dmem%i_phi(i)%i(jp)
                if(jstate==0)then
                    cycle
                endif
                jp_=jp-dmem%idx_l(dmem%norb-ival)+1
                rho(ip_,jp_)=rho(ip_,jp_)+dmem%v(istate)*  &
#ifndef real_version
                        &conjg  &
#endif
                        &(dmem%v(jstate))
            enddo
        enddo
    enddo
    return


end subroutine calc_variational_reduced_rho_nblk

