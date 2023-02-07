module localstder
    use gprec, only: q
    use localstore
    use ghdf5
    use gconstant, only: zi
    use gutil, only: uhau,disimix,dpisimix,pfa_pa,dsimix,dpsimix,  &
            &calc_2ndderivative_sroot_entropyn_at_h,int_to_str
    implicit none
    private

    type,extends(localstore_ob),public::localstder_ob
#ifdef real_version
        real(q),  &
#else
        complex(q),  &
#endif
                &pointer :: pdmpr(:,:)=>null(), &
                &pdmpl(:,:)=>null(), &  ! density matrix derivatives
                &pd0pl(:,:)=>null(), pd0pr(:,:)=>null(), &
                &pdpl(:,:)=>null(), pdpr(:,:)=>null(), &
                &plcpl(:,:)=>null(), plcpr(:,:)=>null(), &
                &pd_coefpl(:,:)=>null(), pd_coefpr(:,:)=>null(), &
                &pd0_coefpl(:,:)=>null(), pd0_coefpr(:,:)=>null(), &
                &pncvpl(:,:)=>null(), pncvpr(:,:)=>null(), &
                &pr0pl(:,:)=>null(), pr0pr(:,:)=>null(), &
                &prpl(:,:)=>null(), prpr(:,:)=>null(), &
                &pr_coefpl(:,:)=>null(), pr_coefpr(:,:)=>null()
        real(q),pointer :: plc_coefpl(:,:)=>null(),plc_coefpr(:,:)=>null(), &
                &pncv_coefpl(:,:)=>null(), pncv_coefpr(:,:)=>null(), &
                &pdm_coefpl(:,:)=>null(), pdm_coefpr(:,:)=>null()
        integer::nxl,nxr,nx

        contains

        procedure::init_der=>init_localstder
        procedure::dfah
        procedure::big_to_bkmat
        procedure::calc_pdm_pp
        procedure::calc_pdplr
        procedure::calc_plcplr
        procedure::calc_pdlc_pp
        procedure::h5save_pdlc
        procedure::setup2_pnrplr
    end type

    contains


    subroutine set_nxlr(this)
    class(localstder_ob)::this

    this%nxl=this%hm_l%dimhst
    this%nxr=this%hm_r%dimhst*this%r_factor
    this%nx=this%hm%dimhst
    return

    end subroutine set_nxlr

        
    subroutine init_localstder(this)
    class(localstder_ob)::this

    call set_nxlr(this)
    allocate(this%pdmpl(this%na2112,this%nxl), &
            &this%pdmpr(this%na2112,this%nxr), &
            &this%pd0pl(this%na2112,this%nxl), &
            &this%pd0pr(this%na2112,this%nxr), &
            &this%pdpl(this%na2112,this%nxl), &
            &this%pdpr(this%na2112,this%nxr),  &
            &this%plcpl(this%na2112,this%nxl),  &
            &this%plcpr(this%na2112,this%nxr))
    this%pdmpl=0; this%pdmpr=0; this%pd0pl=0; this%pd0pr=0
    this%pdpl=0; this%pdpr=0; this%plcpl=0; this%plcpr=0
    return

    end subroutine init_localstder


    subroutine setup2_pnrplr(this,gh5,io)
    class(localstder_ob)::this
    class(hdf5_ob) :: gh5
    integer,intent(in) :: io

    call set_nxlr(this)
    call h5load_pdlc(this,gh5)
    call h5load_pdr_impwise(this,gh5)
    call calc_pnr0plr(this)
    call calc_prplr(this)
    call calc_pnr_pp(this,io)
    return

    end subroutine setup2_pnrplr


    subroutine dfah(this,v,loemat,h,pfph,n,ltrans,wt)
    ! calculate p f(A+t*h) / p t | t=0, given the eigen-vectors v, 
    ! and the loewner matrix.
    class(localstder_ob)::this
    integer,intent(in)::n
    logical,intent(in)::ltrans
    complex(q),intent(in)::v(n,n)
    real(q),intent(in)::loemat(n,n),wt
    complex(q),intent(inout)::h(n,n)
#ifdef real_version
    real(q),  &
#else
    complex(q),  &
#endif
            &intent(inout)::pfph(this%na2112)

    ! rotate to eigen-reprsentation of h_qp
    call uhau(h,v,n,n)
    h=h*loemat
    ! rotate back
    call uhau(h,v,n,n,trul='n',trur='c')
    call big_to_bkmat(this,h,pfph,n,ltrans,wt)
    return

    end subroutine dfah


    subroutine big_to_bkmat(this,a,a_bk,n,ltrans,wt)
    ! set block matrices a_bk from big matrix a.
    class(localstder_ob)::this
    integer,intent(in)::n
    logical,intent(in)::ltrans
    real(q),intent(in)::wt
    complex(q),intent(in)::a(n,n)
#ifdef real_version
    real(q),  &
#else
    complex(q),  &
#endif
            &target,intent(inout)::a_bk(this%na2112)

    integer i,ibase,ijbase,na2,na22
#ifdef real_version
    real(q),  &
#else
    complex(q),  &
#endif
            &pointer::p_a(:,:)

    ibase=0
    ijbase=0
    do i=1,this%num_imp
        na2=this%na2_list(i)
        na22=na2*na2
        p_a(1:na2,1:na2)=>a_bk(ijbase+1:ijbase+na22)
        if(ltrans)then
            p_a=p_a+transpose(a(ibase+1:ibase+na2,ibase+1:ibase+na2))*wt
        else
            p_a=p_a+a(ibase+1:ibase+na2,ibase+1:ibase+na2)*wt
        endif
        ibase=ibase+na2
        ijbase=ijbase+na22
    enddo
    return

    end subroutine big_to_bkmat


    subroutine calc_pdm_pp(this,io)
    class(localstder_ob)::this
    integer,intent(in)::io

    allocate(this%pd0_coefpl(this%hm_r%dimhst,this%nxl), &
            &this%pd0_coefpr(this%hm_r%dimhst,this%nxr),&
            &this%pdm_coefpl(this%hm_l%dimhst,this%nxl),  &
            &this%pdm_coefpr(this%hm_l%dimhst,this%nxr))
    this%pd0_coefpl=0; this%pd0_coefpr=0
    this%pdm_coefpl=0; this%pdm_coefpr=0
    call symm_herm_matrices_list_c(this,this%pdmpl,this%hm_l, &
            &this%nxl,io,'pdmpl',.true.,this%pdm_coefpl)
    call symm_herm_matrices_list_c(this,this%pdmpr,this%hm_l, &
            &this%nxr,io,'pdmpr',.true.,this%pdm_coefpr)
    call symm_matrices_list(this,this%pd0pl,this%hm_r, &
            &this%nxl,io,'pd0pl',.false.,.false.,this%pd0_coefpl)
    call symm_matrices_list(this,this%pd0pr,this%hm_r, &
            &this%nxr,io,'pd0pr',.false.,.false.,this%pd0_coefpr)
    return

    end subroutine calc_pdm_pp


    subroutine calc_pdlc_pp(this,io)
    class(localstder_ob)::this
    integer,intent(in)::io

    allocate(this%pd_coefpl(this%hm_r%dimhst,this%nxl), &
            &this%pd_coefpr(this%hm_r%dimhst,this%nxr),&
            &this%plc_coefpl(this%hm%dimhst,this%nxl),  &
            &this%plc_coefpr(this%hm%dimhst,this%nxr))
    this%pd_coefpl=0; this%pd_coefpr=0
    this%plc_coefpl=0; this%plc_coefpr=0
    call symm_matrices_list(this,this%pdpl,this%hm_r, &
            &this%nxl,io,'pdpl',.false.,.false.,this%pd_coefpl)
    call symm_matrices_list(this,this%pdpr,this%hm_r, &
            &this%nxr,io,'pdpr',.false.,.false.,this%pd_coefpr)
    call symm_herm_matrices_list_c(this,this%plcpl,this%hm, &
            &this%nxl,io,'plcpl',.false.,this%plc_coefpl)
    call symm_herm_matrices_list_c(this,this%plcpr,this%hm, &
            &this%nxr,io,'plcpr',.false.,this%plc_coefpr)
    return

    end subroutine calc_pdlc_pp


    subroutine calc_pnr_pp(this,io)
    class(localstder_ob)::this
    integer,intent(in)::io

    allocate(this%pdm_coefpl(this%hm_l%dimhst,this%nxl), &
            &this%pdm_coefpr(this%hm_l%dimhst,this%nxr), &
            &this%pncv_coefpl(this%hm_l%dimhst,this%nxl), &
            &this%pncv_coefpr(this%hm_l%dimhst,this%nxr), &
            &this%pr_coefpl(this%hm_r%dimhst,this%nxl), &
            &this%pr_coefpr(this%hm_r%dimhst,this%nxr))
    call symm_matrices_list(this,this%prpl,this%hm_r, &
            &this%nxl,io,'prpl',.false.,.false.,this%pr_coefpl)
    call symm_matrices_list(this,this%prpr,this%hm_r, &
            &this%nxr,io,'prpr',.false.,.false.,this%pr_coefpr)
    call symm_herm_matrices_list_c(this,this%pncvpl,this%hm_l, &
            &this%nxl,io,'pncvpl',.true.,this%pncv_coefpl)
    call symm_herm_matrices_list_c(this,this%pncvpr,this%hm_l, &
            &this%nxr,io,'pncvpr',.true.,this%pncv_coefpr)
    call symm_herm_matrices_list_c(this,this%pdmpl,this%hm_l, &
            &this%nxl,io,'pdmpl',.true.,this%pdm_coefpl)
    call symm_herm_matrices_list_c(this,this%pdmpr,this%hm_l, &
            &this%nxr,io,'pdmpr',.true.,this%pdm_coefpr)
    return

    end subroutine calc_pnr_pp


    subroutine symm_matrices_list(this,a_list,mb,n,io,dname,ltrans,lherm,c)
    class(localstder_ob)::this
    integer,intent(in)::io,n
    type(matrix_basis),intent(in)::mb
    character(*),intent(in)::dname
#ifdef real_version
    real(q),  &
#else
    complex(q),  &
#endif
            &intent(inout)::a_list(this%na2112,n)
#ifdef real_version
    real(q),  &
#else
    complex(q),  &
#endif
            &optional,intent(out)::c(mb%dimhst,n)
    logical,intent(in)::ltrans,lherm

    integer i
    real(q)::maxerr=1.e-8_q
#ifdef real_version
    real(q)  &
#else
    complex(q)  &
#endif
            &copy(this%na2112)

    do i=1,n
        copy=a_list(:,i)
        call this%symm_dm_across_atoms(a_list(:,i))
        if(present(c))then
            call this%hm_expand_all(a_list(:,i),c(:,i),mb,0,ltrans,lherm)
        else
            call this%hm_expand_all(a_list(:,i),mb,ltrans,lherm)
        endif
        maxerr=max(maxerr,maxval(abs(a_list(:,i)-copy)))
    enddo
    call report_maxerr(maxerr,dname,io)
    this%symerr=max(this%symerr,maxerr)
    return

    end subroutine symm_matrices_list
   

    subroutine symm_herm_matrices_list_c(this,a_list,mb,n,io,dname,ltrans,c)
    class(localstder_ob)::this
    integer,intent(in)::io,n
    type(matrix_basis),intent(in)::mb
    character(*),intent(in)::dname
#ifdef real_version
    real(q),  &
#else
    complex(q),  &
#endif
            &intent(inout)::a_list(this%na2112,n)
    real(q),target,intent(out)::c(mb%dimhst,n)
    logical,intent(in)::ltrans

#ifdef real_version
    real(q),  &
#else
    complex(q),  &
#endif
            &pointer::c_buf(:,:)

#ifdef real_version
    c_buf=>c
#else
    allocate(c_buf(mb%dimhst,n))
#endif
    call symm_matrices_list(this,a_list,mb,n,io,dname,ltrans,.true.,c_buf)
#ifndef real_version
    c=c_buf
#endif
    return

    end subroutine symm_herm_matrices_list_c


    subroutine calc_pdplr(this)
    class(localstder_ob)::this

    call calc_pd(this,this%pd0pl,this%pdmpl,this%pdpl,this%nxl,0)
    call calc_pd(this,this%pd0pr,this%pdmpr,this%pdpr,this%nxr,0)
    return

    end subroutine calc_pdplr


    subroutine calc_prplr(this)
    class(localstder_ob)::this

    allocate(this%prpl(this%na2112,this%nxl), &
            &this%prpr(this%na2112,this%nxr))
    this%prpl=0; this%prpr=0
    call calc_pd(this,this%pr0pl,this%pdmpl,this%prpl,this%nxl,1)
    call calc_pd(this,this%pr0pr,this%pdmpr,this%prpr,this%nxr,1)
    return

    end subroutine calc_prplr


    subroutine calc_pd(this,pd0,pdm,pd,n,mode)
    ! get pdpl/r from pd0(m)pl/r for mode=0, otherwise,
    ! prpl/r from pr0(m)pl/r
    class(localstder_ob)::this
    integer,intent(in)::n,mode
#ifdef real_version
    real(q),  &
#else
    complex(q),  &
#endif
            &intent(in)::pd0(this%na2112,n),pdm(this%na2112,n)
#ifdef real_version
    real(q),  &
#else
    complex(q),  &
#endif
            &target,intent(out)::pd(this%na2112,n)

    integer i,ibase,na2,na22,j
#ifdef real_version
    real(q),  &
#else
    complex(q),  &
#endif
            &target::zbuf(this%na2max**2),zbuf1(this%na2max**2)
#ifdef real_version
    real(q),  &
#else
    complex(q),  &
#endif
            &pointer::pn(:,:),ph(:,:),p_d(:,:)

    pd=0
    ibase=0
    do i=1,this%num_imp
        na2=this%na2_list(i)
        na22=na2*na2
        pn(1:na2,1:na2)=>zbuf(1:na22)
        ph(1:na2,1:na2)=>zbuf1(1:na22)
        do j=1,n
            p_d(1:na2,1:na2)=>pd(ibase+1:ibase+na22,j)
            ! get pd0 component
            zbuf(1:na22)=pd0(ibase+1:ibase+na22,j)
            ! derivative part 1
            if(mode==0)then
                p_d=matmul(this%co(i)%isimix,pn)
            else
                p_d=matmul(transpose(this%co(i)%isimix),pn)
            endif
            ! get pdm component
            zbuf1(1:na22)=pdm(ibase+1:ibase+na22,j)
            call pfa_pa(this%co(i)%nks,pn,ph,na2,disimix,dpisimix)
            ! p isimix / p l/r
            ! derivative part 2
            if(mode==0)then
                p_d=p_d+matmul(pn,this%co(i)%d0)
            else
                p_d=p_d+matmul(transpose(pn),this%co(i)%r0)
            endif
        enddo
        ibase=ibase+na22
    enddo
    return

    end subroutine calc_pd


    subroutine calc_plcplr(this)
    ! part 1 of plc / pl(r)
    class(localstder_ob)::this

    integer i,ibase,ihbase,ieq,na2,na22,j,k,k1,k2,i2,j2
    real(q) res
#ifdef real_version
    real(q),  &
#else
    complex(q),  &
#endif
            &target::zbuf(this%na2max**2),zbuf1(this%na2max**2), &
            &zbuf2(this%na2max**2),zbuf3(this%na2max**2)
#ifdef real_version
    real(q),  &
#else
    complex(q),  &
#endif
            &pointer::pn(:,:),ph(:,:),p_lc(:,:),pr(:,:),pd(:,:),prd(:,:)

    ! inequivalent atom index
    ieq=0
    ibase=0
    ihbase=0
    do i=1,this%num_imp
        if(this%imap_list(i)==i)then
            ieq=ieq+1
        endif
        na2=this%na2_list(i)
        na22=na2*na2
        ! from lambda part first
        pn(1:na2,1:na2)=>zbuf(1:na22)
        ph(1:na2,1:na2)=>zbuf1(1:na22)
        pr(1:na2,1:na2)=>zbuf2(1:na22)
        prd(1:na2,1:na2)=>zbuf3(1:na22)
        do j=1,this%co(i)%dim_hs_l
            p_lc(1:na2,1:na2)=>this%plcpl(ibase+1:ibase+na22,ihbase+j)
            p_lc=p_lc-this%co(i)%hs_l(:,:,j)
            ! nks is expanded in terms of h.t
            ph=transpose(this%co(i)%hs_l(:,:,j))
            ! part 1 from r.transpose(d)
            ! ph stays the same.
            call pfa_pa(this%co(i)%nks,pn,ph,na2,dsimix,dpsimix) ! p f / p d_n
            ! p lc / p l
            do k=1,this%nxl
                ! r . transpose(pd/pl), no pr / pl
                pd(1:na2,1:na2)=>this%pdpl(ibase+1:ibase+na22,k)
                prd=matmul(this%co(i)%r,transpose(pd))
                res=-real(sum(pn*prd),q)*2
                p_lc(1:na2,1:na2)=>this%plcpl(ibase+1:ibase+na22,k)
                p_lc=p_lc+this%co(i)%hs_l(:,:,j)*res
            enddo
            ! p lc / p r
            ! unique impurity index
            i2=1
            j2=0
            do k=1,this%nxr
                ! r . transpose(pd/pr)
                pd(1:na2,1:na2)=>this%pdpr(ibase+1:ibase+na22,k)
                prd=matmul(this%co(i)%r,transpose(pd))
                if(j2==this%co(i2)%dim_hs_r*this%r_factor)then
                    i2=i2+1
                    j2=0
                endif
                j2=j2+1
                if(i2==ieq)then
                    ! pr . d
                    if(this%r_factor==2)then  ! complex-r version
                        k1=(j2+1)/2
                        k2=mod(j2,2)
                    else
                        k1=j2
                        k2=1
                    endif
                    pr=this%co(i)%hs_r(:,:,k1)
                    if(k2==0)then
                        pr=pr*zi
                    endif
                    prd=prd+matmul(pr,transpose(this%co(i)%d))
                endif
                res=-real(sum(pn*prd),q)*2
                p_lc(1:na2,1:na2)=>this%plcpr(ibase+1:ibase+na22,k)
                p_lc=p_lc+this%co(i)%hs_l(:,:,j)*res
            enddo
            ! part 2 from sqrt(a(1-n)), pn gets overwritten
            prd=matmul(this%co(i)%r,transpose(this%co(i)%d))
            do k=1,this%nxl
                pd(1:na2,1:na2)=>this%pdmpl(ibase+1:ibase+na22,k)
                call calc_2ndderivative_sroot_entropyn_at_h(this%co(i)%nks,  &
                        &ph,pd,pn,na2)
                res=-real(sum(pn*prd),q)*2
                p_lc(1:na2,1:na2)=>this%plcpl(ibase+1:ibase+na22,k)
                p_lc=p_lc+this%co(i)%hs_l(:,:,j)*res
            enddo
            do k=1,this%nxr
                ! r . transpose(pd/pr)
                pd(1:na2,1:na2)=>this%pdmpr(ibase+1:ibase+na22,k)
                call calc_2ndderivative_sroot_entropyn_at_h(this%co(i)%nks,  &
                        &ph,pd,pn,na2)
                res=-real(sum(pn*prd),q)*2
                p_lc(1:na2,1:na2)=>this%plcpr(ibase+1:ibase+na22,k)
                p_lc=p_lc+this%co(i)%hs_l(:,:,j)*res
            enddo
        enddo
        ibase=ibase+na22
        if(i<this%num_imp)then
            if(this%imap_list(i)<this%imap_list(i+1))then
                ihbase=ihbase+this%co(i)%dim_hs_l
            endif
        endif
    enddo
    return

    end subroutine calc_plcplr


    subroutine h5save_pdlc(this,gh5)
    class(localstder_ob)::this
    class(hdf5_ob) :: gh5

    call gh5%fopen("GDLDeri.h5",1,"w")
    ! write the derivative of dm, d_coed and lambdac_coef wrt r and lambda.
    call gh5%dwrite(this%pdmpl,this%na2112,this%nxl,'/PDMPL',1)
    call gh5%dwrite(this%pdmpr,this%na2112,this%nxr,'/PDMPR',1)
    call gh5%dwrite(this%pd_coefpl,this%hm_r%dimhst,this%nxl,'/PD_COEFPL',1)
    call gh5%dwrite(this%pd_coefpr,this%hm_r%dimhst,this%nxr,'/PD_COEFPR',1)
    call gh5%dwrite(this%plc_coefpl,this%hm%dimhst,this%nxl,'/PLC_COEFPL',1)
    call gh5%dwrite(this%plc_coefpr,this%hm%dimhst,this%nxr,'/PLC_COEFPR',1)
    call gh5%dwrite(this%pd0_coefpl,this%hm_r%dimhst,this%nxl,'/PD0_COEFPL',1)
    call gh5%dwrite(this%pd0_coefpr,this%hm_r%dimhst,this%nxr,'/PD0_COEFPR',1)
    call gh5%dwrite(this%pdm_coefpl,this%hm_l%dimhst,this%nxl,'/PDM_COEFPL',1)
    call gh5%dwrite(this%pdm_coefpr,this%hm_l%dimhst,this%nxr,'/PDM_COEFPR',1)
    call gh5%fclose(1)
    return

    end subroutine h5save_pdlc


    subroutine h5load_pdlc(this,gh5)
    class(localstder_ob)::this
    class(hdf5_ob) :: gh5

    allocate(this%pdmpl(this%na2112,this%nxl),  &
            &this%pdmpr(this%na2112,this%nxr), &
            &this%pd_coefpl(this%hm_r%dimhst,this%nxl), &
            &this%pd_coefpr(this%hm_r%dimhst,this%nxr), &
            &this%plc_coefpl(this%hm%dimhst,this%nxl), &
            &this%plc_coefpr(this%hm%dimhst,this%nxr))
    call gh5%fopen("GDLDeri.h5",1,"r")
    ! write the derivative of dm, d_coed and lambdac_coef wrt r and lambda.
    call gh5%read(this%pdmpl,this%na2112,this%nxl,'/PDMPL',1)
    call gh5%read(this%pdmpr,this%na2112,this%nxr,'/PDMPR',1)
    call gh5%read(this%pd_coefpl,this%hm_r%dimhst,this%nxl,'/PD_COEFPL',1)
    call gh5%read(this%pd_coefpr,this%hm_r%dimhst,this%nxr,'/PD_COEFPR',1)
    call gh5%read(this%plc_coefpl,this%hm%dimhst,this%nxl,'/PLC_COEFPL',1)
    call gh5%read(this%plc_coefpr,this%hm%dimhst,this%nxr,'/PLC_COEFPR',1)
    call gh5%fclose(1)
    return

    end subroutine h5load_pdlc


    subroutine h5load_pdr_impwise(this,gh5)
    class(localstder_ob)::this
    class(hdf5_ob) :: gh5

    integer i,na2,nx,nxr,imap

    call gh5%fopen("HEmbed.h5",1,"r")
    do i=1,this%num_imp
        na2=this%na2_list(i)
        nx=this%hm%dim_hs(i)
        nxr=this%hm_r%dim_hs(i)*this%r_factor
        imap=this%imap_list(i)
        if(imap==i)then
            allocate(this%co(i)%pncvpd(na2,na2,this%nxr), &
                    &this%co(i)%pncvplc(na2,na2,this%nx), &
                    &this%co(i)%pr0pd(na2,na2,this%nxr), &
                    &this%co(i)%pr0plc(na2,na2,this%nx))
            call gh5%read(this%co(i)%pncvpd,na2,na2,this%nxr,'/impurity_'// &
                    &trim(int_to_str(i-1))//"/ans/NCV_DERI_D",1)
            call gh5%read(this%co(i)%pncvplc,na2,na2,this%nx,'/impurity_'// &
                    &trim(int_to_str(i-1))//"/ans/NCV_DERI_LC",1)
            call gh5%read(this%co(i)%pr0pd,na2,na2,this%nxr,'/impurity_'// &
                    &trim(int_to_str(i-1))//"/ans/R0_DERI_D",1)
            call gh5%read(this%co(i)%pr0plc,na2,na2,this%nx,'/impurity_'// &
                    &trim(int_to_str(i-1))//"/ans/R0_DERI_LC",1)
        else
            this%co(i)%pncvpd=>this%co(imap)%pncvpd
            this%co(i)%pncvplc=>this%co(imap)%pncvplc
            this%co(i)%pr0pd=>this%co(imap)%pr0pd
            this%co(i)%pr0plc=>this%co(imap)%pr0plc
        endif
    enddo
    call gh5%fclose(1)
    return

    end subroutine h5load_pdr_impwise


    subroutine calc_pnr0plr(this)
    class(localstder_ob)::this

    integer i,na2,ibase,inbase,irbase
    
    allocate(this%pr0pl(this%na2112,this%nxl), &
            &this%pr0pr(this%na2112,this%nxr), &
            &this%pncvpl(this%na2112,this%nxl), &
            &this%pncvpr(this%na2112,this%nxr))
    this%pr0pl=0; this%pr0pr=0; this%pncvpl=0; this%pncvpr=0
    ibase=0
    inbase=0
    irbase=0
    do i=1,this%num_imp
        na2=this%na2_list(i)
        call calc_pnr0plr1(this,this%pr0pl,this%co(i)%pr0pd, &
                &this%co(i)%pr0plc, &
                &this%pd_coefpl(irbase+1:irbase+this%hm_r%dim_hs(i),:), &
                &this%plc_coefpl(inbase+1:inbase+this%hm%dim_hs(i),:), &
                &i,na2,this%nxl,ibase,this%hm_r%dim_hs(i),this%hm%dim_hs(i))
        call calc_pnr0plr1(this,this%pncvpl,this%co(i)%pncvpd, &
                &this%co(i)%pncvplc, &
                &this%pd_coefpl(irbase+1:irbase+this%hm_r%dim_hs(i),:), &
                &this%plc_coefpl(inbase+1:inbase+this%hm%dim_hs(i),:), &
                &i,na2,this%nxl,ibase,this%hm_r%dim_hs(i),this%hm%dim_hs(i))
        call calc_pnr0plr1(this,this%pr0pr,this%co(i)%pr0pd, &
                &this%co(i)%pr0plc, &
                &this%pd_coefpr(irbase+1:irbase+this%hm_r%dim_hs(i),:), &
                &this%plc_coefpr(inbase+1:inbase+this%hm%dim_hs(i),:), &
                &i,na2,this%nxr,ibase,this%hm_r%dim_hs(i),this%hm%dim_hs(i))
        call calc_pnr0plr1(this,this%pncvpr,this%co(i)%pncvpd, &
                &this%co(i)%pncvplc, &
                &this%pd_coefpr(irbase+1:irbase+this%hm_r%dim_hs(i),:), &
                &this%plc_coefpr(inbase+1:inbase+this%hm%dim_hs(i),:), &
                &i,na2,this%nxr,ibase,this%hm_r%dim_hs(i),this%hm%dim_hs(i))
        ibase=ibase+na2*na2
        if(i<this%num_imp)then
            if(this%imap_list(i+1)==i+1)then
                irbase=irbase+this%hm_r%dim_hs(i)
                inbase=inbase+this%hm%dim_hs(i)
            endif
        endif
    enddo
    return

    end subroutine calc_pnr0plr


    subroutine calc_pnr0plr1(this,dm_list,pdmpd,pdmplc,pd_coef,plc_coef, &
            &i,na2,nx,ibase,njd,njl)
    class(localstder_ob)::this
    integer,intent(in)::na2,i,nx,ibase,njd,njl
#ifdef real_version
    real(q),  &
#else
    complex(q),  &
#endif
            &target,intent(inout)::dm_list(this%na2112,nx)
#ifdef real_version
    real(q),  &
#else
    complex(q),  &
#endif
            &intent(in)::pdmpd(na2,na2,njd), &
            &pdmplc(na2,na2,njl), &
            &pd_coef(njd,nx)
    real(q),intent(in)::plc_coef(njl,nx)

    integer j,ix,na22
#ifdef real_version
    real(q),  &
#else
    complex(q),  &
#endif
            &pointer::p_buf(:,:)

    na22=na2*na2
    do ix=1,nx
        p_buf(1:na2,1:na2)=>dm_list(ibase+1:ibase+na22,ix)
        do j=1,njd
            p_buf=p_buf+pdmpd(:,:,(j-1)*this%r_factor+1)* &
                    &real(pd_coef(j,ix),q)
#ifndef real_version
                p_buf=p_buf+pdmpd(:,:,(j-1)*this%r_factor+2)* &
                        &aimag(pd_coef(j,ix))
#endif
        enddo
        do j=1,njl
            p_buf=p_buf+pdmplc(:,:,j)*plc_coef(j,ix)
        enddo
    enddo
    return

    end subroutine calc_pnr0plr1


end module localstder
