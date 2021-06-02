from pygrisb.run.timing import timeit

import warnings
warnings.filterwarnings("ignore")


@timeit
def driver(
        imp=0,
        iembeddiag=-1,
        path="./",
        mpiexec=["mpirun", "-np", "1"],
        nval=5,
        nval_bot=0,
        nval_top=None,
        edinit=0,
        cygtail="",
        jac=None,
        analysis=False,
        ):

    if iembeddiag == -10:
        from pygrisb.gsolver.gs_ml import gs_h5ml_fsoc_krr
        gs_h5ml_fsoc_krr(imp=imp, jac=jac)

    elif iembeddiag == -11:
        from pygrisb.gsolver.gs_ml import gs_h5ml_fsoc_nm1d
        gs_h5ml_fsoc_nm1d(nval=nval, imp=imp, jac=jac)

    elif iembeddiag == -13:
        from pygrisb.gsolver.gs_ml import gs_h5ml_fsoc_nm2d
        gs_h5ml_fsoc_nm2d(nval=nval, imp=imp, jac=jac)

    elif iembeddiag == -14:
        from pygrisb.gsolver.gs_ml import gs_h5ml_fsoc_nm3d
        gs_h5ml_fsoc_nm3d(nval=nval, imp=imp, jac=jac)

    elif iembeddiag in [-1, -2, -3, -4, -201]:
        import subprocess
        siembeddiag = str(-iembeddiag)
        if analysis:
            siembeddiag = "-" + siembeddiag
        cmd = [f'{path}/exe_spci{cygtail.lower()}', '-i', str(imp+1), '-m',
                siembeddiag, '-e', str(edinit)]
        subprocess.run(cmd, check=True)

    elif iembeddiag == 1:
        from pygrisb.gsolver.gsolver import gsolver_h5trans_ed
        gs = gsolver_h5trans_ed(imp=imp, mpiexec=mpiexec, path=path,
                nval_bot=nval_bot, nval_top=nval_top)
        gs.driver()

    elif iembeddiag == 101:
        from pygrisb.gsolver.gsolver import gsolver_h5ed
        gs = gsolver_h5ed(imp=imp, mpiexec=mpiexec, path=path,
                nval_bot=nval_bot, nval_top=nval_top)
        gs.driver()

    elif iembeddiag == 50:
        from pygrisb.gsolver.vhs.vhs import gsolver_vhs
        gs = gsolver_vhs(
                imp=imp,
                read_v2e=True,
                )
        gs.driver()

    else:
        raise NotImplementedError( \
                f"iembeddiag = {iembeddiag} not available!")
