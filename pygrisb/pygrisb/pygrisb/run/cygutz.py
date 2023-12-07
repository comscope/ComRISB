import sys, subprocess, argparse, os
import numpy as np
from scipy import optimize
import pygrisb.run.info as info
from pygrisb.run.config import get_options
from pygrisb.io.h5io import h5open


def get_groot():
    path = os.environ.get('WIEN_GUTZ_ROOT',
            os.environ.get('COMRISB_BIN', "./"))
    return path


def get_cygtail():
    cygtail = ""
    if os.path.isfile("GParam.h5"):
        with h5open("GParam.h5", "r") as f:
            val = f["/impurity_0/V2E"][0,0,0,0]
            if val.dtype == np.float:
                cygtail = "R"
    return cygtail

# todo
# solve the list of inequivalent embedding Hamiltonian.
#   iembeddiag = choice embedding h solver (int). 10 mean Hartree-Fock.
#       if iembeddiag=10=Hartree-Fock then do nothing
def solve_hembed_list(
        mpi=["mpirun", "-np", "1"],
        path=get_groot(),
        iembeddiag=None,
        edinit=0,
        cygtail=get_cygtail(),
        jac=None,
        nval_list=None,
        save_hemb_history=False,
        **kwargs,
        ):
    
    from pygrisb.gsolver.gs_driver import driver as gs_driver
    
    # _iembeddiag = GParam.h5/giembeddiag
    with h5open("GParam.h5", "r") as f:
        imap_list = f["/"].attrs["imap_list"]
        nval_bot_list = f["/"].attrs["nval_bot_list"]
        nval_top_list = f["/"].attrs["nval_top_list"]
        _iembeddiag = f["/"].attrs["giembeddiag"]
        
    if iembeddiag is None:
        iembeddiag = _iembeddiag # _iembeddiag = GParam.h5/giembeddiag
        
    if iembeddiag == 10: # if iembeddiag=10=Hartree-Fock then do nothing
        # hartree-fock
        return
    
    if nval_list is None:
        # default for ml solver.
        nval_list = [5 for i in imap_list] # imap_list = GParam.h5/imap_list
    elif len(nval_list) == 1:
        nval = nval_list[0]
        nval_list = [nval for i in imap_list]
    else:
        assert(len(nval_list) == len(imap_list)), \
                f"nval_list = {nval_list}: inconsistent length!"
                
    print(f" nval_list for ml solver: {nval_list}")
    
    for i, imap in enumerate(imap_list): # imap_list = GParam.h5/imap_list
        if i==imap:
            
            # driver calls one of the following methods:
            # iembeddiag == -10: gs_ml.gs_h5ml_fsoc_krr(imp=imp, jac=jac) -> kernel ridge machine learning
            # iembeddiag == -11: gs_ml.gs_h5ml_fsoc_nm1d(nval=nval, imp=imp, jac=jac) -> machine learning with normal-mode expansion at 1st order for multivariate interpolation
            # iembeddiag == -13: gs_ml.gs_h5ml_fsoc_nm2d(nval=nval, imp=imp, jac=jac) -> machine learning with normal-mode expansion at 2nd order for multivariate interpolation
            # iembeddiag == -14: gs_ml.gs_h5ml_fsoc_nm3d(nval=nval, imp=imp, jac=jac) -> machine learning with normal-mode expansion at 3rd order for multivariate interpolation
            # iembeddiag in [-1, -2, -3, -4, -201]: runs exe_spci{cygtail.lower()} with some arguments
            # iembeddiag == 1: gsolver_h5trans_ed(imp=imp, mpiexec=mpiexec, path=path, nval_bot=nval_bot, nval_top=nval_top) 
            #       -> runs exe_ed, parallel exact diagonalization solver with hamiltonian transformed
            #           into complex spherical hamrmonics basis with spin-faster index.
            # iembeddiag == 101: gsolver_h5ed(imp=imp, mpiexec=mpiexec, path=path, nval_bot=nval_bot, nval_top=nval_top)
            #       -> runs exe_ed, parallel exact diagonalization solver
            # iembeddiag == 50: gsolver_vhs(imp=imp,]read_v2e=True,) -> I can't find the source code.
            # iembeddiag == 10: Hartree-Fock, so do nothing.
            gs_driver(imp=i,
                    iembeddiag=iembeddiag,
                    path=path,
                    mpiexec=mpi,
                    nval=nval_list[i],
                    nval_bot=nval_bot_list[i],# nval_bot_list = GParam.h5/nval_bot_list

                    nval_top=nval_top_list[i],# nval_top_list = GParam.h5/nval_top_list
                    edinit=edinit,
                    cygtail=cygtail,
                    jac=jac,
                    )

    with h5open("HEmbed.h5", "a") as f:
        for i, imap in enumerate(imap_list): # imap_list = GParam.h5/imap_list
            if i == imap:
                if "ans" in f[f"/impurity_{i}"]:
                    del f[f"/impurity_{i}/ans"]
                with h5open(f"HEmbedRes_{i}.h5", "r") as fp:
                    f.copy(fp["/"], f"/impurity_{i}/ans")
        if save_hemb_history:
            with h5open("HEmbed_history.h5", "a") as fh:
                for i, imap in enumerate(imap_list):
                    if i == imap:
                        # to c-convention
                        d_array = [f[f"/impurity_{i}/D"][()].T]
                        if f"/impurity_{i}" in fh:
                            lcombine = True
                            d_history = fh[f"/impurity_{i}/d_history"][()]
                            d_array = np.concatenate((d_history, d_array),
                                    axis=0)
                            del fh[f"/impurity_{i}/d_history"]
                        else:
                            lcombine = False
                        fh[f"/impurity_{i}/d_history"] = d_array
                        # to c-convention
                        lambda_arrary = [f[f"/impurity_{i}/LAMBDA"][()].T]
                        if lcombine:
                            lambda_history = fh[
                                    f"/impurity_{i}/lambda_history"][()]
                            lambda_arrary = np.concatenate(
                                    (lambda_history, lambda_arrary), axis=0)
                            del fh[f"/impurity_{i}/lambda_history"]
                        fh[f"/impurity_{i}/lambda_history"] = lambda_arrary
                        # to c-convention
                        dm_array = [f[f"/impurity_{i}/ans/DM"][()].T]
                        if lcombine:
                            dm_history = fh[
                                    f"/impurity_{i}/dm_history"][()]
                            dm_array = np.concatenate(
                                    (dm_history, dm_array), axis=0)
                            del fh[f"/impurity_{i}/dm_history"]
                        fh[f"/impurity_{i}/dm_history"] = dm_array
                        # to c-convention
                        h1e_arrary = [f[f"/impurity_{i}/H1E"][()].T]
                        if lcombine:
                            h1e_history = fh[
                                    f"/impurity_{i}/h1e_history"][()]
                            h1e_arrary = np.concatenate(
                                    (h1e_history, h1e_arrary), axis=0)
                            del fh[f"/impurity_{i}/h1e_history"]
                        fh[f"/impurity_{i}/h1e_history"] = h1e_arrary

                        emol_array = [f[f"/impurity_{i}/ans"].attrs["emol"]]
                        if lcombine:
                            emol_history = fh[
                                    f"/impurity_{i}/emol_history"][()]
                            emol_array = np.concatenate(
                                    (emol_history, emol_array))
                            del fh[f"/impurity_{i}/emol_history"]
                        fh[f"/impurity_{i}/emol_history"] = emol_array

    # exec_cygband updates quasiparticle bands:
#   - Sets up GIter.h5 : (a) writes x to GIter.h5/x_val, 
#                        (b) increments GIter.h5/iter
#                        (c) removes GIter.h5/v_err
#   - Runs CyGBand{cygtail} -w 0, and if "jac" = False then it adds -d 1                        
def exec_cygband(
        x,
        mpi=["mpirun", "-np", "1"],
        path=".",
        cygtail="",
        **kwargs,
        ):
    
    # set the new trial solution vector
    with h5open("GIter.h5", "a") as f:
        if x is not None:
            f["/x_val"][()] = x
        # iteration counter
        f["/"].attrs["iter"] += 1
        if "/v_err" in f:
            del f["/v_err"]
            
    # calculate the updated qp bands and qp local quantities.
    cmd = mpi + [f"{path}/CyGBand{cygtail}", "-w", "0"]
    
    if kwargs.get("jac", False):
        cmd += ["-d", "1"]
        
    print(" ".join(cmd))
    
    subprocess.run(cmd, check=True)


def exec_cygerr(
        mpi=["mpirun", "-np", "1"],
        path=".",
        cygtail="",
        **kwargs,
        ):
    # calculate error vector
    if len(mpi) > 0:
        cmd = mpi[:-1] + ["1"]
    else:
        cmd = []
    cmd +=  [f"{path}/CyGErr{cygtail}", "-w", "0"]
    if kwargs.get("jac", False):
        cmd += ["-d", "1"]
    # get the error vector
    print(" ".join(cmd))
    subprocess.run(cmd, check=True)
    with h5open("GIter.h5", "r") as f:
        verr = f["/v_err"][()]
        if kwargs.get("jac", False):
            jac = f["/pfpx"][()]
        else:
            jac = None
    maxerr = np.max(np.abs(verr))
    print(f" iter = {info.iter} maxerr = {maxerr:.2e}")
    sys.stdout.flush()
    info.iter += 1
    if jac is None:
        return verr
    else:
        return verr, jac


def gfun_solver(x,
        args={"mpi":["mpirun", "-np", "1"],
                "path":".",
                "edinit":0,
                "cygtail":""},
        ):
    '''Gutzwiller vector error function. Using sparse matrix full-CI.
    '''
    # exec_cygband updates quasiparticle bands:
    #   - Sets up GIter.h5 : (a) writes x to GIter.h5/x_val, 
    #                        (b) increments GIter.h5/iter
    #                        (c) removes GIter.h5/v_err
    #   - Runs CyGBand{cygtail} -w 0, and if "jac" = False then it adds -d 1                        
    exec_cygband(x, **args)
    
    # solve the list of inequivalent embedding Hamiltonian.
    solve_hembed_list(**args)
    
    # calculate error vector
    verr = exec_cygerr(**args)
    # adding absolute convergence criteria
    maxerr = np.max(np.abs(verr))
    # if total energy needed each step.
    if(maxerr < 100):
        cmd = args["mpi"] + [f"{args['path']}/CyGFinal{args['cygtail']}",
                "-r", "0"]
        print(" ".join(cmd))
        subprocess.run(cmd, check=True)

    # dot not mess up jacobian
    if args.get("adjust_verr", True):
        if maxerr < 6.e-6:
            verr *= 0
        elif maxerr < 5.e-5 and info.iter > 300:
            verr *= 0
        elif maxerr < 1.e-4 and info.iter > 400:
            verr *= 0
        elif maxerr < 1.e-3 and info.iter > 600:
            verr *= 0
    return verr


def ndlc_fun(x,
        args={"mpi":["mpirun", "-np", "1"],
                "path":"./",
                "cygtail":""},
        ):
    
    # exec_cygband updates quasiparticle bands:
    #   - Sets up GIter.h5 : (a) writes x to GIter.h5/x_val, 
    #                        (b) increments GIter.h5/iter
    #                        (c) removes GIter.h5/v_err
    #   - Runs CyGBand{cygtail} -w 0, and if "jac" = False then it adds -d 1                            
    exec_cygband(x, **args)
    
    dsetname1 = args.get("dsetname1", "/d_coef")
    dsetname2 = args.get("dsetname2", "/lc_coef")
    with h5open("HEmbed.h5", "r") as f:
        dset1 = f[dsetname1][()]
        dset2 = f[dsetname2][()]
    verr = np.concatenate((dset1, dset2))
    return verr

# todo
def get_jacobian_analytical(
        x,
        args={"mpi":["mpirun", "-np", "1"],
                "path":"./",
                "edinit":0,
                "cygtail":""},
        ):
    
    # exec_cygband updates quasiparticle bands:
    #   - Sets up GIter.h5 : (a) writes x to GIter.h5/x_val, 
    #                        (b) increments GIter.h5/iter
    #                        (c) removes GIter.h5/v_err
    #   - Runs CyGBand{cygtail} -w 0                           
    exec_cygband(x, **args, jac=True)
    
    # solve the list of inequivalent embedding Hamiltonian.
    solve_hembed_list(**args, jac=True)
    
    # calculate error vector
    _, jac = exec_cygerr(**args, jac=True)
    
    print(" analytical jacobian called!")
    return jac


options_dict = {"hybr": {'eps':1.e-6, 'xtol':1.e-5, "factor": 1},
        "excitingmixing": {
                "maxiter": 10,
                "fatol": 0.00001,  # control the energy convegence.
                "xtol": 0.001,
                "jac_options":{"alpha": 0.2,
                        # "alphamax": 5.0,
                        "alphamax": 5.0,
                        }},
        }


def simplemixing(
         fun,
         x0,
         args={"mpi":["mpirun", "-np", "1"], "path":"."},
         tol=0.001,
         maxiter=50,
         alpha=0,  # mixing factor
         ):
    x = np.copy(x0)
    for i in range(maxiter):
        verr = fun(x, args=args)
        x += verr*alpha
    return verr



def gsolve_nleqns(
        fun,
        method="hybr",
        args={"mpi":["mpirun", "-np", "1"], "path":"."},
        jac=None,
        tol=1.e-6,
        ):
    '''main driver to solve the set of Gutzwiller nonlinear equations.
    '''
    info.iter = 0
    # get initil solution vector
    with h5open("GIter.h5", "r") as f:
        if "x_val" in f:
            x0 = f["/x_val"][()]
        else:
            x0 = None
    verr = fun(x0, args=args)
    gerr = np.max(np.abs(verr))
    if gerr > tol:
        conf = get_options()

        if method == "simplemixing":
            verr = simplemixing(
                    fun,
                    x0,
                    args=args,
                    tol=tol,
                    maxiter=conf.get("maxiter", 50),
                    alpha=conf.get("alpha", 0.2),
                    )
        else:
            if gerr < 1.e-4:
                # close enough, seem better
                method = "hybr"
                args["edinit"] = 0

            options = conf.get(method, options_dict.get(method, {}))
            if len(options) > 0:
                print("min-options:", options)

            res = optimize.root(
                    fun,
                    x0,
                    args=args,
                    method=method,
                    jac=jac,
                    options=options,
                    tol=tol,
                    )
            verr = res.fun
    return verr


# get_jacobian_numerical returns jac^T, where
#   jac = (vec_fun(x_list+delta)-vec_fun(x_list))/delta
#   The default values are delta = 1.e-4 and vec_fun=gfun_solver
def get_jacobian_numerical(
        x_list,
        args={"mpi":["mpirun", "-np", "1"],
        "path":"./",},
        delta=1.e-4,
        vec_fun=gfun_solver,
        ):
    
    jac = []
    
    args["adjust_verr"] = False
    
    vref = vec_fun(x_list, args=args)
    
    for i,_ in enumerate(x_list):
        xs = x_list.copy()
        xs[i] += delta
        v = vec_fun(xs, args=args)
        jac.append((v-vref)/delta)
        
    # recover the intial point
    with h5open("GIter.h5", "a") as f:
       f["/x_val"][()] = x_list
       
    print(" numerical jacobian called!")
    
    return np.asarray(jac).T

# get_inline_args (1) reads the following things off the command line:
#       (2) If jac = -1 then it sets args.jac = get_jacobian_numerical
#           If jac =  1 then it sets args.jac = get_jacobian_analytical
#           Otherwise        it sets args.jac = None
#   path = executible path
#   iembeddiag = choice embedding h solver (int)
#   edinit = ed solver init guess (int)
#   mpi = mpirun prefix (str)
#   rmethod = nonlinear solver method (str)
#   jac = jacobian method: 0->no, 1->analytical, -1: numerical
#   tol = tolerance (float)
#   iupdaterho = update rho (int): 0, default: no; 1: kswt; 2: no kswt
def get_inline_args():
    g_root = get_groot()
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--path", type=str, default=g_root,
            help="executible path (str)")
    parser.add_argument("-d", "--iembeddiag", type=int, default=None,
            help="choice embedding h solver (int)")
    parser.add_argument("-e", "--edinit", type=int, default=0,
            help="ed solver init guess (int)")
    parser.add_argument("-m", "--mpi", type=str, default="mpirun -np 1",
            help="mpirun prefix (str)")
    parser.add_argument("-r", "--rmethod", type=str, default="hybr",
            help="nonlinear solver method (str)")
    parser.add_argument("-j","--jac", type=int, default=0,
            help="jacobian method: 0->no, 1->analytical, -1: numerical.")
    parser.add_argument("-t","--tol", type=float, default=1e-6,
            help="tolerance (float).")
    parser.add_argument("-u", "--iupdaterho", type=int, default=0,
            help="update rho (int): 0, default: no; 1: kswt; 2: no kswt")
    args = parser.parse_args()
    if args.jac == -1:
        args.jac = get_jacobian_numerical
    elif args.jac == 1:
        args.jac = get_jacobian_analytical
    else:
        args.jac = None
    args.mpi = args.mpi.split()
    return args


def run_solve_hembed_list():
    
    # get_inline_args (1) reads the following things off the command line:
    #       (2) If jac = -1 then it sets args.jac = get_jacobian_numerical
    #           If jac =  1 then it sets args.jac = get_jacobian_analytical
    #           Otherwise        it sets args.jac = None
    #   path = executible path
    #   iembeddiag = choice embedding h solver (int)
    #   edinit = ed solver init guess (int)
    #   mpi = mpirun prefix (str)
    #   rmethod = nonlinear solver method (str)
    #   jac = jacobian method: 0->no, 1->analytical, -1: numerical
    #   tol = tolerance (float)
    #   iupdaterho = update rho (int): 0, default: no; 1: kswt; 2: no kswt
    args = get_inline_args()

    # solve the list of inequivalent embedding Hamiltonian.    
    solve_hembed_list(
            path=args.path,
            iembeddiag=args.iembeddiag,
            edinit=args.edinit,
            mpi=args.mpi,
            cygtail=get_cygtail(),
            )
    
def run_cygutz(
        path=get_groot(),
        rmethod="hybr",
        jac=0,
        mpi=["mpirun", "-np", "1"], # overridden by comdmft
        iembeddiag=None,
        edinit=0,
        tol=1e-6,
        cmdlargs=True,
        iupdaterho=0, # iupdaterho=2 when called from comdmft
        ):
    
    if cmdlargs:
        
        # get_inline_args (1) reads the following things off the command line:
        #       (2) If jac = -1 then it sets args.jac = get_jacobian_numerical
        #           If jac =  1 then it sets args.jac = get_jacobian_analytical
        #           Otherwise        it sets args.jac = None
        #   path = executible path
        #   iembeddiag = choice embedding h solver (int)
        #   edinit = ed solver init guess (int)
        #   mpi = mpirun prefix (str)
        #   rmethod = nonlinear solver method (str)
        #   jac = jacobian method: 0->no, 1->analytical, -1: numerical
        #   tol = tolerance (float)
        #   iupdaterho = update rho (int): 0, default: no; 1: kswt; 2: no kswt
        args = get_inline_args()
        path = args.path
        rmethod = args.rmethod
        jac = args.jac
        mpi = args.mpi
        iembeddiag = args.iembeddiag
        edinit = args.edinit
        tol = args.tol
        iupdaterho=args.iupdaterho

    driver_cygutz(
            path=path,
            rmethod=rmethod,
            jac=jac,
            mpi=mpi,
            iembeddiag=iembeddiag,
            edinit=edinit,
            cygtail=get_cygtail(),
            tol=tol,
            iupdaterho=iupdaterho,
            )


def driver_cygutz(
        path=".",
        rmethod="hybr",
        jac=None,
        mpi=["mpirun", "-np", "1"],
        iembeddiag=None,
        edinit=0,
        cygtail=get_cygtail(),
        tol=1e-6,
        iupdaterho=0,
        ):
    # set initial guess
    cmd = mpi + [f"{path}/CyGInit{cygtail}"]
    print(" ".join(cmd))
    print(f" nonlinear solver: {rmethod}")
    subprocess.run(cmd, check=True)
    # at least run info.min_iter times
    for i in range(info.min_iter):
        verr = gsolve_nleqns(gfun_solver,
                method=rmethod,
                jac=jac,
                args={
                        "mpi":mpi,
                        "path":path,
                        "iembeddiag":iembeddiag,
                        "edinit":edinit,
                        "cygtail":cygtail,
                        },
                tol=tol,
                )
        if info.iter > info.min_iter:
            break
    maxerr = np.max(np.abs(verr))
    print(f" gerr = {maxerr}.")
    cmd = mpi + [f"{path}/CyGFinal{cygtail}", "-r", str(iupdaterho)]
    print(" ".join(cmd))
    subprocess.run(cmd, check=True)
    return maxerr


def check_jacobian():
    
    # get_inline_args (1) reads the following things off the command line:
    #       (2) If jac = -1 then it sets args.jac = get_jacobian_numerical
    #           If jac =  1 then it sets args.jac = get_jacobian_analytical
    #           Otherwise        it sets args.jac = None
    #   path = executible path
    #   iembeddiag = choice embedding h solver (int)
    #   edinit = ed solver init guess (int)
    #   mpi = mpirun prefix (str)
    #   rmethod = nonlinear solver method (str)
    #   jac = jacobian method: 0->no, 1->analytical, -1: numerical
    #   tol = tolerance (float)
    #   iupdaterho = update rho (int): 0, default: no; 1: kswt; 2: no kswt    
    args = get_inline_args()
    
    # get the point around which the jacobian is calculated.
    with h5open("GIter.h5", "r") as f:
        x0 = f["/x_val"][()]
    args_ = {
            "mpi":args.mpi,
            "path":args.path,
            "iembeddiag":args.iembeddiag,
            "edinit":args.edinit,
            "cygtail":get_cygtail(),
            }
    
    # get_jacobian_numerical returns jac^T, where
    #   jac = (vec_fun(x_list+delta)-vec_fun(x_list))/delta
    #   The default values are delta = 1.e-4 and vec_fun=gfun_solver
    jac_num = get_jacobian_numerical(
            x0.copy(),
            args=args_,
            delta=1.e-5,
            )
    print(" Jacobian-numeric:")
    for row in jac_num:
        print("".join(f"{a:8.4f} " for a in row))

    jac_ana = get_jacobian_analytical(
            x0.copy(),
            args=args_,
            )
    print(" Jacobian-analytical:")
    for row in jac_ana:
        print("".join(f"{a:8.4f} " for a in row))

    diff = jac_num - jac_ana
    print(" Jacobian difference:")
    for row in diff:
        print("".join(f"{a:8.4f} " for a in row))
    print(f" max diff = {np.max(np.abs(diff)):8.4f}")


def check_jac_ndlc(dsetnames1=["/d_coef", "/lc_coef"],
        dsetnames2=["/PD_COEFPR", "/PD_COEFPL",
                "/PLC_COEFPR", "/PLC_COEFPL"]):
            
    # get_inline_args (1) reads the following things off the command line:
    #       (2) If jac = -1 then it sets args.jac = get_jacobian_numerical
    #           If jac =  1 then it sets args.jac = get_jacobian_analytical
    #           Otherwise        it sets args.jac = None
    #   path = executible path
    #   iembeddiag = choice embedding h solver (int)
    #   edinit = ed solver init guess (int)
    #   mpi = mpirun prefix (str)
    #   rmethod = nonlinear solver method (str)
    #   jac = jacobian method: 0->no, 1->analytical, -1: numerical
    #   tol = tolerance (float)
    #   iupdaterho = update rho (int): 0, default: no; 1: kswt; 2: no kswt            
    args = get_inline_args()
    
    # get the point around which the jacobian is calculated.
    with h5open("GIter.h5", "r") as f:
        x0 = f["/x_val"][()]
    args_ = {
            "mpi":args.mpi,
            "path":args.path,
            "iembeddiag":args.iembeddiag,
            "cygtail":get_cygtail(),
            "dsetname1":dsetnames1[0],
            "dsetname2":dsetnames1[1],
            }
    
    # get_jacobian_numerical returns jac^T, where
    #   jac = (vec_fun(x_list+delta)-vec_fun(x_list))/delta
    #   The default values are delta = 1.e-4 and vec_fun=gfun_solver
    jac_num = get_jacobian_numerical(
            x0.copy(),
            args=args_,
            vec_fun=ndlc_fun,
            delta=1.e-4,
            )
    print(" Jacobian-numeric:")
    print(f" max: {np.max(jac_num):8.2e}, min {np.min(jac_num):8.2e}")
    for row in jac_num:
        print("".join(f"{a:8.4f} " for a in row))

    # analytical jacobian
    # exec_cygband updates quasiparticle bands:
    #   - Sets up GIter.h5 : (a) writes x to GIter.h5/x_val, 
    #                        (b) increments GIter.h5/iter
    #                        (c) removes GIter.h5/v_err
    #   - Runs CyGBand{cygtail} -w 0, and if "jac" = False then it adds -d 1                            
    exec_cygband(x0, **args_, jac=True)
    
    with h5open("GDLDeri.h5", "r") as f:
        dset1 = f[dsetnames2[0]][()].T
        dset2 = f[dsetnames2[1]][()].T
        dset3 = f[dsetnames2[2]][()].T
        dset4 = f[dsetnames2[3]][()].T
    pd = np.concatenate((dset1, dset2), axis=1)
    plc = np.concatenate((dset3, dset4), axis=1)
    jac_ana = np.concatenate((pd, plc))
    print(" Jacobian-analytical:")
    print(f" max: {np.max(jac_ana):8.2e}, min {np.min(jac_ana):8.2e}")
    for row in jac_ana:
        print("".join(f"{a:8.4f} " for a in row))
    diff = jac_num - jac_ana
    print(" Jacobian difference:")
    print(f" max: {np.max(diff):8.2e}, min {np.min(diff):8.2e}")
    for row in diff:
        print("".join(f"{a:8.4f} " for a in row))


def check_jac_dlc():
    check_jac_ndlc()


def check_jac_dm():
    check_jac_ndlc(
            dsetnames1=["/d0_coef", "/nks_coef"],
            dsetnames2=["/PD0_COEFPR", "/PD0_COEFPL",
                    "/PDM_COEFPR", "/PDM_COEFPL"])
