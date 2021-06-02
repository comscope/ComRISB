import os, sys, glob, h5py, socket, shutil, time, re, \
        subprocess, warnings, argparse
import numpy as np
from collections import deque
from pygrisb.io.fio import file_exists
import pygrisb.run.environ as env
from pygrisb.run.cygutz import gfun_solver, gsolve_nleqns, get_cygtail
from pygrisb.dft.util import get_dir_list


def get_file_info(fname, unit, idmf, case, scratch, so, para, cmplx, _band,
        updn, dnup):
    '''help function to setup informations in def file.
    '''
    if 'in2' == fname:
        return [unit, "'{}.in2{}'".format(case, cmplx), "'old'",
                "'formatted'", 0]
    elif 'inso' == fname:
        return [unit, "'{}.inso'".format(case), "'unknown'", "'formatted'", 0]
    elif 'indmfl' == fname:
        return [unit, "'{}.indmfl'".format(case), "'old'", "'formatted'", 0]
    elif 'outputdmfupdn' == fname:
        return [unit, "'{}.outputdmf{}{}'".format(case, idmf, updn), \
                "'unknown'", "'formatted'", 0]
    elif 'in1c' == fname:
        return [unit, "'{}.in1c'".format(case), "'unknown'", "'formatted'", 0]
    elif 'vectorupdn' == fname:
        return [unit, "'{}/{}.vector{}{}{}'".format(scratch, case, so, \
                updn, para), "'unknown'","'unformatted'",9000]
    elif 'vectordnup' == fname:
        return [unit, "'{}/{}.vector{}{}{}'".format(scratch, case, so, \
                dnup, para), "'unknown'","'unformatted'",9000]
    elif 'klist' == fname:
        return [unit, "'{}.klist{}'".format(case, _band), "'old'", \
                "'formatted'", 0]
    elif 'kgen' == fname:
        return [unit, "'{}.kgen'".format(case), "'unknown'", "'formatted'", 0]
    elif 'vspupdn' == fname:
        return [unit, "'{}.vsp{}'".format(case, updn), "'old'", \
                "'formatted'", 0]
    elif 'vspdnup' == fname:
        return [unit, "'{}.vsp{}'".format(case, dnup), "'unknown'", \
                "'formatted'", 0]
    elif 'struct' == fname:
        return [unit, "'{}.struct'".format(case), "'old'", "'formatted'", 0]
    elif 'rotlm' == fname:
        return [unit, "'{}.rotlm'".format(case), "'unknown'", "'formatted'", 0]
    elif 'energysodum' == fname:
        if so == 'so':
            sodum = 'dum'
        else:
            sodum = dnup
        return [unit, "'{}.energy{}'".format(case, sodum), \
               "'unknown'", "'formatted'", 0]
    elif 'energyupdn' == fname:
        return [unit, "'{}.energy{}{}{}'".format(case, so, updn, para), \
                "'unknown'", "'formatted'", 0]
    elif 'energydnup' == fname:
        return [unit, "'{}.energy{}{}{}'".format(case, so, dnup, para), \
                "'unknown'", "'formatted'", 0]
    elif 'clmval' == fname:
        return [unit, "'{}.clmval{}'".format(case, updn), "'unknown'", \
                "'formatted'", 0]
    elif 'recprlist' == fname:
        return [unit, "'{}.recprlist'".format(case), "'unknown'", \
                "'formatted'", 9000]
    elif 'scf2updn' == fname:
        return [unit, "'{}.scf2{}'".format(case, updn), \
                "'unknown'", "'formatted'", 0]
    elif 'normupdn' == fname:
        if so == "so" and updn == "":
            _updn = "up"
        else:
            _updn = updn
        return [unit, "'{}.norm{}{}{}'".format(case, so, _updn, para), \
                "'unknown'", "'formatted'", 0]
    elif 'normdnup' == fname:
        return [unit, "'{}.norm{}{}{}'".format(case, so, dnup, para), \
                "'unknown'", "'formatted'", 0]
    else:
        raise ValueError('No matching file name {}!'.format(fname))


def fcreate_def_gwien(case, scratch='.', so='', para='', idmf='1', cmplx='',
        _band='', updn='', dnup='dn'):
    '''create gwien1/2.def file.
    '''
    fdef = open('gwien{}{}.def'.format(idmf,updn), 'w')
    if idmf == '1':
        fname_list = ['in2', 'inso', 'indmfl', 'outputdmfupdn', \
                'in1c', 'vectorupdn', 'vectordnup', 'klist', \
                'kgen', 'vspupdn', 'vspdnup', 'struct', \
                'rotlm', 'energydnup', 'energyupdn', 'normupdn', \
                'normdnup']
        unit_list = [3, 4, 5, 6, \
                7, 9, 10, 13, \
                14, 18, 19, 20, \
                22, 59, 60, 12, \
                11]
    elif idmf == '2':
        fname_list = ['in1c', 'inso', 'in2', 'outputdmfupdn', 'indmfl', \
                'clmval', 'vectorupdn', 'vectordnup', 'recprlist', 'kgen', \
                'vspupdn', 'struct', 'scf2updn', 'rotlm', 'energyupdn', \
                'normupdn', 'normdnup']
        unit_list = [3, 4, 5, 6, 7, \
                8, 9, 10, 13, 14, \
                18, 20, 21, 22, 30, \
                12, 11]

    for fname, unit in zip(fname_list, unit_list):
        fdef.write("{:3d}, {:<15s}, {:<10s}, {:<13s}, {:<4d}\n".format(\
                *get_file_info(fname, unit, idmf, case, scratch, so, \
                para, cmplx, _band, updn, dnup)))
    fdef.close()


def onestep(fday, case, exec_name, w_root, para="", so="", \
        band=None, updn=None):
    '''wien2k steps.
    '''
    time_start = time.strftime("%H:%M:%S")
    cmd = ['{}/x'.format(w_root), exec_name]
    if exec_name == "dstart":
        cmd.append("-lcore")
    cmd += ['-f', case]
    if para != "":
        cmd.append(para)
    if band == '-band':
        cmd.append(band)
        if not os.path.isfile('EFLDA.INP'):
            shutil.copy2('EFLDA.OUT', 'EFLDA.INP')
    if updn in ["-up", "-dn"]:
        cmd.append(updn)
    if so == "so":
        cmd.extend(["-c", "-so"])

    print(' '.join(x for x in cmd))
    process = subprocess.run(cmd, check=True, capture_output=True)
    fday.write('>{:<10s} ({}) {}\n'.format(exec_name, time_start, \
            process.stdout[:-1].decode()))
    fday.flush()
    for f in glob.glob('{}.error*'.format(exec_name)):
        if os.path.getsize(f) > 0:
            print('error in {} from file: {}'.format(
                    f, open(f, 'r').readlines()))
            sys.exit(1)


def gonestep(
        fday,
        exec_name,
        mpi,
        updn="",
        gpath=os.environ['WIEN_GUTZ_ROOT'],
        args=None,
        ):
    '''gwien1, CyGutz and gwien2 steps.
    '''
    time_start = time.strftime("%H:%M:%S")
    short_name = exec_name.split("/")[-1]
    with open(':log', 'a') as f:
        f.write('{}>   {}\n'.format(time.strftime("%a %b %d %H:%M:%S %Z %Y"), \
                short_name))

    cmd = ['/usr/bin/time']
    if mpi != '':
        cmd.extend(mpi)
    cmd.append('{}'.format(exec_name))
    if 'gwien' in short_name:
        cmd.append('{}{}.def'.format(short_name, updn))
    if args is not None:
        cmd.extend(args)

    print(' '.join(x for x in cmd))
    process = subprocess.run(cmd, check=True, capture_output=True)
    with open('{}_info.out'.format(short_name), 'w') as f:
        f.write(process.stdout.decode())
    fday.write('>{:<10s} ({}) {}\n'.format(short_name, time_start, \
            process.stderr.splitlines()[-2].decode()))
    fday.flush()
    for f in glob.glob('{}.error*'.format(short_name)):
        if os.path.getsize(f) > 0:
            print('error in {} from file: {}'.format(
                    f, open(f, 'r').readlines()))
            sys.exit(1)


def get_file_content(fname):
    if os.path.exists(fname):
        data = '\n------- {} --------\n'.format(fname)
        with open(fname, 'r') as f:
            data += f.read()
        return data
    else:
        return ''


def scf(case, spinpol):
    # scf file content
    if spinpol:
        f_list = ['{}.scf{}'.format(case, i) for i in ['0', \
                '1up', '1dn', 'so', '2up', '2dn', \
                '1s', '2s', 'cup', 'cdn']]
    else:
        f_list = ['{}.scf{}'.format(case, i) for i in ['0', \
                '1', 'so', '2', '1s', '2s', 'c']]

    data = ''.join(get_file_content(f) for f in f_list)

    with open('{}.scf'.format(case), 'a') as f:
        f.write(data)

    # files saved for mixing.
    if spinpol:
        f_list = ['clmsum', 'vspup', 'vspdn', 'vnsup', 'vnsdn', 'vrespsum',
                'clmdn', 'clmup']
    else:
        f_list = ['clmsum', 'vsp', 'vns', 'vrespsum']

    for i in f_list:
        name = '{}.{}'.format(case, i)
        if file_exists(name):
            shutil.copy2(name, '{}_old'.format(name))


def scfm(case):
    f_scf = '{}.scfm'.format(case)
    data = get_file_content(f_scf)

    with open('{}.scf'.format(case), 'a') as f:
        f.write(data)


def diff(fday, case, mix_dc, avg_dc, gskip):
    e_que = deque([], 2)
    with open('{}.scf'.format(case), 'r')  as f:
        for line in f:
            if ':DIS' in line:
                d_rho = float(line.split()[-1])
            if ':ENE' in line:
                e_que.append(float(line.split()[-1]))
    if len(e_que) == 2:
        d_etot = np.abs(e_que[1] - e_que[0])
    else:
        d_etot = 0.0

    dcv_err = 0.
    if not gskip:
        with h5py.File("GParam.h5", 'a') as f:
            ldc = f["/"].attrs["dc_mode"]
            if os.path.isfile("GLog.h5"):
                with h5py.File("GLog.h5", 'r') as fp:
                    nelf_list_inp = fp["/"].attrs["dc_nelf_list_inp"]
                    nelf_list_out = fp["/"].attrs["dc_nelf_list_out"]
                nelf_diff_list = nelf_list_out - nelf_list_inp
                nelf_list_mix = nelf_list_inp + mix_dc*nelf_diff_list
                if avg_dc:
                    ntot_avg = np.sum(nelf_list_mix)/nelf_list_mix.shape[0]
                    nelf_list_mix = [n/np.sum(n)*ntot_avg
                            for n in nelf_list_mix]
                if ldc == 12:
                    dcv_err = np.max(np.abs(nelf_diff_list))
                    f["/"].attrs["dc_nelf_list"] = nelf_list_mix

    fday.write(':ENERGY convergence: {}\n'.format(d_etot))
    fday.write(':CHARGE convergence: {}\n'.format(d_rho))
    fday.write(':VDC convergence: {}\n'.format(dcv_err))
    return d_rho, d_etot, dcv_err


def processes_convert(so, updn):
    if not file_exists('.processes'):
        print('.processes file not present. It must be a serial run.')
        return
    lines = open('.processes', "r").readlines()
    work = {}
    nkstart = 0
    for line in lines:
        data = line.split(':')
        if data[0].strip().isdigit():
            vecn = ["" for i in range(6)]
            i, nkp, nprc = map(int,data[::2])
            if not so:
                fdef = open('{}lapw1_{}.def'.format(updn,i), 'r')
                for line in fdef:
                    data = line.split(',')
                    data0 = int(data[0])
                    if data0 == 10 or data0 == 11:
                        data0 = data0 % 10
                        m = re.search('.*[\'|\"](.*)_(\d+)', data[1])
                        assert m is not None, 'vector file to macth ' + \
                                ' lapw1.def not found!'
                        name = m.group(1)
                        indx = m.group(2)
                        vecn[data0*2] = f'{name}_{indx}'
                        if updn == 'up':
                            name = name[:-2] + 'dn'
                        vecn[data0*2+1] = f'{name}_{indx}'
                fdef.close()
                prefix = vecn[0][:vecn[0].rfind(".")]
                vecn[4] = f"{prefix}.norm"
                vecn[5] = f"{prefix}.normdn"
            else:
                fdef = open('{}lapwso_{}.def'.format(updn,i), 'r')
                for line in fdef:
                    data = line.split(',')
                    if int(data[0])==42:
                        vecn[0]=data[1].split("'")[1]
                    elif int(data[0])==41:
                        vecn[1]=data[1].split("'")[1]
                    elif int(data[0])==52:
                        vecn[2]=data[1].split("'")[1]
                    elif int(data[0])==51:
                        vecn[3]=data[1].split("'")[1]
                    elif int(data[0])==46:
                        vecn[4]=data[1].split("'")[1]
                    elif int(data[0])==45:
                        vecn[5]=data[1].split("'")[1]
                fdef.close()

            if nprc in work:
                work[nprc].append((i, nkp, nkstart, vecn))
            else:
                work[nprc]=[(i, nkp, nkstart, vecn)]
            nkstart += nkp

    for prc in sorted(work.keys()):
        fo = open('_processes_{}'.format(prc-1), 'w')
        for (i, nkp, nkstart, vecn) in work[prc]:
            fo.write('{} {} {} "{}" "{}" "{}" "{}" "{}" "{}"\n'.format(\
                    i, nkp, nkstart, *vecn))


def run_gwien(
        startp='lapw0',
        endp="",
        cc=1.e-3,
        ec=1.e-5,
        vc=1.e-2,
        maxiter=100,
        mix_dc=0.2,
        band=False,
        avgdc=False,
        spinpol=False,
        spinorbit=False,
        dynamicslurm=False,
        dft=False,
        embhistory=False,
        rmethod="hybr",
        edinit=0,
        jac=0,
        nvals=None,
        ):
    parser = argparse.ArgumentParser(description=
            "Driver for wien2k + gutzwiller-slave-boson job.")
    parser.add_argument("-s", "--startp", type=str, default=startp,
            help="starting program (str)")
    parser.add_argument("-e", "--endp", type=str, default=endp,
            help="ending program (str)")
    parser.add_argument("-cc", "--cc", type=float, default=cc,
            help="charge density convergence criterion (float)")
    parser.add_argument("-ec", "--ec", type=float, default=ec,
            help="energy convergence criterion (float)")
    parser.add_argument("-vc", "--vc", type=float, default=vc,
            help="dc potential convergence criterion (float)")
    parser.add_argument("-n", "--maxiter", type=int, default=maxiter,
            help="maximal iterations (int)")
    parser.add_argument("-amix", "--mix_dc", type=float, default=mix_dc,
            help="dc mixing parameter (float)")
    parser.add_argument("-band", action="store_true", default=band,
            help="band structure calculation")
    parser.add_argument("-avgdc", action="store_true", default=avgdc,
            help="average dc across all the correlated sites.")
    parser.add_argument("-sp", "--spinpol", action="store_true",
            default=spinpol,
            help="perform spin-polarized calculation")
    parser.add_argument("-so", "--spinorbit", action="store_true",
            default=spinorbit,
            help="include spin-orbit interaction")
    parser.add_argument("-ds", "--dynamicslurm", action="store_true",
            default=dynamicslurm,
            help="dynamically switch off slurm")
    parser.add_argument("-dft", action="store_true", default=dft,
            help="pure dft calculation")
    parser.add_argument("--embhistory", action="store_true",
            default=embhistory,
            help="record embedding hamiltonian history")
    parser.add_argument("-r", "--rmethod", type=str, default=rmethod,
            help="nleq solver: hybr, excitingmixing, df-sane (str)")
    parser.add_argument("-edinit", type=int, default=edinit,
            help="ed with init v from previous run (1) or not (0) (int)")
    parser.add_argument("-j","--jac", type=int, default=jac,
            help="jacobian method: 0->no, 1->analytical, -1: numerical.")
    parser.add_argument("--nvals", type=str, default=nvals,
            help="valences list: e.g. '5 5 5' for 3 impurities.")
    args = parser.parse_args()
    if args.jac == -1:
        from pygrisb.run.cygutz import get_jacobian_numerical
        args.jac = get_jacobian_numerical
    elif args.jac == 1:
        from pygrisb.run.cygutz import get_jacobian_analytical
        args.jac = get_jacobian_analytical
    else:
        args.jac = None
    if args.nvals is not None:
        args.nvals = list(map(int, args.nvals.split()))

    if args.band:
        band ='-band'
        _band = '_band'
        args.maxiter = 1
    else:
        band = _band = ""
    cygtail = get_cygtail()
    para = ''
    _para = ''
    if file_exists('.machines') :
        para = ' -p'
        _para = '_x'

    clean_list = glob.glob('*.scf*') + glob.glob('*.error*') + \
            glob.glob('*.outputdmf?.*')
    for f in clean_list:
        os.remove(f)

    struct_file = glob.glob('*.struct')
    if len(struct_file) != 1:
        raise ValueError('{} struct files present while only one must exist!'. \
                format(len(struct_file)))
    w_case = struct_file[0].split('.')[0]
    w_root = os.environ['WIENROOT']
    w_scratch = os.environ['SCRATCH']
    g_root = os.environ['WIEN_GUTZ_ROOT']

    # infomation file
    fday = open(w_case + '.dayfile', 'w')
    fday.write('Calculating {} in {} \non {} with PID {}\n'.format(\
            w_case, os.getcwd(), socket.gethostname(), os.getpid()))

    if os.path.isfile(w_case + '.inso') and \
            os.path.getsize(w_case + '.inso') > 0 and not args.spinorbit:
        msg = 'spin-orbit is off while {} file is present!'.format(\
                w_case + '.inso')+'\nappending -so to turn it on.'
        warnings.warn(msg)

    if args.spinorbit:
        so='so'
        cmplx = 'c'
    else:
        so = ''
        cmplx = ''

    # In addition, check in1c file
    if file_exists(w_case+'.in1c'):
        cmplx = 'c'

    f_mpi = 'mpi_prefix.dat'
    if os.path.isfile(f_mpi):
        with open(f_mpi, 'r') as f:
            lines = f.readlines()
            mpi = lines[0].split()
            if len(lines) == 1 or len(lines[1].split())==0:
                mpi2 = mpi
            else:
                mpi2 = lines[1].split()
        print('{} exists -- running in parallel mode.'.format(f_mpi))
        print(' '.join(x for x in mpi2)+" for gsolver.")
        print(' '.join(x for x in mpi)+" for others.")
    else:
        if para != '':
            raise ValueError('missing mpi_prefix.dat with .machines present!')
        mpi = mpi2 = ["mpirun", "-np", "1"]
        print('{} not available -- running in serial mode.'.format(f_mpi))

    if not args.dft:
        # create gwien1/2.def files
        if args.spinpol:
            fcreate_def_gwien(w_case, scratch=w_scratch, \
                    so=so, para=_para, \
                    idmf='1', cmplx=cmplx, _band=_band, \
                    updn="up", dnup='dn')
            if not args.spinorbit:
                fcreate_def_gwien(w_case, scratch=w_scratch, \
                        so=so, para=_para,
                        idmf='1', cmplx=cmplx, _band=_band, \
                                updn="dn", dnup='up')

            fcreate_def_gwien(w_case, scratch=w_scratch, \
                    so=so, para=_para, \
                    idmf='2', cmplx=cmplx, _band=_band, \
                    updn="up", dnup='dn')
            fcreate_def_gwien(w_case, scratch=w_scratch, \
                    so=so, para=_para, \
                    idmf='2', cmplx=cmplx, _band=_band, \
                    updn="dn", dnup='up')
        else:
            fcreate_def_gwien(w_case, scratch=w_scratch, \
                    so=so, para=_para, \
                    idmf='1', cmplx=cmplx, _band=_band)
            fcreate_def_gwien(w_case, scratch=w_scratch, \
                    so=so, para=_para, \
                    idmf='2', cmplx=cmplx, _band=_band)

    # save SLURM_ related environment variables
    slurm_envs = env.get_env_dict(key="SLURM_")

    if args.maxiter > 0:
        if not os.path.isfile('{}.clmsum'.format(w_case)):
            err_msg = 'no {}.clmsum file found--necessary for lapw0!'.\
                    format(w_case)
            print(err_msg)
            fday.print(err_msg+'\n')
            sys.exit(1)
        for f in glob.glob('*.broyd*'):
            os.remove(f)

        fday.write('   start at {} with {} \n   1/{} to go.\n'.format(
                time.asctime(), args.startp, args.maxiter))

    # Major charge density loop
    for icycle in range(args.maxiter):

        # unset slurm environment variables for wien2k type run
        if args.dynamicslurm:
            env.unset_environ(slurm_envs)

        if icycle > 0 or (args.startp in 'lapw0'):
            onestep(fday, w_case, 'lapw0', w_root, para=para)
        if icycle > 0 or (args.startp in 'lapw0 lapw1'):
            if args.spinpol:
                onestep(fday, w_case, 'lapw1', w_root, para=para, band=band, \
                        updn="-up")
                onestep(fday, w_case, 'lapw1', w_root, para=para, band=band, \
                        updn="-dn")
            else:
                onestep(fday, w_case, 'lapw1', w_root, para=para, band=band)
        if (icycle > 0 or (args.startp in 'lapw0 lapw1 lapwso')) \
                and args.spinorbit:
            if args.spinpol:
                onestep(fday, w_case, 'lapwso', w_root, para=para, \
                        band=band, updn="-up")
            else:
                onestep(fday, w_case, 'lapwso', w_root, para=para, band=band)

        if icycle==0 and para != '':
            if args.spinpol:
                processes_convert(args.spinorbit, updn="up")
            else:
                processes_convert(args.spinorbit, updn="")

        #set slurm environment variables for mpi run
        if args.dynamicslurm:
            env.set_environ(slurm_envs)

        if args.dft:
            # run dft only
            if icycle > 0 or (args.startp in 'lapw0 lapw1 lapwso lapw2'):
                if args.spinpol:
                    onestep(fday, w_case, 'lapw2', w_root, para=para, \
                            updn="-up", so=so)
                    onestep(fday, w_case, 'lapw2', w_root, para=para, \
                            updn="-dn", so=so)
                else:
                    onestep(fday, w_case, 'lapw2', w_root, para=para, so=so)
            gerr = 0.
        else:
            if icycle > 0 or (args.startp in 'lapw0 lapw1 lapwso gwien1'):
                if args.spinpol:
                    gonestep(fday, g_root+'/gwien1', mpi, updn="up")
                    if not args.spinorbit:
                        gonestep(fday, g_root+'/gwien1', mpi, updn="dn")
                else:
                    gonestep(fday, g_root+'/gwien1', mpi)
            if args.endp == 'gwien1':
                sys.exit(0)

            gonestep(fday, f"{g_root}/CyGInit{cygtail}", mpi)
            if band == '-band' or args.endp == 'CyGBand':
                gonestep(fday, f"{g_root}/CyGBand{cygtail}", mpi,
                        args=["-w", "1"])
                sys.exit(0)

            # solve the Gutzwiller nonlinear equations.
            verr = gsolve_nleqns(
                    gfun_solver,
                    jac=args.jac,
                    method=args.rmethod,
                    args={
                        "mpi":mpi,
                        "mpi2":mpi2,
                        "path":g_root,
                        "edinit":args.edinit,
                        "cygtail":cygtail,
                        "nval_list":args.nvals,
                        "save_hemb_history":args.embhistory,
                        },
                    )
            gerr = np.max(np.abs(verr))

            cmd = mpi + [f"{g_root}/CyGFinal{cygtail}", "-r", "1"]
            print(" ".join(cmd))
            subprocess.run(cmd, check=True)

            if args.endp == 'CyGutz':
                sys.exit(0)

            # copy for monitoring
            shutil.copyfile('Gutz.log', 'old_Gutz.log')

            if args.spinpol:
                gonestep(fday, g_root+'/gwien2', mpi, updn="up")
                gonestep(fday, g_root+'/gwien2', mpi, updn="dn")
            else:
                gonestep(fday, g_root+'/gwien2', mpi)

            if args.endp == 'gwien2':
                sys.exit(0)

        # unset slurm environment variables for wien2k type run
        if args.dynamicslurm:
            env.unset_environ(slurm_envs)

        if args.spinpol:
            onestep(fday, w_case, 'lcore', w_root, para='', updn="-up")
            onestep(fday, w_case, 'lcore', w_root, para='', updn="-dn")
            if os.path.isfile('.lcore'):
                onestep(fday, w_case, 'dstart', w_root, para='', updn="-up")
                onestep(fday, w_case, 'dstart', w_root, para='', updn="-dn")
                os.remove(f"{w_case}.clmcorup")
                os.remove(f"{w_case}.clmcordn")
        else:
            onestep(fday, w_case, 'lcore', w_root, para='')
            if os.path.isfile('.lcore'):
                onestep(fday, w_case, 'dstart', w_root, para='')
                os.remove(f"{w_case}.clmcor")

        scf(w_case, args.spinpol)
        onestep(fday, w_case, 'mixer', w_root, para='')
        scfm(w_case)
        drho, dene, dvdc = diff(fday, w_case, args.mix_dc,
                args.avgdc, args.dft)

        print(('dc={:.1e}, cc={:.1e} -> {:.0e}, ec={:.1e} ' + \
                '-> {:.0e}, gc={:.1e} icycle={}').format(
                dvdc, drho, args.cc, dene, args.ec, gerr, icycle))
        if drho < args.cc and dene < args.ec and dvdc < args.vc:
            sys.exit(0)


def batch_init_ga(dir_template='./template'):
    '''Loop over all the directories to initialize CyGutz calculations
     -- actually, since the CyGutz input files remain the same for different
    volumes, it simply copy the input files in template directory to
    each folder.
    '''
    cwd = os.getcwd()+'/'
    for dname in get_dir_list():
        print(f"copy grisb init file to {dname}.")
        os.chdir(dname+'/case')
        shutil.copy(cwd+'/'+dir_template+'/ginit.json', './')
        shutil.copy(cwd+'/'+dir_template+'/GParam.h5', './')
        shutil.copy(cwd+'/'+dir_template+'/case.indmfl', './')
        os.chdir(cwd)


def batch_init_mott(dir_template='./template'):
    '''Loop over all the directories to initialize CyGutz-Mott calculations
     -- actually, since the CyGutz input files remain the same for different
    volumes, it simply copy the input files in template directory to
    each folder.
    '''
    cwd = os.getcwd()+'/'
    for dname in get_dir_list():
        os.chdir(dname+'/case')
        shutil.copy(cwd+'/'+dir_template+'/GMOTT.h5', './')
        os.chdir(cwd)


def batch_modify_ga_setup(args, nproc=1):
    '''Loop over all the directories to modify CyGutz set up file.
    '''
    cwd = os.getcwd()+'/'
    cmd = [os.environ['WIEN_GUTZ_ROOT']+'/switch_gparam.py'] + args
    if '-p' in sys.argv:
        nproc = int(sys.argv[sys.argv.index('-p')+1])
    for i, dname in enumerate(get_dir_list()):
        os.chdir(dname+'/case')
        proc = subprocess.Popen(cmd)
        os.chdir(cwd)
        if (i+1) % nproc == 0:
            proc.communicate()


def batch_job_slurm(jname="dftg",
        dir_template='./template',
        dir_work='./',
        subdir='case',
        jobfile="job.slurm"):
    '''copy template/job.slurm file to each working directory and submit jobs.
    '''
    cwd = os.getcwd()+'/'
    if '-n' in sys.argv:
        jname = sys.argv[sys.argv.index('-n')+1]

    if '-w' in sys.argv:
        dir_work = sys.argv[sys.argv.index('-w')+1]

    if '-j' in sys.argv:
        jobfile = sys.argv[sys.argv.index('-j')+1]

    # command to submit job.
    cmd_s = ['sbatch', './job.slurm']

    for dname in get_dir_list():
        os.chdir(f"{dname}/{subdir}/{dir_work}")
        print(f"work in dir {dname}.")
        # get job.slurm file
        with open(f"{cwd}/{dir_template}/{jobfile}", 'r') as fin:
            with open('./job.slurm', 'w') as fout:
                for line in fin:
                    fout.write(line.replace('VV', dname). \
                            replace('UJ', jname))
        # submit job
        subprocess.run(cmd_s, check=True)
        os.chdir(cwd)


def run_ga(nproc=1):
    '''Loop over all the directories to run_ga using nproc processors.
    '''
    cmd = [os.environ['WIEN_GUTZ_ROOT']+'/run_gwien.py']
    if '-so' in sys.argv:
        cmd += ['-so']
    if '-sp' in sys.argv:
        cmd += ['-sp']
    cwd = os.getcwd()+'/'
    if '-p' in sys.argv:
        nproc = int(sys.argv[sys.argv.index('-p')+1])
    for i,dname in enumerate(get_dir_list()):
        os.chdir(dname+'/case')
        proc = subprocess.Popen(cmd)
        os.chdir(cwd)
        if (i+1) % nproc == 0:
            proc.communicate()


def batch_gsave(sdir='ldag', args=['-f']):
    '''Loop over all the directories to save_lapw.
    '''
    if '-d' in sys.argv:
        sdir = sys.argv[sys.argv.index('-d')+1]
    cmd = [os.environ['WIEN_GUTZ_ROOT']+'/save_ldag', '-d'] + [sdir] + args
    cwd = os.getcwd()+'/'
    for dname in get_dir_list():
        os.chdir(dname+'/case')
        subprocess.run(cmd)
        os.chdir(cwd)



if __name__=="__main__":
    run_gwien()
