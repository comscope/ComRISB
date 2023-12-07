import numpy as np
import os, shutil, subprocess, glob, h5py, pickle, time, argparse, json
from tabulate import tabulate
from pygrisb.run.cygutz import driver_cygutz
from pygrisb.basic.units import rydberg_to_ev


def open_h_log(control):
    control['h_log'] = open('./cmd.log', 'w')
    control['h_log'].write(
            "\n*********************************\n" +
            "             ComRISB\n" +
            "*********************************\n\n")

# init_comrisb:
#   - reads in  comrisb.ini
#   - loads ginit.json into  control['impurity_problem'] and control['impurity_solution']
#   - parses its arguments: 
#       -s or --startp sets control['start_prog']
#       -e or --endp sets control['end_prog']
#       -c or --continuous sets control['restart']
def init_comrisb():
    vglobl = {}
    vlocal = {}
    exec(open('comrisb.ini').read(), vglobl, vlocal)
    control = vlocal['control']
    wan_hmat = vlocal['wan_hmat']

    control['name'] = 'control'
    wan_hmat['name'] = 'wan_hmat'
    imp = {'name': 'imp'}

    # timing
    control['t_rspflapw'] = 0.
    control['t_comwann'] = 0.
    control['t_interface'] = 0.
    control['t_cygutz'] = 0.
    control['t_total'] = 0.

    # directory
    control['top_dir'] = os.path.abspath('./')
    control['wannier_directory'] = os.path.abspath('./wannier')
    control['lowh_directory'] = os.path.abspath('./lowh')

    # compatibility
    if 'initial_dft_dir' in control:
        control['initial_lattice_dir'] = control['initial_dft_dir']

    if control['method'] == "lqsgw+risb":
        control['initial_lattice_dir'] = os.path.abspath(
                control.get('initial_lattice_dir', '../lqsgw'))
        control['lattice_directory'] = control['initial_lattice_dir']
    else:
        control['initial_lattice_dir'] = os.path.abspath(
                control.get('initial_lattice_dir', '../dft'))
        control['lattice_directory'] = os.path.abspath('./lattice')

    control['allfile'] = find_allfile(control['initial_lattice_dir'])

    # log file
    open_h_log(control)

    if 'comsuitedir' not in control:
        control['comsuitedir'] = os.environ.get('COMRISB_BIN', None)
    if control['comsuitedir'] is None:
        raise OSError("please set comsuitedir either in environ var or "+\
                "comrisb.ini!")

    # convergence table
    control['conv_table'] = []

    # load ginit.json into  control['impurity_problem'] and control['impurity_solution']
    ginit = json.load(open(f"{control['lowh_directory']}/ginit.json", "r"))
    control["spin_orbit"] = ginit["gchoice"]["iso"] == 2
    # get impurity_problem
    impurities = []
    for i, symbol in enumerate(ginit["struct"]["symbols"]):
        if symbol in ginit["gchoice"]["unique_cor_symbol_list"]:
            iun = ginit["gchoice"]["unique_cor_symbol_list"].index(symbol)
            # one-based here
            impurities.append([i+1, ginit["gchoice"]["unique_cor_orb_list"] \
                    [iun]])
    control['impurity_problem'] = impurities
    control['impurity_solution'] = control.get('impurity_solution', 1)
    
    # fermi-dirac, only option in rspflapw.
    control['ismear'] = -1
    
    # smearing factor. default 300K to eV   
    with open(f"{control['initial_lattice_dir']}/ini", "r") as f:
        for line in f:
            if 'temperature' in line:
                delta = float(line.replace('=', '').split()[1])
    # in unit of ev
    control['delta'] = delta/11604.52500617

    control['max_iter_num_outer'] = control.get('max_iter_num_outer', 20)
    control['h_log'].write('allfile = {}\n'.format(control['allfile']))
    control["mpi_prefix"] = control.get('mpi_prefix', "")
    control["mpi_prefix_wannier"] = control.get('mpi_prefix_wannier', \
            control["mpi_prefix"])
    control["mpi_prefix_lattice"] = control.get('mpi_prefix_lattice', \
            control["mpi_prefix"])
    control['restart'] = control.get('restart', True)
    control['cc'] = control.get('cc', 1.e-8)
    control['ec'] = control.get('ec', 1.e-4)

# The arguments are parsed:
#   -s or --startp sets control['start_prog']
#   -e or --endp sets control['end_prog']
#   -c or --continuous sets control['restart']
    parser = argparse.ArgumentParser(description="Driver for comrisb job.")
    parser.add_argument("-s", "--startp", type=str, default='initcopy',
            help="starting program (str) [initcopy, wannier, gwann, dft]")
    parser.add_argument("-e", "--endp", type=str, default='',
            help="ending program (str)")
    parser.add_argument("-c", "--continuous", action="store_true",
            default=control['restart'],
            help="continuous calculation")
    args = parser.parse_args()
    control['restart'] = args.continuous
    control['start_prog'] = args.startp
    control['end_prog'] = args.endp

    # in wan_hmat
    check_key_in_string('froz_win_max', wan_hmat)
    wan_hmat['dis_win_max'] = wan_hmat.get('dis_win_max',
            wan_hmat['froz_win_max']+40.0)
    wan_hmat['froz_win_min'] = wan_hmat.get('froz_win_min', None)
    control['proj_win_max'] = control.get('proj_win_max',
            wan_hmat['froz_win_max'])
    wan_hmat['num_iter'] = wan_hmat.get('num_iter', 0)
    wan_hmat['dis_num_iter'] = wan_hmat.get('dis_num_iter', 100)
    wan_hmat['cut_low']=wan_hmat.get('cut_low', 0.4)
    wan_hmat['cut_froz']=wan_hmat.get('cut_froz', 0.1)
    wan_hmat['cut_total']=wan_hmat.get('cut_total', 0.0)

    control['h_log'].write('top_dir {}\n'.format(control['top_dir']))
    control['h_log'].write('lattice_directory {}\n'.
            format(control['lattice_directory']))
    control['h_log'].write('wannier_directory {}\n'.
            format(control['wannier_directory']))
    control['h_log'].write('lowh_directory {}\n'
            .format(control['lowh_directory']))

    return control, wan_hmat, imp

# find_impurity_wan sets up control['impurity_wan'], based on the 
#   contents of wan_hmat['basis']
#   wan_hmat['basis'] is populated by read_wan_hmat_basis from wannier.inip, and it
#       it is used only here and in get_locrot_list. 
#  control['impurity_wan'] is used only to populate params["corbs_list"] ,
#       which is passed to the gutzwiller solver.  
#   find_impurity_wan is implemented only for d,f, and not s,p shells.
def find_impurity_wan(control, wan_hmat, tol=1.e-6):
    num_wann = np.shape(wan_hmat['basis'])[0]
    control['impurity_wan'] = []
    l_ch = {'f': 3, 'd': 2}
    for ip in range(np.shape(control['impurity_problem'])[0]):
        l = l_ch[control['impurity_problem'][ip][1].lower()]
        if control['spin_orbit']:
            control['impurity_wan'].append([0]*(2*l + 1)*2)
            for iwan in range(num_wann):
                # one-based index
                if wan_hmat['basis'][iwan]['atom'] == \
                        control['impurity_problem'][ip][0] \
                        and wan_hmat['basis'][iwan]['l'] == l:
                    s = wan_hmat['basis'][iwan]['i']
                    ishft = int(s + 0.51)*int((l-0.5)*2+1.1)
                    idx = int(wan_hmat['basis'][iwan]['m'] + (l+s) + 0.1)
                    idx += ishft
                    control['impurity_wan'][ip][idx] = \
                            wan_hmat['basis'][iwan]['ind']
        else:
            control['impurity_wan'].append([0]*(2*l + 1))
            for iwan in range(num_wann):
                if wan_hmat['basis'][iwan]['atom'] == \
                        control['impurity_problem'][ip][0]\
                        and wan_hmat['basis'][iwan]['l'] == l:
                    idx = wan_hmat['basis'][iwan]['m']+l
                    control['impurity_wan'][ip][idx] = \
                            wan_hmat['basis'][iwan]['ind']
        if (control['impurity_wan'][ip].count(0) != 0):
            print(f"impurity {ip} orbital map: {control['impurity_wan'][ip]}")
            raise ValueError('Something wrong in find_impurity_wan.')


def initial_file_directory_setup(control):
    for tempdir in [control['wannier_directory'],
            control['lattice_directory'],
            control['lowh_directory']]:
        os.makedirs(tempdir, exist_ok=True)


def initial_lattice_directory_setup(control):
    os.chdir(control['lattice_directory'])
    files = glob.iglob(control['initial_lattice_dir']+"/*.rst")
    for filename in files:
        shutil.copy(filename, './')
    files = glob.iglob(control['initial_lattice_dir']+"/*el_density")
    for filename in files:
        shutil.copy(filename, './')
    if os.path.exists(control['initial_lattice_dir']+'/kpath'):
        shutil.copy(control['initial_lattice_dir']+'/kpath', './')
    shutil.copy(control['initial_lattice_dir']+'/wannier_win.dat', './')
    shutil.copy(control['initial_lattice_dir']+'/ini', './')

    iter_string = '_'+str(control['iter_num_outer'])
    shutil.copy(control['initial_lattice_dir']+'/'+control['allfile']+'.out',
            control['allfile']+iter_string+'.out')

    control['h_log'].write("initial dft directory setup done.\n")
    os.chdir(control['top_dir'])


# create_comwann_ini creates comwann.ini.  If not doing a qsgw calculation,
#   it uses information about the band energies and fermi level to decide
#   the values of wan_hmat['froz_win_min'] and wan_hmat['dis_win_min'] .
def create_comwann_ini(control, wan_hmat):

    # Begin: This block gets el, eh, and efermi, and is not used if doing a qsgw calculation.

    # get ncoremax, a count of bands, from wannier_win.dat.    
    line = open(control['lattice_directory']+ \
            "/wannier_win.dat","rb").readline()    
    ncoremax = int(line.split()[-1])
    
    # read in ebands and efermi
    with h5py.File(f"{control['lattice_directory']}/{control['allfile']}.rst",
            "r") as f:
        if "qsgw" in control["method"]:
            ebands = f[f"/{control['allfile']}/e_qp/qp"][()]
            efermi = f["/chemical_potential/qp/chem_pot"][0]
        else:
            ebands = f[f"/{control['allfile']}/e_bnd/dft"][()]
            efermi = f["/chemical_potential/dft/chem_pot"][0]
            
    # Using ebands, decide the highest and lowest band energies, el and eh.
    el, eh = np.max(ebands[0, :, ncoremax-1]), np.min(ebands[0, :, ncoremax])

    # End: This block gets el, eh, and efermi, and is not used if doing a qsgw calculation.
    
    # Decide  the values of wan_hmat['froz_win_min'] and wan_hmat['dis_win_min']
    #   If not doing qsgw, these values are decided by using el, eh, and efermi.
    if "qsgw" in control["method"]:
        # the qp energies not gapped some how.
        # can be specified in comrisb.ini.
        if wan_hmat['froz_win_min'] is None:
            wan_hmat['froz_win_min'] = -10
        if 'dis_win_min' not in wan_hmat:
            wan_hmat['dis_win_min'] = wan_hmat['froz_win_min'] - 40
    else:
        if eh <= el:
            raise ValueError("no gap between core and valence states.")
        emin = ((el + eh)/2 - efermi)*rydberg_to_ev
        # prefixed win_min
        wan_hmat['dis_win_min'] = emin
        # can be specified in comrisb.ini.
        if wan_hmat['froz_win_min'] is None:
            wan_hmat['froz_win_min'] = emin

    with open('comwann.ini', 'w') as f:
        f.write(control['lattice_directory']+'\n')
        if "qsgw" in control["method"]:
            f.write('qp\n')
        else:
            f.write('dft\n')
        f.write(str(wan_hmat['dis_win_max'])+'\n')
        f.write(str(wan_hmat['dis_win_min'])+'\n')
        f.write(str(wan_hmat['froz_win_max'])+'\n')
        f.write(str(wan_hmat['froz_win_min'])+'\n')
        f.write(str(wan_hmat['num_iter'])+'\n')
        f.write(str(wan_hmat['dis_num_iter'])+'\n')
        f.write('0\n') # wan_hmat['write_wan'] is ignored and hard-coded to 0.
        f.write(str(wan_hmat['cut_low'])+'\n')
        f.write(str(wan_hmat['cut_froz'])+'\n')
        f.write(str(wan_hmat['cut_total'])+'\n')
#        f.write(str(wan_hmat['radfac'])+'\n')
#        wan_hmat['froz_win_max_fac']=wan_hmat.get('froz_win_max_fac', 0.7)    
#        f.write(str(wan_hmat['froz_win_max_fac'])+'\n')        
        


# read_wan_hmat_basis reads wannier.inip into the returned variable
def read_wan_hmat_basis(control):
    # in the wannier directory
    inip = np.loadtxt(control['wannier_directory']+'/wannier.inip')
    basis_info = []

    if (control['spin_orbit']):
        for ii in range(np.shape(inip)[0]):
            basis_info.append({'atom': int(inip[ii, 0]), \
                    'l': int(inip[ii, 1]), 'i': inip[ii, 2], \
                    'm': inip[ii, 3], 'xaxis': inip[ii, 4:7], \
                    'zaxis': inip[ii, 7:10], 'ind': ii})
    else:
        for ii in range(np.shape(inip)[0]):
            basis_info.append({'atom': int(inip[ii, 0]), \
                    'l': int(inip[ii, 1]), 'm': int(inip[ii, 2]), \
                    'xaxis': inip[ii, 3:6], 'zaxis': inip[ii, 6:9], \
                    'ind': ii})
    control['h_log'].write(\
            'reading wannier.inip to get basis information.\n')
    return basis_info


def check_key_in_string(key, dictionary):
    if key not in dictionary:
        raise ValueError('missing \''+key+'\' in '+dictionary['name'])


def labeling_file(filename, iter_string):
    dirname = os.path.abspath(os.path.dirname(filename))
    filenameonly = os.path.basename(filename)
    temp = filenameonly.split('.')
    shutil.copy(dirname+'/'+filenameonly, dirname+"/" +
                '.'.join(temp[0:-1])+iter_string+'.'+temp[-1])

# write_conv gathers info from rspflapw's out file, from GLog.h5,
#   and from GIter.h5, and writes it out in convergence.log
def write_conv(control):
    os.chdir(control['lattice_directory'])
    iter_string = '_' + str(control['iter_num_outer'])
    
    # read in delta_rho and ETOT from rspflapw's out file
    with open('./dft'+iter_string+'.out', "r") as f:
        for line in f:
            if "charge density" in line:
                line = line.split()
                delta_rho=float(line[3])
            elif "ETOT" in line:
                line = line.split()
                etot = float(line[4])
                
    # read in mu=efermi, egamma_dc, ebands, ebands_bare, RMAT from GLog.h5
    with h5py.File(control['lowh_directory']+"/GLog.h5", "r") as f:
        mu = f["/"].attrs["efermi"]
        egamma_dc = f["/"].attrs["egamma_dc"]
        ebands = f["/"].attrs["ebands"]
        ebands_bare = f["/"].attrs["ebands_bare"]
        rmat2 = f["/"].attrs["RMAT"].swapaxes(1,2)

    # calculate etot
    # save control['iter_num_outer'], delta_rho, etot in convergence.log
    etot += (egamma_dc + sum(ebands-ebands_bare))/rydberg_to_ev   
    control['conv_table'].append(
            [control['iter_num_outer'],
            delta_rho,
            etot])

    # read in v_err from GIter.h5
    # save mu=efermi, err_risb in convergence.log
    with h5py.File(control['lowh_directory']+"/GIter.h5", "r") as f:
        v_err = f["/v_err"][()]
        err_risb = v_err[abs(v_err).argmax()]
        control['conv_table'][-1].extend([mu, err_risb])
        
    # check quasiparticle weight
    # zmat = rmat^\dagger * rmat
    # save zmin in convergence.log
    zmin = 100.
    for rmat in rmat2:
        zmat = rmat.T.conj().dot(rmat)
        w = np.linalg.eigvalsh(zmat)
        zmin = min(zmin, np.amin(w))        
    control['conv_table'][-1].append(zmin)
    
    # save out 
    os.chdir(control['top_dir'])
    with open("convergence.log", "w") as f:
        f.write(tabulate(control['conv_table'], \
                headers=['i_outer', 'delta_rho', 'etot', "mu", \
                'err_risb', "min_z"], \
                numalign="right",  floatfmt=".8f"))
        f.write('\n')

# check_wannier_function_input creates comwann.ini
def check_wannier_function_input(control, wan_hmat):
    os.chdir(control['wannier_directory'])
    
    # create_comwann_ini creates comwann.ini.  If not doing a qsgw calculation,
#   it uses information about the band energies and fermi level to decide
#   the values of wan_hmat['froz_win_min'] and wan_hmat['dis_win_min'] .
    create_comwann_ini(control, wan_hmat)
    
    if os.path.exists(control['top_dir']+'/local_axis.dat'):
        shutil.copy(control['top_dir']+'/local_axis.dat', './')
        
    os.chdir(control['top_dir'])

# run_dft runs rspflapw
def run_dft(control):
    
    os.chdir(control['lattice_directory'])
    
    iter_string = '_'+str(control['iter_num_outer'])
    cmd = control['mpi_prefix_lattice'] + ' ' + \
            control['comsuitedir'] + "/rspflapw.exe"
    control['h_log'].write(cmd+"\n")
    hlog_time(control['h_log'], "rspflapw start")
    t_start = time.time()
    with open(control['lattice_directory']+'/dft.out', 'w') as logfile:
        subprocess.run(cmd,
                shell=True,
                stdout=logfile,
                stderr=logfile,
                check=True,
                )
    hlog_time(control['h_log'], "rspflapw end")
    t_end = time.time()
    control['t_rspflapw'] += t_end-t_start
    allfile=control['allfile']
    labeling_file('./'+allfile+'.out',iter_string)
    shutil.move('./dft.out', './dft'+iter_string+'.out')
    control['h_log'].write("dft calculation done.\n")
    os.chdir(control['top_dir'])

# prepare_dft_input makes sure that  wannier_den_matrix.dat is in 
#       the right place for rspflapw to use it.
def prepare_dft_input(control):
    shutil.copy(control['lowh_directory']+"/wannier_den_matrix.dat",
            control['lattice_directory'])
    control['h_log'].write("prepare_dft_input done.\n")


def hlog_time(f, prestr, endl=""):
    f.write("{} at {} {}".format(prestr, time.strftime("%d/%m/%Y %H:%M:%S"), \
            endl))
    f.flush()


# wannier_run:
#   - runs ComWann, but only if the fullrun argument is True
#   - reads wannier.inip, and uses it to set up control['impurity_wan']
#       control['impurity_wan'] is used only to populate params["corbs_list"] ,
#           which is passed to the gutzwiller solver.  
# wannier_run is implemented only for d,f, and not s,p shells.   
def wannier_run(control,wan_hmat,fullrun=True):
    os.chdir(control['wannier_directory'])
    if fullrun:
        cmd = control['mpi_prefix_wannier']+' '+\
                control['comsuitedir']+"/ComWann"
        control['h_log'].write(cmd+"\n")

        hlog_time(control['h_log'], "comwann start")
        t_start = time.time()
        with open(control['wannier_directory']+'/comwann.out', 'w') \
                as logfile:
            subprocess.run(cmd,
                    shell=True,
                    stdout=logfile,
                    stderr=logfile,
                    check=True,
                    )
        hlog_time(control['h_log'], "end", endl="\n")
        t_end = time.time()
        control['t_comwann'] += t_end-t_start

        iter_string = '_' + str(control['iter_num_outer'])
        shutil.move('./wannier.wout', './wannier'+iter_string+'.wout')

    # read_wan_hmat_basis reads wannier.inip into wan_hmat['basis']
    wan_hmat['basis'] = read_wan_hmat_basis(control)
    
    # find_impurity_wan sets up control['impurity_wan'], based on the 
    #   contents of wan_hmat['basis']
    #   wan_hmat['basis'] is populated by read_wan_hmat_basis from wannier.inip, and it
    #       it is used only here and in get_locrot_list. 
    #  control['impurity_wan'] is used only to populate params["corbs_list"] ,
    #       which is passed to the gutzwiller solver.  
    #   find_impurity_wan is implemented only for d,f, and not s,p shells.
    find_impurity_wan(control, wan_hmat)
    
    control['h_log'].write("control['impurity_wan']: {}\n".format(\
            control['impurity_wan']))
    os.chdir(control['top_dir'])

# get_locrot_list calculates a rotation matrix based on the xaxis and yaxis
#   supplied by wannier.inip
#   wan_hmat['basis'] is populated by read_wan_hmat_basis from wannier.inip, and it
#       it is used only here and in find_impurity_wan. 
def get_locrot_list(wan_hmat):
    iatm = -1
    lrot_list = []
    for basis in wan_hmat["basis"]:
        if basis["atom"] != iatm:
            x = np.asarray([float(e) for e in basis["xaxis"]])
            z = np.asarray([float(e) for e in basis["zaxis"]])
            y = np.cross(z, x)
            lrot_list.append(np.asarray([x,y,z]).T)
            iatm = basis["atom"]
    return lrot_list

# gwannier_run:
#   - saves out various parameters into gwannier_params.pkl
#   - runs gwannier.py
#   - runs pygrisb.run.cygutz.driver_cygutz, not as an external process
#   - runs gwannden.py
def gwannier_run(control, wan_hmat, imp, icycle):
        
    os.chdir(control['lowh_directory'])

    # get_locrot_list calculates a rotation matrix based on the xaxis and yaxis
    #   supplied by wannier.inip
    #   wan_hmat['basis'] is populated by read_wan_hmat_basis from wannier.inip, and it
    #       it is used only here and in find_impurity_wan.     
    lrot_list = get_locrot_list(wan_hmat)
    
    # Save all of these variables in gwannier_params.pkl
    params = {}
    params["corbs_list"] = control['impurity_wan']
    params["wpath"] = control['wannier_directory']
    params["lrot_list"] = lrot_list
    params["icycle"] = icycle
    params["ismear"] = control['ismear']
    params["delta"] = control['delta']
    if control['spin_orbit']:
        params["iso"] = 2
    params["lprefix"] = control['allfile']
    params["method"] = control['method']

    # prepare parameter file for mpirun
    with open("gwannier_params.pkl", "wb") as f:
        pickle.dump(params, f, protocol=pickle.HIGHEST_PROTOCOL)

    # run gwannier.py
    cmd = control['mpi_prefix_wannier'] + ' ' + \
            control['comsuitedir'] + "/gwannier.py"
    control['h_log'].write(cmd+"\n")
    hlog_time(control['h_log'], "gwannier start")
    t_start = time.time()
    with open(control['lowh_directory']+'/gwannier.out', 'w') as f:
        subprocess.run(cmd,
                shell=True,
                stdout=f,
                stderr=f,
                check=True,
                )
    hlog_time(control['h_log'], "end", endl="\n")
    t_end = time.time()
    control['t_interface'] += t_end-t_start

    # run pygrisb.run.cygutz.driver_cygutz, not as an external process
    # cygutz calculation
    cmd =  control['comsuitedir'] + "/run_cygutz.py -u 2"
    control['h_log'].write(cmd+"\n")
    hlog_time(control['h_log'], "cygutz start")
    t_start = time.time()
    
    if control['impurity_solution'] == 0:
        if os.path.isfile("GLog.h5"):
            os.remove("GLog.h5")
            
    err = driver_cygutz(path=control['comsuitedir'],
            rmethod="hybr",
            mpi=control['mpi_prefix'].split(),
            tol=1e-6,
            iupdaterho=2,
            )
    
    hlog_time(control['h_log'], "end", endl="\n")
    t_end = time.time()
    control['t_cygutz'] += t_end-t_start
    
    if err > 0.1:
        raise ValueError("cygutz failed to converge.")

    shutil.copy("./Gutz.log", "save_Gutz.log")

    # run gwannden.py
    cmd = control['mpi_prefix_wannier'] + ' ' + \
            control['comsuitedir'] + "/gwannden.py"
    control['h_log'].write(cmd+"\n")
    hlog_time(control['h_log'], "gwannden start")
    t_start = time.time()
    with open(control['lowh_directory']+'/gwannden.out', 'w') as f:
        subprocess.run(cmd,
                shell=True,
                stdout=f,
                stderr=f,
                check=True,
                )
    hlog_time(control['h_log'], "end", endl="\n")
    t_end = time.time()
    control['t_interface'] += t_end-t_start
    os.chdir(control['top_dir'])


def find_allfile(dft_dir):
    files = glob.iglob(dft_dir+"/*.rst")
    for filename in files:
        temp = filename[:-4].split('/')[-1].split('_')
        if temp[0] != 'info':
            return temp[0]
    raise OSError("failed to locate *.rst files in {}!".format(\
            dft_dir))


def dft_risb(control, wan_hmat, imp):
    
    control['h_log'].write("\n\n")
    control['iter_num_outer'] = 0
    
    while control['iter_num_outer'] < control['max_iter_num_outer']:
        
        control['h_log'].write(\
                "************************************************\n"+ \
                "iteration: {}\n".format(control['iter_num_outer'])+ \
                "************************************************\n")
        
        if control['iter_num_outer'] == 0 and \
                control['start_prog'] in ["initcopy"]:
            if not control['restart'] or \
                    len(os.listdir(control['lattice_directory']))==0:
                initial_lattice_directory_setup(control)

        fullrun = control['iter_num_outer'] > 0 or \
                control['start_prog'] in ["initcopy", "wannier"]
        
        # check_wannier_function_input creates comwann.ini        
        check_wannier_function_input(control, wan_hmat)
        
        # wannier_run:
        #   - runs ComWann, but only if the fullrun argument it True
        #   - reads wannier.inip, and uses it to set up control['impurity_wan']
        #       control['impurity_wan'] is used only to populate params["corbs_list"] ,
        #           which is passed to the gutzwiller solver.  
        # wannier_run is implemented only for d,f, and not s,p shells.   
        wannier_run(control, wan_hmat, fullrun)
        
        if control['iter_num_outer'] > 1 or \
                control['start_prog'] in ["initcopy", "wannier", "gwann"]:
                    
            # gwannier_run:
            #   - saves out various parameters into gwannier_params.pkl
            #   - runs gwannier.py
            #   - runs run_cygutz.py
            #   - runs gwannden.py
            gwannier_run(control, wan_hmat, imp, control['iter_num_outer'])
            
        if control['iter_num_outer'] > 1 or \
                control['start_prog'] in ["initcopy", "wannier", \
                "gwann", "dft"]:
            # prepare_dft_input makes sure that  wannier_den_matrix.dat is in 
            #       the right place for rspflapw to use it.
            prepare_dft_input(control)
            
            # run_dft runs rspflapw
            run_dft(control)
            
        if control['end_prog'] in ["dft"]:
            break

        # write_conv gathers info from rspflapw's out file, from GLog.h5,
        #   and from GIter.h5, and writes it out in convergence.log        
        write_conv(control)

        # strict charge density convergence
        # or combination of charge density and total energy
        if control['conv_table'][-1][1] < 1.e-8  \
                or (control['iter_num_outer'] > 1
                and control['conv_table'][-1][1] < control['cc'] and \
                abs(control['conv_table'][-1][2]- \
                control['conv_table'][-2][2]) < control['ec']):
            break
        
        control['iter_num_outer']=control['iter_num_outer']+1


def lqsgw_risb(control, wan_hmat, imp):
    
    control['h_log'].write("\n\n")
    control['iter_num_outer'] = 0
    
    # check_wannier_function_input creates comwann.ini
    check_wannier_function_input(control, wan_hmat)
    
    # wannier_run:
    #   - runs ComWann
    #   - reads wannier.inip, and uses it to set up control['impurity_wan']
    #       control['impurity_wan'] is used only to populate params["corbs_list"] ,
    #           which is passed to the gutzwiller solver.  
    # wannier_run is implemented only for d,f, and not s,p shells.   
    wannier_run(control, wan_hmat)
    
    # gwannier_run:
    #   - saves out various parameters into gwannier_params.pkl
    #   - runs gwannier.py, which calls if_wannier in ifwannier.py
    #   - runs run_cygutz.py, which calls run_cygutz from pygrisb.run.cygutz
    #   - runs gwannden.py, which calls wannier_den_matrix in ifwannier.py
    gwannier_run(control, wan_hmat, imp, control['iter_num_outer'])


# Some notes on the five fortran programs/exe's in  CyGUTZ:
   
# CyGErr, CyGErrR: program gerr
#cygutz/exec_cygerr  - CyGErr
#gerr: Function: calculate vector function error of Gutzwiller nonlinear equations.
#
#CyGInit, CyGInitR:  program ginit
#cygutz/driver_cygutz - CyGInit 
#ginit: 
#! Function: initial processing, including necessary transformation 
#!           of the bare Hamiltonian and append it 
#!           to the "GBareH.h5" files, 
#!           and generating initial solution vectors.
#
# CyGFinal, CyGFinalR: 
#cygutz/driver_cygutz  - CyGFinal
#cygutz/gfun_solver -  CyGFinal
#gfinal:
#! Function: calculate once more quasi-particle (qp) band structure,
#!           total energy of the model, and optionally the renormalized
#!           electron density.
#
# CyGBand, CyGBandR
#cygutz/exec_cygband  - CyGBand
#gbands:
#! Function: calculate quasi-particle (qp) band structure
#!           and qp local quantities.
#
#program exespci calls gspci_gh5exe, an exe named exe_spcir, exe_spci, 
#   seems to be called by gs_driver/driver, which is called by cygutz/solve_hembed_list



def log_time(control):
    control['h_log'].write("\ntimings:\n")
    control['h_log'].write("  rspflapw : {:.1f}s\n".format(\
            control['t_rspflapw']))
    control['h_log'].write("  comwann  : {:.1f}s\n".format(\
            control['t_comwann']))
    control['h_log'].write("  interface: {:.1f}s\n".format(\
            control['t_interface']))
    control['h_log'].write("  cygutz   : {:.1f}s\n".format(\
            control['t_cygutz']))
    control['h_log'].write("  subtotal : {:.1f}s\n".format(\
            control['t_rspflapw']+control['t_comwann']+\
            control['t_interface']+control['t_cygutz']))
    control['h_log'].write("  total    : {:.1f}s\n".format(\
            control['t_total']))


def driver():
    t_start = time.time()
    
    # delete venviron of GWIEN
    if "WIEN_GUTZ_ROOT" in os.environ:
        del os.environ["WIEN_GUTZ_ROOT"]

    # init_comrisb:
    #   - reads in  comrisb.ini
    #   - loads ginit.json into  control['impurity_problem'] and control['impurity_solution']
    #   - parses its arguments: 
    #       -s or --startp sets control['start_prog']
    #       -e or --endp sets control['end_prog']
    #       -c or --continuous sets control['restart']        
    control,wan_hmat,imp = init_comrisb()
    
    initial_file_directory_setup(control)
    
    if control["method"] == "lqsgw+risb":
        lqsgw_risb(control, wan_hmat, imp)
    else:
        dft_risb(control, wan_hmat, imp)
        
    t_end = time.time()
    control['t_total'] = t_end-t_start
    log_time(control)
    control['h_log'].close()


if __name__ == '__main__':
    driver()
