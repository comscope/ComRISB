import numpy as np
import sys, os, glob, subprocess, shutil, fileinput, h5py, warnings
import pygrisb.basic.splot as splot
import pygrisb.basic.units as units
from pygrisb.dft.util import get_dir_list, batch_job_slurm
from pygrisb.run import gwien


'''help routines for processing wien2k data. Assume the current work directory
has a subfolder template with reference case.struct and possibly case.inso.
'''


def get_rotations(case_fname):
    '''get the rotation list in struct file.
    '''
    with open(str(case_fname), 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if 'SYMMETRY OPERATIONS' in line:
                break
        nops = int(lines[i].split()[0])
        lines = lines[i+1:]
        rot_list = []
        for i in range(0,nops*4,4):
            rot = []
            lok = True
            for j in range(3):
                line = lines[i+j]
                rot.append([int(line[:2]), int(line[2:4]), int(line[4:6])])
                if abs(float(line[6:17])) > 1.e-6:
                    lok = False
                    break
            if lok:
                if abs(np.linalg.det(rot)-1)<1.e-6:
                    rot_list.append(rot)
    return rot_list


def get_equivalent_atom_indices(case_fname):
    '''get equivalent atomic indices according to the atom labels in
    case.struct file.
    '''
    with open(str(case_fname), 'r') as f:
        idx_equivalent_atoms = []
        ibase = 0
        for line in f:
            if 'MULT=' in line:
                mult = int(line[line.index("MULT=")+5:].split()[0])
                for i in range(mult):
                    idx_equivalent_atoms.append(ibase)
                ibase += mult
    return idx_equivalent_atoms


def get_local_rotations(case_fname):
    '''get list of locrot in struct file.
    '''
    with open(str(case_fname), 'r') as f:
        locrot_list = []
        readline = 100
        for line in f:
            if 'MULT=' in line:
                mult = int(line.split()[1])
            elif 'LOCAL ROT MATRIX' in line:
                readline = 0
                locrot = []
            if readline < 3:
                locrot.append(list(map(float, [line[20:30], line[30:40], \
                        line[40:50]])))
                readline += 1
                if readline == 3:
                    # wien2k convention
                    locrot = np.asarray(locrot).T.tolist()
                    for i in range(mult):
                        # maybe there are problems here.
                        locrot_list.append(locrot)
    return locrot_list


def get_scf_data(file_scf='./case.scf', name=":VOL"):
    '''Get data with name from scf file f_scf.
    '''
    with open(file_scf, 'r') as f:
        res = None
        for line in f:
            if name in line:
                line = line.split('=')
                res = float(line[1])
        if res is None:
            raise ValueError(f"{name} not found!")
        return res


def get_ca_data(fname='./case.struct'):
    '''get c/a data from struct file.
    '''
    with open(fname, 'r') as f:
        line = f.readlines()[3]
        a_list = []
        for j in range(3):
            a_list.append(float(line[j*10:(j+1)*10]))
        if abs(a_list[0]-a_list[1]) > 1.e-2:
            warnings.warn("a,b not even close while checking c/a!")
        return a_list[2]/a_list[0]


def create_struct(dir_template, abc_scale, case='case'):
    '''Create a case.struct file with lattice constants scaled by abc_scale
    of the reference struct file in dir_template.
    '''
    f_struct = glob.glob(dir_template+'/*struct')
    if len(f_struct) == 0:
        raise IOError('No struct file existing in {}!'.format(dir_template))
    f_struct = f_struct[0]
    abc_scale = np.asarray(abc_scale)
    with open(f_struct, 'r') as fs:
        with open(case+'.struct', 'w') as fd:
            for i, line in enumerate(fs.readlines()):
                if i == 3:
                    a_list = []
                    for j in range(3):
                        a_list.append(float(line[j*10:(j+1)*10]))
                    a_list = np.array(a_list)
                    a_list *= abc_scale
                    line = '{:10.6f}{:10.6f}{:10.6f}'.format(*a_list) + \
                            line[30:]
                fd.write(line)


def update_case_ref(vfrac_min):
    '''Update case.struct based on the initial run of the Vmin sample.
    '''
    cwd = os.getcwd()+'/'
    os.chdir('./template')
    a_scale = 1./vfrac_min**(1./3)
    create_struct('../Ref_Vmin/case', [a_scale, a_scale, a_scale])
    os.chdir(cwd)


def create_dir1(abc_scale, dname, dir_template='./template'):
    '''Create a subdirectory dname in the current directory with
    volumes changing by a vfrac wth repect to the reference structure
    in the dir_template.
    '''
    # create it if necessary
    if not os.path.isdir(dname):
        os.mkdir(dname)
    if not os.path.isdir(dname+'/case'):
        os.mkdir(dname+'/case')

    cwd = os.getcwd()+'/'
    # go to the subdirectory
    os.chdir(dname+'/case')

    # cretae the case.struct file
    create_struct(cwd+dir_template, abc_scale)

    # return to the upper work directory
    os.chdir(cwd)


def create_dir_vlist(vfrac_min=0.7, vfrac_max=1.3, vfrac_step=0.05, \
        dir_template='./template'):
    '''Create a list of subdirectories in the current directory with
    volumes changing by a fraction of vfrac_step in the range of
    [vfrac_min, vfrac_max] wth repect to the reference structure
    in the dir_template.
    '''

    if not os.path.isdir(dir_template):
        raise IOError('{} does not exist!'.format(dir_template))
    # Get volume fraction list
    vfrac_list = np.arange(vfrac_min, vfrac_max+1.e-5, vfrac_step)

    # loop over to create subdirectory list
    with open('vol_record.txt', 'w') as f:
        for i, vfrac in enumerate(vfrac_list):
            # subdirectory name
            dname = 'V{}'.format(i)
            f.write(f'{dname} {vfrac:6.3f}\n')
            a_scale = vfrac**(1./3)
            create_dir1([a_scale, a_scale, a_scale], dname,
                    dir_template=dir_template)


def create_dir_calist(ca_min=0.8, ca_max=1.2, ca_step=0.05, \
        dir_template='./'):
    '''Create a list of subdirectories in the current directory with
    c/a ratio changing by a ca_step in the range of
    [ca_min, ca_max] wth repect to the reference structure
    in the dir_template.
    '''
    if not os.path.isdir(f"{dir_template}/case"):
        raise IOError('{} does not exist!'.format(dir_template))
    # Get volume fraction list
    ca_list = np.arange(ca_min, ca_max+1.e-5, ca_step)

    prefix = dir_template.split("/")[-1]
    # loop over to create subdirectory list
    with open('ca_record.txt', 'w') as f:
        for i, ca in enumerate(ca_list):
            # subdirectory name
            dname = f'{prefix}_{i}'
            f.write(f'{dname} {ca:6.3f}\n')
            if np.abs(ca-1) < 1.e-6:
                os.symlink(dir_template, dname)
            else:
                abc_scale = [ca**(-1/3), ca**(-1/3), ca**(2/3)]
                create_dir1(abc_scale, dname,
                        dir_template=f"{dir_template}/case")


def batch_init_lapw():
    '''Loop over all the directories to run init_lapw.
    '''
    cmd = [os.environ['WIENROOT']+'/init_lapw']
    fin = open("init_lapw.inp", "r")
    line = fin.readlines()[0]
    if "-" in line:
        # batch mode
        args = line.split()
        cmd += ['-b'] + args
        print(f"init_lapw -b parameters: {args}")
    else:
        os.environ["EDITOR"] = ""
        print("init_lapw with input from init_lapw.inp.")

    f = open('binit_lapw.log', 'w')
    cwd = os.getcwd()+'/'

    for dname in get_dir_list():
        os.chdir(dname+'/case')
        print(f"initilising {dname}")
        # back to the very begining.
        fin.seek(0)
        res = subprocess.run(cmd, stdin=fin, stdout=f)
        assert(res.returncode == 0), "init_lapw error!"
        os.chdir(cwd)
    f.close()
    fin.close()


def modify_emax_case_in1(emax=7.5):
    '''Modify the emax value in case_in1 for the case with
    spin-orbit calculation.
    '''
    files = glob.glob("case.in1*")
    for case_in1 in ["case.in1", "case.in1c"]:
        if case_in1 in files:
            for line in fileinput.input(files=case_in1,
                    inplace=True):
                if 'emax' in line:
                    line = line.replace(line[33:38],
                            '{:5.1f}'.format(emax))
                print(line, end='')


def batch_initso_lapw(dir_template='./template', emax=7.5):
    '''Loop over all the directories to run initso_lapw -- actually,
    because there is no batch mode provided, it simply copy the case.inso
    from dir_template and modify the EMAX value to 7.5 in case.in1.
    '''
    cwd = os.getcwd()+'/'
    for dname in get_dir_list():
        os.chdir(dname+'/case')
        shutil.copy(cwd+'/'+dir_template+'/case.inso', './')
        shutil.copy(cwd+'/'+dir_template+'/case.in2c', './')
        modify_emax_case_in1(emax=emax)
        os.chdir(cwd)


def h5wrt_data(f, path, data):
    if path in f:
        del f[path]
    f[path] = data


def h2get_data_volume_list(path='lapw', d_name=":ENE", d_scale=1.,
        list_name="e_list", dlabel="E (eV/primitive cell)",
        fname="ev"):
    '''Loop over all the directories to get the data vs volume data,
    and save it to the metadata file results.h5 at path.
    '''
    v_list = []
    e_list = []
    cwd = os.getcwd()+'/'
    for dname in get_dir_list(chk_link=False):
        try:
            os.chdir(dname+'/case/'+path)
            v = get_scf_data(name=":VOL")
            e = get_scf_data(name=d_name)
            v_list.append(v)
            e_list.append(e)
        except FileNotFoundError:
            warnings.warn(f"dir {dname} is missing!")
            pass
        os.chdir(cwd)

    # Ryd/Bohr units to eV/A units
    v_list = np.array(v_list)*units.bohr_to_angstrom**3
    e_list = np.array(e_list)*d_scale

    with h5py.File('results.h5', 'a') as f:
        h5wrt_data(f, f"/{path}/v_list", v_list)
        h5wrt_data(f, f"/{path}/{list_name}", e_list)

    splot.xy_plot(v_list, e_list, xlabel='V ($\AA^{3}$/primitive cell)',
            ylabel=dlabel, fsave=f'{fname}_{path}.pdf')


def h2get_data_ca_list(path='lapw', d_name=":ENE", d_scale=1.,
        list_name="e_list", dlabel="E (eV/primitive cell)",
        fname="eca"):
    '''Loop over all the directories to get the data vs c/a,
    and save it to the metadata file results.h5 at path.
    '''
    ca_list = []
    e_list = []
    cwd = os.getcwd()+'/'
    for dname in get_dir_list(chk_link=False):
        try:
            os.chdir(dname+'/case/'+path)
            ca = get_ca_data()
            e = get_scf_data(name=d_name)
            ca_list.append(ca)
            e_list.append(e)
        except FileNotFoundError:
            warnings.warn(f"dir {dname} is missing!")
            pass
        os.chdir(cwd)

    e_list = np.array(e_list)*d_scale
    if "_" in dname:
        prefix = f"{dname.split('_')[0]}"
    else:
        prefix = ""
    with h5py.File('results.h5', 'a') as f:
        h5wrt_data(f, f"/{prefix}/{path}/ca_list", ca_list)
        h5wrt_data(f, f"/{prefix}/{path}/{list_name}", e_list)

    splot.xy_plot(ca_list, e_list, xlabel='c/a',
            ylabel=dlabel, fsave=f'{prefix}{fname}_{path}.pdf')


def compare_data_volume_plot(path_list, fres='results.h5',
        list_name="e_list",
        d_label="E (eV/primitive cell)",
        fname="ev"):
    '''Compare lapw/lapwso/lapwsog in a plot.
    '''
    v_list = []
    e_list = []
    pattern_list = []
    label_list = []
    with h5py.File(fres, 'r') as f:
        postfix=''
        for path in path_list:
            if path in f:
                label_list.append(path)
                e_list.append(f[f"/{path}/{list_name}"][()])
                v_list.append(f[f"/{path}/v_list"][()])
                pattern_list.append('-o')
                postfix += path.split('/')[0]
    if v_list != []:
        splot.xy2_plot(v_list, e_list, pattern_list, label_list,
                xlabel='V ($\AA^{3}$/primitive cell)',
                ylabel=d_label,
                fsave=f'{fname}_{postfix}.pdf')


def batch_save_lapw(sdir='lapw', args=['-f']):
    '''Loop over all the directories to save_lapw.
    '''
    if '-d' in sys.argv:
        sdir = sys.argv[sys.argv.index('-d')+1]
    cmd = [os.environ['WIENROOT']+'/save_lapw', '-d'] + [sdir] + args
    cwd = os.getcwd()+'/'
    for dname in get_dir_list():
        os.chdir(dname+'/case')
        subprocess.run(cmd)
        os.chdir(cwd)


def batch_struct2cif():
    '''Loop over all the directories to convert struct file to cif.
    '''
    cmd = [os.environ['WIENROOT']+'/x', 'struct2cif', '-f', 'case']
    cwd = os.getcwd()+'/'
    for dname in get_dir_list():
        os.chdir(dname+'/case')
        subprocess.run(cmd)
        os.chdir(cwd)


def run_lapw(args=['-i', '70'], nproc=1):
    '''Loop over all the directories to run_lapw using nproc processors
    with arguments provided in args.
    '''
    if '-sp' in sys.argv:
        cmd = [os.environ['WIENROOT']+'/runsp_lapw'] + args
    if '-fsm' in sys.argv:
        mm = sys.argv[sys.argv.index('-fsm')+1]
        cmd = [os.environ['WIENROOT']+'/runfsm_lapw', "-m", mm] + args
    else:
        cmd = [os.environ['WIENROOT']+'/run_lapw'] + args
    cwd = os.getcwd()+'/'
    if '-p' in sys.argv:
        nproc = int(sys.argv[sys.argv.index('-p')+1])
    for i,dname in enumerate(get_dir_list()):
        os.chdir(dname+'/case')
        proc = subprocess.Popen(cmd)
        os.chdir(cwd)
        if (i+1) % nproc == 0:
            proc.communicate()


def batch_copy_todir(src, dst):
    '''Loop over all the directories and copy subfolder src to dst.
    '''
    cwd = os.getcwd()+'/'
    for dname in get_dir_list():
        os.chdir(dname+'/case')
        if not os.path.lexists(dst):
            os.mkdir(dst)
        if os.path.isfile(src):
            shutil.copy(src, dst)
        else:
            for f in os.listdir(src):
                shutil.copy(src+'/'+f, dst)
        os.chdir(cwd)


def batch_rename_dir(src, dst):
    '''Loop over all the directories and rename subfolder.
    '''
    cwd = os.getcwd()+'/'
    for dname in get_dir_list():
        os.chdir(dname+'/case')
        if os.path.exists(src):
            os.rename(src, dst)
        os.chdir(cwd)


def steps(vfrac_min=0.70, vfrac_max=1.30, vfrac_step=0.05, \
        dir_template='./template'):
    if len(sys.argv) == 1 or '-h' in sys.argv[1]:
        print(('Please provide with inline argument chosen from below:\n'
                '  Vmin -- setup init. calculation of the min. vol. point;\n'
                '  update_case_ref -- update esp. RMT values of ref. case;\n'
                '  Vlist -- generate directories for a range of volumes;\n'
                '           must provide vrange.inp file. \n'
                '  calist -d V* -- generate directories for a range of c/a;\n'
                '           must provide carange.inp file. \n'
                '  batch_struct2cif -- covert all struct file to cif. \n'
                '  batch_init_lapw -- init_lapw all the directories;\n'
                '           must provide init_lapw.inp file. \n'
                '  batch_initso_lapw -- initso_lapw all the directories;\n'
                '  batch_init_ga -- init ga for all the directories;\n'
                '  batch_init_mott -- init ga-mott for all directories;\n'
                '  batch_setuj "u1 u2 .." "j1 j2 ..." -- set uj for dirs;\n'
                '  batch_modify_ga ... -- modify ga for all directories;\n'
                '  batch_run_lapw [-sp, -fsm m] -- run_lapw all directories;\n'
                '  batch_run_lapwso [-sp ]-- run_lapw -so all directories;\n'
                '  batch_run_ga -- run_ga all the directories; \n'
                '  batch_save [-d lapw] -- save_lapw all the directories; \n'
                '  batch_gsave [-d uxjx]  -- save DFT+G results to dirs;\n'
                '  -cp src dst -- copy subfolder src/or file to dst; \n'
                '  -rename-dir src dst -- rename subfilder src to dst; \n'
                '  -jobq -n name -w dir -j jobfile -- batch submit jobs; \n'
                '  -ev dirname -- save e-v data for dirname calc.\n'
                '  -eca dirname -- save e-c/a data for dirname calc.\n'
                '  -mv dirname -- save magmom-v data for dirname calc.\n'
                '  -mca dirname -- save magmom-c/a data for dirname calc.\n'
                '  -pv dirname -- Murnaghan EOS fir for dirname results.\n'
                '  -ev-cmp path1 path2 ... -- Compare e-v data.\n'
                '  -pv-cmp path1 path2 ... -- Compare p-v data.\n'
                '  -mv-cmp path1 path2 ... -- Compare m-v data.\n'
                '  You may append "-p nproc" to specify # of procs in use.\n'
                '  -dkey V5_ to specify the directory list.'))
        sys.exit('Please choose proper inline argument!')

    with open("vrange.inp", "r") as f:
        line = f.readline()
        line = line.split()
    vfrac_min, vfrac_max, vfrac_step = map(float, line)

    if 'Vmin' in sys.argv[1]:
        a_scale = vfrac_min**(1./3)
        create_dir1([a_scale, a_scale, a_scale], 'Ref_Vmin',
                dir_template=dir_template)
        print('Please goto Ref_Vmin/case and finish manual test.')
    elif 'update_case_ref' == sys.argv[1]:
        update_case_ref(vfrac_min)
    elif 'Vlist' == sys.argv[1]:
        create_dir_vlist(vfrac_min=vfrac_min, vfrac_max=vfrac_max,
                vfrac_step=vfrac_step, dir_template=dir_template)
    elif 'calist' in sys.argv[1]:
        with open("carange.inp", "r") as f:
            line = f.readline()
            line = line.split()
        ca_min, ca_max, ca_step = map(float, line)
        dir_template = sys.argv[sys.argv.index("-d") + 1]
        create_dir_calist(ca_min=ca_min, ca_max=ca_max, ca_step=ca_step, \
                dir_template=dir_template)
    elif 'batch_struct2cif' == sys.argv[1]:
        batch_struct2cif()
    elif 'batch_init_lapw' == sys.argv[1]:
        batch_init_lapw()
    elif 'batch_run_lapw' == sys.argv[1]:
        run_lapw(args=['-cc', '0.0001'])
    elif 'batch_run_lapwso' == sys.argv[1]:
        run_lapw(args=['-so', '-i', '70', '-cc', '0.0001'])
    elif 'batch_save' in sys.argv[1]:
        batch_save_lapw()
    elif 'batch_gsave' in sys.argv[1]:
        gwien.batch_gsave()
    elif '-cp' == sys.argv[1]:
        batch_copy_todir(sys.argv[2], sys.argv[3])
    elif '-rename-dir' == sys.argv[1]:
        batch_rename_dir(sys.argv[2], sys.argv[3])
    elif '-ev' == sys.argv[1]:
        h2get_data_volume_list(path=sys.argv[2],
                d_name=":ENE",
                d_scale=units.rydberg_to_ev,
                list_name="e_list",
                dlabel="E (eV/primitive cell)",
                fname="ev")
    elif '-eca' == sys.argv[1]:
        h2get_data_ca_list(path=sys.argv[2],
                d_name=":ENE",
                d_scale=units.rydberg_to_ev,
                list_name="e_list",
                dlabel="E (eV/primitive cell)",
                fname="eca")
    elif '-mca' == sys.argv[1]:
        h2get_data_ca_list(path=sys.argv[2],
                d_name=":MMTOT",
                d_scale=1.,
                list_name="m_list",
                dlabel="$M_{tot}$",
                fname="mca")
    elif '-mv' == sys.argv[1]:
        h2get_data_volume_list(path=sys.argv[2],
                d_name=":MMTOT",
                d_scale=1.,
                list_name="m_list",
                dlabel="$M_{tot}$",
                fname="mv")
    elif 'batch_initso_lapw' == sys.argv[1]:
        batch_initso_lapw()
    elif '-pv' == sys.argv[1]:
        from pygrisb.dft.eos import h5get_mfit_ev
        h5get_mfit_ev(path=sys.argv[2])
    elif '-ev-cmp' == sys.argv[1]:
        compare_data_volume_plot(sys.argv[2:],
                list_name="e_list",
                d_label="E (eV/primitive cell)",
                fname="ev")
    elif '-pv-cmp' == sys.argv[1]:
        compare_data_volume_plot(sys.argv[2:],
                list_name="p_list",
                d_label="P (GPa)",
                fname="pv")
    elif '-mv-cmp' == sys.argv[1]:
        compare_data_volume_plot(sys.argv[2:],
                list_name="m_list",
                d_label="$M_{tot}$",
                fname="mv")
    elif 'batch_init_ga' == sys.argv[1]:
        gwien.batch_init_ga()
    elif 'batch_init_mott' == sys.argv[1]:
        gwien.batch_init_mott()
    elif 'batch_setuj' in sys.argv[1]:
        args = ['-unique_u_ev', sys.argv[2], '-unique_j_ev', sys.argv[3]]
        gwien.batch_modify_ga_setup(args)
    elif 'batch_modify_ga' == sys.argv[1]:
        args = sys.argv[2:]
        gwien.batch_modify_ga_setup(args)
    elif '-jobq' == sys.argv[1]:
        batch_job_slurm()
    elif 'batch_run_ga' == sys.argv[1]:
        gwien.run_ga()
    else:
        raise ValueError('Inline option not defined!')



if __name__=='__main__':
    steps(vfrac_step=0.05)
