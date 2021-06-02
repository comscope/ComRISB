import numpy as np
import sys, os, subprocess, shutil, h5py, warnings
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
                locrot.append(map(float, [line[20:30], line[30:40], \
                        line[40:50]]))
                readline += 1
                if readline == 3:
                    # wien2k convention
                    locrot = np.asarray(locrot).T
                    for i in range(mult):
                        # maybe there are problems here.
                        locrot_list.append(locrot)
    return locrot_list


def get_scf_data(file_scf='./case.out', name="volume", startpos=15):
    '''Get data with name from scf file file_scf.
    '''
    with open(file_scf, 'r') as f:
        res = None
        for line in f:
            if name in line:
                line = line[startpos:-2]
                res = float(line)
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

    # load init script
    cmd0 = [s.strip() for line in open("init.inp", "r").readlines() for s in
            line[:-1].split(";") if s != ""]

    # loop over to create subdirectory list
    with open('vol_record.txt', 'w') as f:
        for i, vfrac in enumerate(vfrac_list):
            vdir = f"V{i}/dft"
            f.write(f'V{i} {vfrac:6.3f}\n')
            os.makedirs(vdir, exist_ok=True)
            os.chdir(vdir)
            print(f"working at {vdir}")
            cmd = cmd0 + [f"{vfrac:.2f}"]
            print(f"cmd: {cmd}")
            subprocess.run(cmd)
            os.chdir("../..")


def h5wrt_data(f, path, data):
    if path in f:
        del f[path]
    f[path] = data


def h2get_data_volume_list(path='dft',
        d_name="ETOT ",
        startpos=13,
        d_scale=1.,
        list_name="e_list",
        dlabel="E (eV/primitive cell)",
        fname="ev",
        ):
    '''Loop over all the directories to get the data vs volume data,
    and save it to the metadata file results.h5 at path.
    '''
    v_list = []
    e_list = []
    cwd = os.getcwd()+'/'
    for dname in get_dir_list(chk_link=False):
        try:
            os.chdir(f"{dname}/{path}")
            v = get_scf_data()
            e = get_scf_data(name=d_name, startpos=startpos)
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
                '  Vlist -- generate directories for a range of volumes;\n'
                '           must provide vrange.inp file. \n'
                '  batch_init_ga -- init ga for all the directories;\n'
                '  batch_init_mott -- init ga-mott for all directories;\n'
                '  batch_setuj "u1 u2 .." "j1 j2 ..." -- set uj for dirs;\n'
                '  batch_modify_ga ... -- modify ga for all directories;\n'
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

    if 'Vlist' == sys.argv[1]:
        create_dir_vlist(vfrac_min=vfrac_min, vfrac_max=vfrac_max,
                vfrac_step=vfrac_step, dir_template=dir_template)
    elif 'batch_modify_ga' == sys.argv[1]:
        args = sys.argv[2:]
        gwien.batch_modify_ga_setup(args)
    elif '-jobq' == sys.argv[1]:
        batch_job_slurm(dir_work="dft")
    elif '-ev' == sys.argv[1]:
        h2get_data_volume_list(path=sys.argv[2],
                d_name="ETOT ",
                d_scale=units.rydberg_to_ev,
                list_name="e_list",
                dlabel="E (eV/primitive cell)",
                fname="ev")
    else:
        raise ValueError('Inline option not defined!')



if __name__=='__main__':
    steps(vfrac_step=0.05)
