import glob, sys, os, subprocess


def get_dir_list(key="V", chk_link=True):
    if "-dkey" in sys.argv:
        key = sys.argv[sys.argv.index("-dkey")+1]
    dir_list = glob.glob(f'{key}*')
    v_list = []
    for i, v in enumerate(dir_list):
        if "_" not in key and "_" in v:
            continue
        if chk_link and os.path.islink(v):
            continue
        v_list.append(int(v.split(key)[1]))
    v_list.sort()
    dir_list = [f'{key}{v}' for v in v_list]
    print(f"work in dirs {dir_list}.")
    if "-silent" not in sys.argv:
        ans = input("proceed? [Y/n]")
        if "n" in ans.lower():
            sys.exit()
    return dir_list


def batch_job_slurm(jname="dftg",
        dir_template='./template',
        dir_work='case',
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
        os.chdir(f"{dname}/{dir_work}")
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
