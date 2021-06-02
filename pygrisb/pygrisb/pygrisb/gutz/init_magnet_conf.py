#!/usr/bin/env python
import h5py, numpy, os, inquirer
from pygrisb.basic.units import rydberg_to_ev
import pygrisb.mbody.double_counting as dc


def get_b_field_list():
    '''get the list of magnetic field by q&a.
    '''
    question = [
            inquirer.List('unit',
                    message="choose unit used in CyGutz calculation",
                    choices=[("Rydberg (13.6)", rydberg_to_ev),
                            ("eV (1)", 1.)],),
            inquirer.List('ivext',
                    message="which way to apply vext",
                    choices=[("initial step only", 0), \
                            ("all iterations", 1)],),
            ]

    ans = inquirer.prompt(question)
    unit = ans["unit"]
    givext = ans["ivext"]

    with h5py.File('GParam.h5', 'r') as f:
        ispin = f["/"].attrs["ispin"]
        if ispin == 1:
            print("ispin = 1 in GParam.h5: no magnetic calculation.")
            return None, None
        num_imp = f["/"].attrs["num_imp"]
        imap_list = f['/'].attrs['imap_list']

    print(' total {} impurities with equivalence indices \n {}'.format( \
            num_imp, imap_list))
    bvec_list = []
    for imp in range(num_imp):
        if imap_list[imp] == imp:
            print('\n impurity {}'.format(imp))
            # field direction
            def direction_accepted(_, x):
                try:
                    x = list(map(float, x.split()))
                    return len(x) == 3
                except:
                    return False
            def float_accepted(_, x):
                try:
                    x = float(x)
                    return x >= 0
                except:
                    return False

            question = [
                    inquirer.Text('b_direction',
                            message="enter field direction x y z " + \
                                    "seperated by space. (e.g., 0 0 1)",
                            validate=direction_accepted),
                    inquirer.Text('b_magnitude',
                            message="enter b field magnitude" + \
                                    " (eV/Bohr magneton)",
                            validate=float_accepted),
                    ]
            ans = inquirer.prompt(question)
            vec = numpy.array(list(map(float, ans["b_direction"].split())))
            bm = float(ans["b_magnitude"])/unit
            bvec_list.append(bm*vec)
        else:
            bvec_list.append(bvec_list[imap_list[imp]])
    return bvec_list, givext


def get_sym_2darray(a, imp=0):
    '''get symmetrized 2d array with error reporting.
    '''
    with h5py.File('GParam.h5', 'r') as f:
        path = f'/impurity_{imp}/matrix_basis'
        h_list = f[path][()]
        a_sym = numpy.zeros_like(a)
        for h in h_list:
            a_sym += h.T.conj().dot(a).trace()*h
        sym_err = numpy.max(numpy.abs(a-a_sym))
    return a_sym, sym_err


def get_vext_list(bvec_list):
    '''get the list of magnetic potential matrix in the single particle space,
    given the list of local magnetic field.
    '''
    vext_list = []
    max_sym_err = 0.
    with h5py.File('GParam.h5', 'r') as f:
        iso = f["/"].attrs["iso"]
        for imp, bvec in enumerate(bvec_list):
            prepath = f"/impurity_{imp}"
            sx = f[prepath+'/sx'][()]
            sy = f[prepath+'/sy'][()]
            sz = f[prepath+'/sz'][()]
            vext = -(bvec[0]*sx + bvec[1]*sy + bvec[2]*sz)*2
            if iso == 2: # spin-orbit coupling is present
                lx = f[prepath+'/lx'][()]
                ly = f[prepath+'/ly'][()]
                lz = f[prepath+'/lz'][()]
                vext -= bvec[0]*lx + bvec[1]*ly + bvec[2]*lz
            vext_sym, sym_err = get_sym_2darray(vext, imp)
            max_sym_err = max(max_sym_err, sym_err)
            vext_list.append(vext_sym)
    if max_sym_err > 1.e-5:
        print(" Warning:", end='')
    print(f' maximal symmetrization error of vext = {max_sym_err:.2e}')
    return vext_list


def get_vext_given_1pdm_list(dm_list):
    '''get the external potential for lda+u calculation
    given the initial one-particle density matrix.
    '''
    with h5py.File("GParam.h5", "r") as f:
        javg_list = f["/dc_j_avg"][()]
        uavg_list = f["/dc_u_avg"][()]
        v2e_list = []
        for i in range(f["/"].attrs["num_imp"]):
            v2e_list.append(f[f"/impurity_{i}/V2E"][()].T)

    vext_list = [-dc.get_vdc_hf(v2e, dm)+ \
            dc.get_vdc_fll(uavg, javg, numpy.trace(dm))*\
            numpy.eye(dm.shape[0]) \
            for v2e, dm, uavg, javg in zip(v2e_list, dm_list, \
            uavg_list, javg_list)]
    return vext_list


def chk_local_one_body(vext_list):
    if not os.path.isfile("GPARAMBANDS.h5"):
        return
    db2sab_list = []
    sx_list = []
    sy_list = []
    sz_list = []
    lx_list = []
    ly_list = []
    lz_list = []
    with h5py.File('GParam.h5', 'r') as f:
        num_imp = f["/"].attrs["num_imp"]
        for imp in range(num_imp):
            prepath = f"/impurity_{imp}"
            db2sab_list.append(f[prepath+"/db_to_sab"][()])
            sx_list.append(f[prepath+"/sx"][()])
            sy_list.append(f[prepath+"/sy"][()])
            sz_list.append(f[prepath+"/sz"][()])
            lx_list.append(f[prepath+"/lx"][()])
            ly_list.append(f[prepath+"/ly"][()])
            lz_list.append(f[prepath+"/lz"][()])
    # start actually analysis
    numpy.set_printoptions(precision=2, suppress=True)
    for imp in range(num_imp):
        h = vext_list[imp]
        w, v = numpy.linalg.eigh(h)
        dm = v[:,0:1].dot(v[:,0:1].T.conj())
        dm_sym, _ = get_sym_2darray(dm, imp)
        print(" <Sx> = {}".format(numpy.sum(dm_sym*sx_list[imp])))
        print(" <Sy> = {}".format(numpy.sum(dm_sym*sy_list[imp])))
        print(" <Sz> = {}".format(numpy.sum(dm_sym*sz_list[imp])))
        print(" <Lx> = {}".format(numpy.sum(dm_sym*lx_list[imp])))
        print(" <Ly> = {}".format(numpy.sum(dm_sym*ly_list[imp])))
        print(" <Lz> = {}".format(numpy.sum(dm_sym*lz_list[imp])))


def h5wrt_gmagnet(vext_list, g_ivext):
    with h5py.File("GParam.h5", 'a') as f:
        if "/vext" in f:
            del f["/vext"]
        for imp, vext in enumerate(vext_list):
            f[f'/vext/impurity_{imp}/v'] = vext
        f['/vext'].attrs['givext'] = g_ivext


def init_magnet_conf():
    '''
    initialize the the magnetic configuration for magnetic calculation.
    '''
    bvec_list, g_ivext = get_b_field_list()
    if bvec_list is None:
        return
    vext_list = get_vext_list(bvec_list)
    chk_local_one_body(vext_list)
    h5wrt_gmagnet(vext_list, g_ivext)


def init_magnet_conf_with_init_dm(dm_list):
    '''
    initialize the the magnetic configuration for magnetic calculation
    based on given initial one-particle density-matrix.
    '''
    vext_list = get_vext_given_1pdm_list(dm_list)
    chk_local_one_body(vext_list)
    h5wrt_gmagnet(vext_list, g_ivext=0)


if __name__ == "__main__":
    init_magnet_conf()
