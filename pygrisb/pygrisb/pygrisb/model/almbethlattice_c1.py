# Author: Yongxin Yao

'''Anderson lattice model with conduction site in beth lattice.
'''

import numpy as np
import pygrisb.gutz.gparam as gp
import pygrisb.gutz.atoms as gatoms
import pygrisb.model.tbASE as tb
import pygrisb.model.special as special
import pygrisb.mbody.coulomb_matrix as cm
import pygrisb.symm.angular_momentum as am


def gutz_model_setup(
        u=0.0,
        v=0.4,
        eta_f=0.,
        eps_c=0.,
        spindeg=True,
        iembeddiag=-3,
        delta=0.01,
        nmesh=5000,
        dtype=np.complex,
        ):
    '''Set up Gutzwiller calculations for periodic anderson mode with
    semi-circular DOS for conduction band.

    Parameters:

    * u: real number
      Hubbard U.
    * nmesh: interger number
      number of energy mesh

    Result:

    Create all the necessary input file of ``GParam.h5``, ``GBareH.h5``
    for *CyGutz* calculation.
    '''
    # get semi-circular class, predifined in pygrisb.model.special
    sc = special.semicircular()

    # get semi-circular energy mesh
    e_list = sc.get_e_list_of_uniform_wt(nmesh=nmesh)
    # shift conduction band center
    e_list += eps_c

    # get atomic structure with 1 atom per unit cell.
    a = tb.AtomsTB("H", [(0, 0, 0)], cell=(1, 1, 1))

    # specify single s-orbital and orbital dimensions.
    a.set_orbitals_spindeg(orbitals=[('f','s')], spindeg=spindeg)

    # set maximal number of bands (here we have two bands.)
    num_band_max = 2

    eps_f = (eta_f-1)*u/2.
    hk_list = [np.array([[eps_f, v],[v,e+0.j]]) for e in e_list]

    # get the tight-binding model
    aTB = tb.TB(a, hk_list=hk_list)

    # here the energy mesh is the same as k-points.
    kps = e_list
    num_k = len(kps)

    # set uniform k-point weight.
    kps_wt = 1.0 / num_k * np.ones((num_k))

    # Include the spin degeneracy to k-point weight.
    if aTB.Atoms.spindeg:
        kps_wt *= 2

    # We will work with grand canonical model, num_e will not be used.
    num_e = 0.0

    # list of one-body part of the local Hamiltonians.
    # here the s-orbital is located at zero.
    # always in complex format
    h1e_list = [np.array([[eps_f]], dtype=np.complex)]

    # create ``GBareH.h5`` file
    # ensemble is set to 1 for grand canonical system.
    gbndinfo = gp.gparam_bnds(kps_wt,
            num_e,
            num_band_max,
            ensemble=1,
            h1e_list=h1e_list,
            delta=delta)
    gbndinfo.h5save()
    aTB.save_bareham(kps)

    # prepare ``GParam.h5`` file.
    # dummy atom and lattice.
    gatms = gatoms.cor_structure(symbols=["H"], scaled_positions=[[0,0,0]],
            cell=np.identity(3))
    gatms.update({
            "dc_mode": 0,
            "imap_list": [0],
            "na2_list": [2],
            "iso": 1,
            "ispin": 1,
            # one spin channel only for now
            "db2sab_list": [np.eye(1, dtype=dtype)],
            "symbol_matrix_list": [np.arange(1,2).reshape(1,1)],
            # no dc
            "dc_u_avg_list": [0.],
            "dc_j_avg_list": [0.],
            "giembeddiag": iembeddiag,
            "nval_bot_list": [0],
            "nval_top_list": [2],
            })
    gatms.set_full_symm_adaptive_basis()
    # local single-particle matrix basis.
    gatms.set_matrix_basis_list(realhemb=(dtype==np.float))
    # set Coulomb matrix, one spin channel for now
    v2e_1 = np.zeros((1,1,1,1), dtype=dtype)
    v2e_1[0,0,0,0] = u
    v2e_2 = cm.add_spin_comp_to_coul_mat(v2e_1)
    v2e_2 = cm.apply_a4_unitary_trans(v2e_2, gatms.ps["db2sab_list"][0])
    v2e_2 = np.asarray(v2e_2, dtype=dtype)
    gatms.update({"v2e_list": [v2e_2]})
    # spin operator
    samop = am.samop_csh2(0.0)
    s_vec = samop.am_op
    gatms.set_s_vec_csh2_list([s_vec])
    gatms.set_sl_vector_sab_list(mode="sz_only")
    gpa = gp.gparam_atms(gatms.ps["na2_list"], gatms.ps)
    gpa.h5save()



if __name__=="__main__":
    gutz_model_setup(u=2.)
