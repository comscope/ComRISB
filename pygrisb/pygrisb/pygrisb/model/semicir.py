# Author: Yongxin Yao

'''One band with semi-circular dos.
'''

import numpy as np
import h5py
import pygrisb.gutz.gparam as gp
import pygrisb.gutz.atoms as gatoms
import pygrisb.model.tbASE as tb
import pygrisb.model.special as special
import pygrisb.symm.angular_momentum as am


def gutz_model_setup(
        u=0.0,
        nmesh=5000,
        mu=0,
        iembeddiag=-2,
        num_e=0.,  # electron filling, not referenced for grand canonical model
        ensemble=1,  # default: 1, grand canonical mode; else: canonical model
        ):
    '''Set up Gutzwiller calculations for 1-band model with semi-circular DOS.

    Parameters:

    * u: real number
      Hubbard U.
    * nmesh: interger number
      number of energy mesh
    * mu: real number
      chemical potential

    Result:

    Create all the necessary input file of ``GParam.h5``, ``GBareH.h5``
    for *CyGutz* calculation.
    '''
    # get semi-circular class, predifined in pygrisb.model.special
    sc = special.semicircular()

    # get semi-circular energy mesh
    e_list = sc.get_e_list_of_uniform_wt(nmesh=nmesh)

    # get atomic structure with 1 atom per unit cell.
    a = tb.AtomsTB("H", [(0, 0, 0)], cell=(1, 1, 1))

    # specify single s-orbital and orbital dimensions.
    a.set_orbitals_spindeg(orbitals=[('s')])
    norb = 1
    norb2 = norb*2

    # get a list of one-body Hamilonian for the enrgy mesh.
    hk_list = [np.array([[e+0.j]]) for e in e_list]

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

    # list of one-body part of the local Hamiltonians.
    # here the s-orbital is located at zero.
    # input is always complex
    h1e_list = [np.zeros((norb, norb), dtype=np.complex)]

    # create ``GBareH.h5`` file
    gbndinfo = gp.gparam_bnds(kps_wt,
            num_e,
            norb,
            ensemble=ensemble,
            h1e_list=h1e_list,
            delta=0.001)
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
            "db2sab_list": [np.eye(2, dtype=np.complex)],
            "symbol_matrix_list": [np.eye(2, dtype=np.int)],
            "dc_u_avg_list": [0.],
            "dc_j_avg_list": [0.],
            "giembeddiag": iembeddiag,
            "nval_bot_list": [0],
            "nval_top_list": [2],
            })

    # set Coulomb matrix
    v2e = np.zeros((norb2,norb2,norb2,norb2))
    v2e[0,0,0,0] = v2e[0,0,1,1] = v2e[1,1,0,0] = v2e[1,1,1,1] = u
    gatms.update({"v2e_list": [v2e]})
    # local single-particle matrix basis.
    gatms.set_matrix_basis_list(realhemb=True)
    # spin operator
    samop = am.samop_csh2(0.0)
    s_vec = samop.am_op
    gatms.set_s_vec_csh2_list([s_vec])
    gatms.set_sl_vector_sab_list(mode="sz_only")
    gpa = gp.gparam_atms(gatms.ps["na2_list"], gatms.ps)
    gpa.h5save()

    # set the potential shift such that the system has particle-hole
    # symmetry with recpect to zero.
    # also include chemical potential here.
    vext_list = [np.zeros((norb2, norb2))]
    vext_list[0][0, 0] = vext_list[0][1, 1] = -u/2
    gpa.add_vext_list(vext_list, iext=-1)
    # update mu or ef
    with h5py.File("GBareH.h5", "a") as f:
        f["/"].attrs["chempot"] = mu



if __name__=="__main__":
    gutz_model_setup(u=2.)
