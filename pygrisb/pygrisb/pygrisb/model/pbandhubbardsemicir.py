# Author: Yongxin Yao

'''degenerate p-band hubbard model with semi-circular dos.
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
        j=0.0,
        nmesh=5000,
        mu=0,
        iembeddiag=-2,
        num_e=3.,  # electron filling, not referenced for grand canonical model
        ensemble=0,  # default: 1, grand canonical mode; else: canonical model
        ):
    '''Set up Gutzwiller calculations for p-band hubbard model
    with semi-circular DOS.

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
    a = tb.AtomsTB("C", [(0, 0, 0)], cell=(1, 1, 1))

    # specify single s-orbital and orbital dimensions.
    a.set_orbitals_spindeg(orbitals=[('px', 'py', 'pz')])
    norb = 3
    norb2 = norb*2

    # get a list of one-body Hamilonian for the enrgy mesh.
    hk_list = [np.identity(norb)*e for e in e_list]

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
    gatms = gatoms.cor_structure(symbols=["C"], scaled_positions=[[0,0,0]],
            cell=np.identity(3))
    from pygrisb.symm.unitary import complx_sph_harm_to_real_harm
    csh2rh = complx_sph_harm_to_real_harm(1).u_trans
    db2sab = np.zeros((norb2, norb2), dtype=np.complex)
    db2sab[:norb, ::2] = db2sab[norb:, 1::2] = csh2rh
    gatms.update({
            "dc_mode": 1,
            "imap_list": [0],
            "na2_list": [norb2],
            "iso": 1,
            "ispin": 1,
            "db2sab_list": [db2sab],
            "symbol_matrix_list": [np.eye(norb2, dtype=np.int)],
            "dc_u_avg_list": [u],
            "dc_j_avg_list": [j],
            "giembeddiag": iembeddiag,
            "nval_bot_list": [0],
            "nval_top_list": [norb2],
            })

    # set Coulomb matrix
    import pygrisb.mbody.coulomb_matrix as c_matrix
    cm = c_matrix.coulomb_matrix_slater.from_uj(1, u, j)
    cm.unitary_transform(db2sab)
    v2e = cm.u_spin_orb.real
    if np.max(np.abs(cm.u_spin_orb.imag)) > 1.e-6:
        raise AssertionError("real version with complex v2e!")
    gatms.update({"v2e_list": [v2e]})
    # local single-particle matrix basis.
    gatms.set_matrix_basis_list(realhemb=True)
    # spin operator
    samop = am.samop_csh2(1.0)
    s_vec = samop.am_op
    gatms.set_s_vec_csh2_list([s_vec])
    gatms.set_sl_vector_sab_list(mode="sz_only")
    gpa = gp.gparam_atms(gatms.ps["na2_list"], gatms.ps)
    gpa.h5save()


if __name__=="__main__":
    gutz_model_setup(u=2.)
