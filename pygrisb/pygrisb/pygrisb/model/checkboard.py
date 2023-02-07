import numpy as np
import pickle
from ase.dft import kpoints
import pygrisb.model.tbASE as tb
import pygrisb.gutz.gparam as gp
import pygrisb.gutz.atoms as gatoms
import pygrisb.mbody.coulomb_matrix as cm
import pygrisb.symm.angular_momentum as am


def gutz_model_setup(u=0.0, spindeg=True, num_e=2., iembeddiag=-2):
    '''Set up Gutzwiller calculations for 2d body-centered square lattice.

    Parameters:

    * u: real number
      Hubbard U.
    * spindeg: boolean number
      whether to keep spin degeneracy or not.
    * num_e: real number
      number of electron per unit cell
    * iembeddiag: integer
      flag for method to solve the embedding Hamiltonian.

      * -3: valence truncation ED with S=0 (spin-singlet) constraint;
      * -1: valence truncation ED;
      * 10: Hartree-Fock.

    Result:

    Create all the necessary input file of ``GBareH.h5`` and ``GParam.h5``
    for *CyGutz* calculation.
    '''

    # two "H" atoms in the square unit cell, one at the corner
    # and one at the center.
    symbols=['H', 'H']
    scaled_positions=[(0, 0, 0), (0.5, 0.5, 0)]
    cell = np.identity(3)
    a = tb.AtomsTB(symbols=symbols, scaled_positions=scaled_positions,
            cell=cell)

    # set spin degeneracy accordingly.
    a.set_orbitals_spindeg(spindeg=spindeg)

    # create a tight-binding model class given the AtomsTB.
    aTB = tb.TB(a)

    # set real space (nearest neighbour) hopping elements.
    t = -1.0
    aTB.set_hop([
            (( 0, 0,0),0,1, t),
            ((-1, 0,0),0,1, t),
            (( 0,-1,0),0,1, t),
            ((-1,-1,0),0,1, t),
            (( 0, 0,0),1,0, t),
            (( 1, 0,0),1,0, t),
            (( 0, 1,0),1,0, t),
            (( 1, 1,0),1,0, t),
            ])

    # set 2d k-mesh
    kps_size = (50, 50, 1)
    kps = kpoints.monkhorst_pack(kps_size)

    # set uniform k-point weight
    num_k = len(kps)
    kps_wt = 1.0 / num_k * np.ones((num_k))
    if aTB.Atoms.spindeg:
        kps_wt *= 2
        ispin = 1
    else:
        ispin = 2

    # se maximal number of bands (here we have two bands.)
    num_band_max = 2

    # set list of one-body part of the local Hamiltonian (trivial here.)
    h1e_list = [np.array([[0,]], dtype=np.complex) for symbol in symbols]

    # create ``GBareH.h5`` file
    gbndinfo = gp.gparam_bnds(kps_wt, num_e, num_band_max, h1e_list=h1e_list)
    gbndinfo.h5save()
    aTB.save_bareham(kps)

    # get ``GParam.h5`` file
    gatms = gatoms.cor_structure(symbols, scaled_positions, cell)
    gatms.update({
            "dc_mode": 0,
            "imap_list": [0,1],
            "na2_list": [2,2],
            "iso": 1,
            "ispin": ispin,
            "db2sab_list": [np.eye(1, dtype=np.complex), \
                    np.eye(1, dtype=np.complex)],
            "symbol_matrix_list": [np.eye(1, dtype=np.int), \
                    np.eye(1, dtype=np.int)],
            "dc_u_avg_list": [0., 0.],
            "dc_j_avg_list": [0., 0.],
            "giembeddiag": iembeddiag,
            "nval_bot_list": [0,0],
            "nval_top_list": [2,2],
            })
    gatms.set_full_symm_adaptive_basis()
    # Coulomb matrix
    v2e_1 = np.zeros((1,1,1,1))
    v2e_1[0,0,0,0] = u
    v2e_2 = cm.add_spin_comp_to_coul_mat(v2e_1)
    gatms.update({"v2e_list": [v2e_2, v2e_2]})
    # local single-particle matrix basis.
    gatms.set_matrix_basis_list(realhemb=True)
    # spin operator
    samop = am.samop_csh2(0.)
    s_vec = samop.am_op
    gatms.set_s_vec_csh2_list([s_vec, s_vec])
    gatms.set_sl_vector_sab_list(mode="sz_only")
    gpa = gp.gparam_atms(gatms.ps["na2_list"], gatms.ps)
    gpa.h5save()

    if ispin == 2:
        # introduce symmetry breaking by vext
        vext_list = [np.zeros((2, 2)), np.zeros((2, 2))]
        # afm configuration
        vext_list[0][0, 0] = vext_list[1][1, 1] = -0.1
        vext_list[0][1, 1] = vext_list[1][0, 0] = 0.1
        gpa.add_vext_list(vext_list, iext=0)

    # save the aTB
    with open('aTB.pckl', 'wb') as f:
        pickle.dump([a, aTB], f)


def get_bard_bands():
    '''get bare band structure.
    '''
    # load aTB
    with open('aTB.pckl', 'rb') as f:
        a, aTB = pickle.load(f)

    # k-point path
    kG = [0.0, 0.0, 0]
    kX = [0.5, 0.0, 0]
    kM = [0.5, 0.5, 0]

    # set up a ase.dft.kpoints kpath object
    kps = kpoints.get_bandpath([kG, kX, kM, kG], a.cell)

    # get band structure of a square lattice
    aTB.get_bandstructure(kps, saveto="bare_bands.dat")


if __name__=='__main__':
    gutz_model_setup(num_e=1.)
    get_bard_bands()
