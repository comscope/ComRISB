import numpy as np
import pickle
from ase.dft import kpoints
import pygrisb.gutz.gparam as gp
import pygrisb.gutz.atoms as gatoms
import pygrisb.model.tbASE as tb
import pygrisb.mbody.coulomb_matrix as cm
import pygrisb.symm.angular_momentum as am


def gutz_model_setup(
        nelect = 1.933878,
        u=12.0,
        joveru=0.2,
        spindeg=True,
        iembeddiag=-2,
        delta=0.01,
        kps_size = (50, 50, 1),
        dtype=np.float,
        mu=0,   #1.45
        ensemble=0,
        ):
    '''Set up Gutzwiller calculations for (dxz, dyz) 2D model
    (PRB 77, 220503, R2008) in square lattice.

    Parameters:

    * u: real number
      Hubbard U.
    * spindeg: boolean number
      whether to keep spin degeneracy or not.
    * iembeddiag: integer
      flag for method to solve the embedding Hamiltonian.

      * -3: valence truncation ED with S=0 (spin-singlet) constraint;
      * -2: valence truncation ED with Sz symmetry;
      * -1: valence truncation ED;
      * 10: Hartree-Fock.

    Result:

    Create all the necessary input file of ``GBareH.h5`` and ``GParam.h5``
    for *CyGutz* calculation.
    '''

    # One pseudo "H" atoms in the square unit cell
    symbols=['Fe']
    scaled_positions=[(0, 0, 0)]
    cell = np.identity(3)
    a = tb.AtomsTB(symbols=symbols, scaled_positions=scaled_positions,
            cell=cell)

    # set spin degeneracy accordingly.
    a.set_orbitals_spindeg(orbitals=[("xz","yz")], spindeg=spindeg)

    # create a tight-binding model class given the AtomsTB.
    aTB = tb.TB(a)

    # set real space (nearest neighbour) hopping elements.
    t1 = -1.; t2 = 1.3; t3 = t4 = -0.85
    # covention difference
    t1 = -t1; t2 = -t2; t3 = -t3; t4 = -t4

    aTB.set_hop([
            (( 0, 0, 0), 0, 0, -mu),
            (( 0, 0, 0), 1, 1, -mu),

            (( 1, 0, 0), 0, 0, t1),
            (( 1, 0, 0), 1, 1, t2),
            ((-1, 0, 0), 0, 0, t1),
            ((-1, 0, 0), 1, 1, t2),
            (( 0, 1, 0), 0, 0, t2),
            (( 0, 1, 0), 1, 1, t1),
            (( 0, -1, 0), 0, 0, t2),
            (( 0, -1, 0), 1, 1, t1),

            (( 1, 1, 0), 0, 0, t3),
            (( 1, 1, 0), 1, 0, t4),
            (( 1, 1, 0), 0, 1, t4),
            (( 1, 1, 0), 1, 1, t3),
            ((-1, -1, 0), 0, 0, t3),
            ((-1, -1, 0), 0, 1, t4),
            ((-1, -1, 0), 1, 0, t4),
            ((-1, -1, 0), 1, 1, t3),

            (( 1, -1, 0), 0, 0, t3),
            (( 1, -1, 0), 0, 1, -t4),
            (( 1, -1, 0), 1, 0, -t4),
            (( 1, -1, 0), 1, 1, t3),
            (( -1, 1, 0), 0, 0, t3),
            (( -1, 1, 0), 0, 1, -t4),
            (( -1, 1, 0), 1, 0, -t4),
            (( -1, 1, 0), 1, 1, t3),
            ])

    # set 2d k-mesh
    kps = kpoints.monkhorst_pack(kps_size)

    # set uniform k-point weight
    num_k = len(kps)
    kps_wt = 1.0 / num_k * np.ones((num_k))
    if aTB.Atoms.spindeg:
        kps_wt *= 2

    # set maximal number of bands (here we have two bands.)
    num_band_max = 2

    # set list of one-body part of the local Hamiltonian (trivial here.)
    # always complex
    h1e_list = [np.array([[0, 0], [0, 0]], dtype=np.complex)
            for symbol in symbols]

    # create ``GBareH.h5`` file
    gbndinfo = gp.gparam_bnds(kps_wt,
            nelect,
            num_band_max,
            h1e_list=h1e_list,
            delta=delta,
            ensemble=ensemble,
            )
    gbndinfo.h5save()

    # adding bare Hamiltonian to ``GBareH.h5`` file.
    aTB.save_bareham(kps)

    # get ``GParam.h5`` file
    gatms = gatoms.cor_structure(symbols, scaled_positions, cell)
    gatms.update({
            "dc_mode": 0,
            "imap_list": [0],
            "na2_list": [4],
            "iso": 1,
            "ispin": 1,
            "db2sab_list": [np.eye(2, dtype=dtype)],
            "symbol_matrix_list": [np.eye(2)],
            "dc_u_avg_list": [0.],
            "dc_j_avg_list": [0.],
            "giembeddiag": iembeddiag,
            "nval_bot_list": [0],
            "nval_top_list": [4],
            })
    gatms.set_full_symm_adaptive_basis()
    # Coulomb matrix
    v2e = cm.coulomb_matrix_kanamori(2, u, joveru*u)
    v2e.evaluate_all()
    v2e_2 = cm.apply_a4_unitary_trans(v2e.u_spin_orb,
            gatms.ps["db2sab_list"][0])
    v2e_2 = np.asarray(v2e_2, dtype=dtype)
    gatms.update({"v2e_list": [v2e_2]})
    # local single-particle matrix basis.
    gatms.set_matrix_basis_list(realhemb=(dtype==np.float))
    # spin operator
    samop = am.samop_csh2(0.5)
    s_vec = samop.am_op
    gatms.set_s_vec_csh2_list([s_vec])
    gatms.set_sl_vector_sab_list(mode="sz_only")
    gpa = gp.gparam_atms(gatms.ps["na2_list"], gatms.ps)
    gpa.h5save()

    # save the aTB
    with open('aTB.pckl', 'wb') as f:
        pickle.dump([a, aTB], f)


def get_bare_bands():
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
    kps = kpoints.bandpath([kG, kX, kM, kG], a.cell, npoints=100)


    # get band structure of a square lattice
    aTB.get_bandstructure(kps, saveto="bare_bands.dat")


def get_bare_dos():
    '''get bare band structure.
    '''
    # load aTB
    with open('aTB.pckl', 'rb') as f:
        a, aTB = pickle.load(f)

    # get dos of a square lattice
    aTB.get_dos(kps_size=[500, 500, 1], saveto="bare_dos.dat")



if __name__=='__main__':
    gutz_model_setup()
    get_bare_bands()
    # get_bare_dos()
