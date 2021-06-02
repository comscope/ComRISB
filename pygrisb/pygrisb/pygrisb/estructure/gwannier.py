import h5py, json, numpy, warnings, os
from scipy.linalg import block_diag
from mpi4py import MPI
from pymatgen.core import Structure
from pymatgen.symmetry.bandstructure import HighSymmKpath
import itertools as it
from pygrisb.iface.wanniertb import w90
from pygrisb.mpi.mpi import get_myk_range as get_myk_range


def get_csh2sab():
    '''get transformation from complex spherical harmonics basis to
    g-risb symmetry adapted basis.
    '''
    csh2sab_list = []
    with h5py.File("GParam.h5", "r") as f:
        imap_list = f["/"].attrs["imap_list"]
        for i, imap in enumerate(imap_list):
            if i == imap:
                csh2sab_list.append(\
                        f[f"/impurity_{i}/db2sab"][()])
            else:
                csh2sab_list.append(csh2sab_list[imap])
    return csh2sab_list


def get_wan2sab():
    '''get transformation from wannier basis to cygutz symmetry adapted basis.
    '''
    csh2sab_list = get_csh2sab()
    with h5py.File("GBareH.h5", "r") as f:
        wan2csh = f["/"].attrs["u_wan2csh"]
    with h5py.File("GParam.h5", "r") as f:
        iso = f["/"].attrs["iso"]

    csh2sab = []
    for u1 in csh2sab_list:
        n1 = u1.shape[0]//(3-iso)
        csh2sab.append(u1[:n1, ::3-iso])
    csh2sab = block_diag(*csh2sab)
    wan2sab = wan2csh.copy()
    n1 = csh2sab.shape[0]
    wan2sab[:,:n1] = wan2csh[:,:n1].dot(csh2sab)
    return wan2sab


def get_gloc_in_wannier_basis():
    '''get the gutzwiller local matrices, including r, lambda, nr, and nphy,
    h1e, from grisb calculation.
    '''
    # read from cygutz output
    with h5py.File("GLog.h5", "r") as f:
        # renormalized local one-body part of the quasiparticle hamiltonian
        lambda_list = f["/"].attrs["LAMAT"].swapaxes(1, 2)
        # quasiparticle renormalization matrix
        r_list = f["/"].attrs["RMAT"].swapaxes(1, 2)
        # local density part to be subtracted in density calculations
        nr_list = f["/"].attrs["NRLMAT"].swapaxes(1, 2)
        # physical density matrix to be added in density calculations
        nphy_list = f["/"].attrs["NPHYMAT"].swapaxes(1, 2)
        # h1e for band structure calculations
        h1e_list = f["/"].attrs["H1E"].swapaxes(1, 2)
        if h1e_list.shape[0] < lambda_list.shape[0]:
            # possibly spin-polarized calculation
            h1e_list = numpy.concatenate((h1e_list, h1e_list))

    # get transformation from wannier basis to cygutz symmetry adapted basis.
    wan2sab = get_wan2sab()
    # convert to wannier basis in each spin block.
    lam2 = []
    r2 = []
    nr2 = []
    nphy2 = []
    h1e2 = []
    n2 = wan2sab.shape[0]
    for rmat, lam, nr, nphy, h1e in zip(r_list,
            lambda_list, nr_list, nphy_list, h1e_list):
        n1 = rmat.shape[0]
        if n2 > n1:
            # rmat adding identity block
            rmat = block_diag(rmat, numpy.eye(n2-n1, dtype=numpy.complex))
            # others adding zeros block
            zeromat = numpy.zeros((n2-n1, n2-n1), numpy.complex)
            lam = block_diag(lam, zeromat)
            nr = block_diag(nr, zeromat)
            nphy = block_diag(nphy, zeromat)
            h1e = block_diag(h1e, zeromat)

        r2.append(wan2sab.dot(rmat).dot(wan2sab.T.conj()))
        lam2.append(wan2sab.dot(lam).dot(wan2sab.T.conj()))
        nr2.append(wan2sab.conj().dot(nr).dot(wan2sab.T))
        nphy2.append(wan2sab.conj().dot(nphy).dot(wan2sab.T))
        h1e2.append(wan2sab.dot(h1e).dot(wan2sab.T.conj()))
    return numpy.asarray(r2), numpy.asarray(lam2), \
            numpy.asarray(nr2), numpy.asarray(nphy2), numpy.asarray(h1e2)


def get_structure():
    data = json.load(open('ginit.json', "r"))
    return Structure(
            lattice=data["struct"]["cell"],
            species=data["struct"]["symbols"],
            coords=data["struct"]["scaled_positions"],
            )


def get_bands(kpoints,
        gmodel=None,
        wfwannier_list=None,
        bnd_es_in=None,
        mode="tb",
        evmode=0,
        #0: bnd_vs is list of eigen-vectrors [v]; otheres: [R^\dag v R]
        ):
    if mode == "risb":
        r_mat, lam_mat, _, _, h1_mat = get_gloc_in_wannier_basis()
        ispin = r_mat.shape[0]
        with h5py.File("GLog.h5", "r") as f:
            efermi = f["/"].attrs["efermi"]
    else:
        ispin = 1
        efermi = 0.
    if evmode != 0:
        # get transformation from wannier basis to cygutz symmetry adapted
        # basis.
        wan2sab = get_wan2sab()

    bnd_es = []
    bnd_vs = []
    for isp in range(ispin):
        if mode == "risb" and wfwannier_list is not None:
            ispp = min(isp, len(wfwannier_list)-1)
        else:
            ispp = min(0, isp)
        bnd_es.append([])
        bnd_vs.append([])
        for ik, kpt in enumerate(kpoints):
            if gmodel is not None:
                hmat = gmodel._gen_ham(kpt,ispp)
            else:
                hmat = wfwannier_list[ispp][ik].T.conj().dot(\
                        numpy.diag(bnd_es_in[ispp][ik])).dot(\
                        wfwannier_list[ispp][ik])
            if mode == "risb":
                hmat -= h1_mat[ispp]
                hmat = r_mat[isp].dot(hmat).dot(r_mat[isp].T.conj())
                hmat += lam_mat[isp]
            evals, evecs = numpy.linalg.eigh(hmat)
            # shift fermi level to zero.
            evals -= efermi
            bnd_es[isp].append(evals)
            if evmode == 0:
                bnd_vs[isp].append(evecs)
            else:
                bnd_vs[isp].append(r_mat[isp].dot(wan2sab).T.conj().dot(evecs))
    bnd_es = numpy.asarray(bnd_es)
    bnd_vs = numpy.asarray(bnd_vs)
    return bnd_es, bnd_vs


def get_symkpath(atol=1.e-6):
    struct = get_structure()
    kpath = HighSymmKpath(struct)
    # check warning and perform transformation if needed.
    if not numpy.allclose(kpath._structure.lattice.matrix,
            kpath.prim.lattice.matrix, atol=atol):
        warnings.warn("Input structure does not match expected standard "
                "primitive! Try k-path transformation.")
        ktrans = kpath.prim.lattice.reciprocal_lattice.matrix.dot(\
                numpy.linalg.inv(kpath._structure.lattice.\
                reciprocal_lattice.matrix))
        for kname in kpath.kpath["kpoints"]:
            kpath.kpath["kpoints"][kname] = \
                    kpath.kpath["kpoints"][kname].dot(ktrans)
    return kpath


def get_gmodel(wpath="../wannier", wprefix="wannier"):
    wannier90 = w90(wpath, wprefix)
    gmodel = wannier90.model()
    return gmodel


def mpiget_bndev(k_list,
        gmodel=None,
        wfwannier_list=None,
        bnd_es_in=None,
        mode="tb",
        ):
    comm = MPI.COMM_WORLD
    nktot = len(k_list)
    iksta, ikend = get_myk_range(nktot)
    kvec_loc = k_list[iksta: ikend]
    # wfwannier_list is always in local k-block
    if bnd_es_in is not None:
        bnd_es_in = bnd_es_in[:, iksta:ikend, :]
    bnd_es, bnd_vs = get_bands(kvec_loc,
            gmodel=gmodel,
            wfwannier_list=wfwannier_list,
            bnd_es_in=bnd_es_in,
            mode=mode,
            )

    # gather bnd_es at root
    bndes_all = comm.gather(bnd_es, root=0)
    if comm.Get_rank() == 0:
        bndes_all = numpy.concatenate(bndes_all, axis=1)
        bnd_es = bndes_all
        assert(bnd_es.shape[1] == nktot), "error in merging bnd_es!"

    return bnd_es, bnd_vs


def get_wannier_den_matrix_risb(bnd_vs, ferwes, wk, nktot):
    r_mat, _, nr_mat, nphy_mat, _ = get_gloc_in_wannier_basis()
    with h5py.File("GBareH.h5", "r") as f:
        ispin_dft = f["/"].attrs["ispin"]
        iso = f["/"].attrs["iso"]
    ispin_risb = r_mat.shape[0]
    # spin factor
    f_ispin = 3-max(iso, ispin_risb)
    wan_den = []
    # total number of electrons to be compared.
    sum_elec1 = numpy.sum(ferwes)
    sum_elec2 = 0.
    for isp in range(ispin_risb):
        if isp < ispin_dft:
            wan_den.append([])
        else:
            wan_den = numpy.asarray(wan_den)
        for ik, bndvk1, ferwek1, wk1 in zip(it.count(),
                bnd_vs[isp], ferwes[isp], wk):
            # notice the convention a bit different from cygutz.
            # <a|psi>f<psi|b>
            afb = bndvk1.dot(numpy.diag(ferwek1/wk1/f_ispin)).dot(\
                    bndvk1.T.conj())
            # R^\dagger_{A,a} * <a|psi>f<psi|b> * R_{b,B}
            rdafbr = r_mat[isp].T.conj().dot(afb).dot(r_mat[isp])
            dmk = rdafbr
            # \rho_{A,B} = R^\dagger_{A,a} * <a|psi>f<psi|b> * R_{b,B}
            #            +(n_phys.^{A,B} - n_{sub.}^{A,B})
            dmk += (nphy_mat[isp]-nr_mat[isp]).T
            sum_elec2 += dmk.trace()*wk1*f_ispin
            if isp < ispin_dft:
                wan_den[-1].append(dmk)
            else:
                wan_den[-1][ik] += dmk
                wan_den[-1][ik] *= 0.5
    sum_elec2 = sum_elec2.real
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    sum_elec_all1 = comm.reduce(sum_elec1)
    sum_elec_all2 = comm.reduce(sum_elec2)
    if rank == 0:
        print(f"sum_ferwt: {sum_elec_all1:.6f} vs kswt: {sum_elec_all2:.6f}!")
        elec_diff = sum_elec_all2 - sum_elec_all1
        if numpy.abs(elec_diff) > 0.1:
            raise ValueError(f"too big charge difference: {elec_diff:.2f}")
    wan_den = numpy.asarray(wan_den)
    # merge wan_den to master node
    wan_den_list = comm.gather(wan_den, root=0)

    if rank == 0:
        wan_den = numpy.concatenate(wan_den_list, axis=1)
        assert(wan_den.shape[1] == nktot), "error in merging wan_den!"

    return wan_den


def get_bands_symkpath(mode="risb", ndiv=24):
    # get tb model
    gmodel = get_gmodel()

    # k-path by high symmetry k-points and labels.
    if os.path.isfile("kpath.dat"):
        kpts = []
        ktick_labels = []
        with open("kpath.dat", "r") as f:
            for i, line in enumerate(f):
                line = line.split()
                if i == 0:
                    # number of points of between high symmetry points.
                    ndiv = int(line[0])
                elif len(line) == 4:
                    kpts.append(list(map(float, line[:3])))
                    ktick_labels.append(line[3])
    else:
        kpath = get_symkpath()
        ktick_labels = list(kpath.kpath["kpoints"].keys())
        kpts = list(kpath.kpath["kpoints"].values())
    kpts = numpy.asarray(kpts)

    print("k-path:")
    for kpt, label in zip(kpts, ktick_labels):
        print(kpt[0], kpt[1], kpt[2], label)

    nk = (len(kpts)-1)*ndiv
    k_vec, k_dist, k_node = gmodel.k_path(
            kpts,
            nk,
            )
    bnd_es, bnd_vs = get_bands(
            k_vec,
            gmodel=gmodel,
            mode=mode,
            evmode=1,  # include R-factor in the eigen-vector
            )
    return bnd_es, bnd_vs, k_dist, k_node, ktick_labels


def get_bandstructure_data():
    bnd_es, bnd_vs, k_dist, k_node, ktick_labels = \
            get_bands_symkpath(mode="risb")
    bnd_vs = bnd_vs.swapaxes(2, 3)
    return k_dist, k_node, ktick_labels, bnd_es, bnd_vs
