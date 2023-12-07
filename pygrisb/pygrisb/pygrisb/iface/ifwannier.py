from mpi4py import MPI
import numpy as np
import h5py
from pygrisb.symm.unitary import get_u_csh2wan_all
from scipy.linalg import block_diag
from scipy.io import FortranFile
from pygrisb.estructure.gwannier import mpiget_bndev, \
        get_wannier_den_matrix_risb
from pygrisb.estructure.fermi import get_fermi_weight, get_fermi_level
from pygrisb.iface.wannierio import mpiget_wannier_data

# if_gwannier reads info from wannier.dat, and also gets  the number of electrons from  wannier_0.out.
# if_gwannier saves the following info to GBareH.h5:
# 'iso' = iso # if control['spin_orbit'] then iso=2, otherwise iso=1
# 'ispin' = ispin # incoming value is 1 - paramagnetic DFT calculation, 
#                                     2 for spin-polarized DFT calculation (maybe not supported)
# 'kptdim' = numk = the total number of k-points    
# 'nbmax' = nbmax = the number of bands. maximal number of band wannier funcitons across all the k-points.
# 'kptwt' = wk_list is a list of the k-point weights
# ['kptx', 'kpty', 'kptz'] = kpts is a list of all the k-points
#   brillouin zone integration method: fermi or gaussian smearing
# 'ismear' = ismear # control['ismear'] This specifies the 
#       integration method: fermi or gaussian smearing.
#       gaussian smearing (ismear=0) or Fermi smearing (-1) or 
#           tetra-hedron method (-5) will be used for the brillouin zone integration.
# 'delta' = delta = inverse temperature,  control['delta']
# 'nelectron' = nelectron is an electron count: sum over momenta of the
#   Fermi-Dirac distribution, multiplied by the k-point weight
#   and the spin factor spin_deg [1 if spin-orbit, 2 otherwise]
#   Then it is rounded: nelectron = int(nelectron+0.5)
#   nelectron is the number of valence electrons in the wannier manifold, and in
#       the model which is solved by CyGUTZ.
# 'symnop' = 1 is the number of symmetry operations
# 'SYMIE' = 1 = the index of identity symmetry operation
#                   the index of the identity operator in the list of symmetry operations.
# 'mode_hk' = 0
# 'u_wan2csh' = u_wan2csh  is a  matrix that sorts orbitals into first wannier orbitals and then non-wannier orbitals.
#       In the wannier portion of u_wan2csh, info is included about the basis
#       transformation from real spherical harmonics to complex spherical harmonics.
# The following have an extra index isp, for spin:
# 'ispin_{isp}/H1E_LIST' = h1e_array
#       h1e_list = hmat / k-point weight, then truncated down to the 
#           wannier part of the basis.
#       one-body part (10 x 10 matrix for d-orbital) of the local Hamiltonian
# 'ispin_{isp}/ne_list' = ne_list[isp]
#       ne_list is the band indices associated with wannier orbital construction.
#           NE_LIST[:,1]: total number of bands at each k-point
#           NE_LIST[:,2/3]: the first/last band WF included to expand the local ‘’correlated’’ orbitals.
# 'ispin_{isp}/e_list' = evals_list[isp]
# 'ispin_{isp}/U_PSIK0_TO_HK0_BASIS_LIST' = psi0b_list[isp] with the axis swapped
#       evals_list=evals,psi0b_list=evecs are the  eigenvalues and eigenvectors of hmat
#           hmat_list  = U^\dagger W^\dagger E W U, where:
#               - W = wannier_list is a list of overlap between band wavefunctions and 
#                       wannier orbitals.
#               - U = u_wan2csh transforms spherical harmonics and re-orders orbitals
#               - E = bnd_es is a list of band energies     
def if_gwannier(corbs_list, # control['impurity_wan']
        delta_charge=0., # incoming valaue is 0.
        wpath="../wannier", # control['wannier_directory']
        lpath="../lattice", # incoming value is "../lattice"
        wprefix="wannier", # incoming value is "wannier"
        lprefix="case", # control['allfile']
        iso=1, # if control['spin_orbit'] then iso=2, otherwise iso=1
        ispin=1, # incoming value is 1
        ismear=0, # control['ismear'] This specifies the brillouin zone 
                    # integration method: fermi or gaussian smearing
        delta=0.02585,  # 300k in eV # control['delta']
        method="lda+risb", # control['method']
        icycle=0, # icycle
        ):
    
# mpiget_wannier_data calls get_wannier_data_binary, which reads data
#    out of wannier.dat . For more info look at the comments on get_wannier_data_binary.
#   It broadcasts reals_lat, recip_lat, kpts, bnd_es, include_bands to all mpi processes.
#   It also splits kpts and wfwannier_list between mpi processes.
#       After the split, the local process' k-points are stored in k_range,
#       and its share of wfwannier_list is still called wfwannier_list.
#   -----------
#   - reals_lat is the real space lattice vectors
#   - recip_lat is the reciprocal space lattice vectors
#   - kpts is a list of all the k-points
#   - include_bands is low energy band indices included in the wannier construction
#   - wfwannier_list is a list of overlap between band wavefunctions and 
#           wannier orbitals [This is split over mpi processes.]
#   - bnd_es list of band energies
#   - k_range is a list of the k-points that are assigned to the current process.
# Things read but not returned directly:
#   - num_bands is the maximal number of bands (controls the size of wfwannier_list,
#            bnd_es)
#   - num_wann is the number of wannier orbitals ( controls the size of wfwannier_list)
#   - ndiv is the k-point mesh (not returned at all, but the number of k-points
#           controls  the size of kpts, wfwannier_list, bnd_es)
    _, _, kpts, include_bands, wfwannier_list, bnd_es, k_range = \
            mpiget_wannier_data(path=wpath)
            
    dft_calc = 'lda' in method or 'dft' in method
    if dft_calc:
        # total number of valence electrons
        # get_total_valence_elec gets the number of electrons by grepping through
        #   file lpath/lprefix_0.out=lattice/wannier_0.out, 
        #   broadcasts it to all nodes, and returns it
        n_elec = get_total_valence_elec("{}/{}_0.out".format(lpath, lprefix))
        
        n_elec -= max(0, (include_bands[0]-1)*(3-iso))

   
    # numk = the total number of k-points    
    numk = kpts.shape[0]
    wk = 1./numk # weight of each k-point
    
    # nbmax = the number of bands
    nbmax = wfwannier_list.shape[3]
    
    # spin degeneracy spin_deg = 1 if spin-orbit, 2 otherwise
    # if control['spin_orbit'] then iso=2, otherwise iso=1
    # ispin: incoming value is 1
    spin_deg = 3 - max(ispin, iso)

    comm = MPI.COMM_WORLD
    myrank = comm.Get_rank()
    if myrank == 0:
        print(f"corbs_list: {corbs_list}")
        
    # get_wann2csh creates u_wan2csh, a matrix.
    # u_wan2csh sorts orbitals into first wannier orbitals and then non-wannier orbitals.
    # In the wannier portion of u_wan2csh, info is included about the basis
    #   transformation from real spherical harmonics to complex spherical harmonics.
    # corbs_list comes from control['impurity_wan']
    # nbmax is the number of bands
    u_wan2csh = get_wann2csh(nbmax, corbs_list)
    
    h1e_all = []
    nelectron = 0
    eband = 0
    
    # input file of dft bare band structure information for cygutz
    hmat_list = []
    evals_list = []
    psi0b_list = []
    for isp in range(wfwannier_list.shape[0]):
        
        hmat_list.append([])
        evals_list.append([])
        psi0b_list.append([])
        h1e_all.append(np.zeros((nbmax, nbmax), dtype=np.complex))
        
        # go through the k-points in the current mpi process
        for ik in range(k_range[0], k_range[1]):
            ikl = ik-k_range[0]
            
            # rescontruct dft hamiltonian matrix
            # in the wannier basis.
            # since it involves a downfolding,
            # some information outside of the frozen energy window
            # will be lost.
            hmat = wfwannier_list[isp][ikl].T.conj().dot( \
                    np.diag(bnd_es[isp][ik])).dot( \
                    wfwannier_list[isp][ikl])
            # from wannier basis to correlated orbital-ordered
            # complex spherical harmonics basis,
            # which is the convention used in the cygutz
            # initialization script.
            hmat = u_wan2csh.T.conj().dot(hmat).dot(u_wan2csh)
            
            # record the onsite one-body part
            h1e_all[isp] += hmat*wk
            hmat_list[isp].append([hmat.copy()])
            
            # get the eigen-value and eigen-vectors
            evals, evecs = np.linalg.eigh(hmat)
            # another way to evaluate total valence electrons
            # according to sangkook.
            ferwes = 1/(np.exp(evals/delta) + 1)*wk*spin_deg 
            nelectron += sum(ferwes)
            
            # eband is the eigenvalue, multiplied by the Fermi-Dirac distribution,
            #       the k-point weight, and the spin factor spin_deg
            eband += sum(ferwes*evals)
            
            # yes, here it is redundant here
            # but for the sake of consistent with wien2k interface.
            # here it implies the downfolding procedure is not necessary.
            evals_list[isp].append(evals)
            psi0b_list[isp].append(evecs)

            
    nelectron = comm.reduce(nelectron, op=MPI.SUM)
    eband = comm.reduce(eband, op=MPI.SUM)
    h1e_all = np.asarray(h1e_all)
    h1e_all = comm.reduce(h1e_all, op=MPI.SUM)
    evals_list = np.asarray(evals_list)
    hmat_list = np.asarray(hmat_list)
    psi0b_list = np.asarray(psi0b_list)
    hmat_list = comm.gather(hmat_list, root=0)
    evals_list = comm.gather(evals_list, root=0)
    psi0b_list = comm.gather(psi0b_list, root=0)

    if myrank == 0: # If on root node, save to GBareH.h5
        
        hmat_list = np.concatenate(hmat_list, axis=1)
        evals_list = np.concatenate(evals_list, axis=1)
        psi0b_list = np.concatenate(psi0b_list, axis=1)
        
        h1e_list = []
        for isp in range(wfwannier_list.shape[0]):
            h1e_list.append([])
            base = 0
            for corbs in corbs_list:
                norbs = len(corbs)
                h1e_list[isp].append(h1e_all[isp][base:base+norbs, \
                        base:base+norbs])
                base += norbs
                
        print(f"estimated wannier valence electrons: {nelectron:.4f}")
        if dft_calc:
            print(f"    exact valence electrons: {n_elec:.4f}")
            if np.abs(nelectron - n_elec) > 0.1:
                warnings.warn("valence electron number inconsistent!")
        nelectron = int(nelectron+0.5)
        
        # print out eband
        print(f"wannier band energy  = {eband:.6f}")

        ne_list = []
        for _ in h1e_list:
            # zero-based indices: total, istart, iend
            ne_list.append([[nbmax, 0, nbmax-1] for k in range(numk)])
        wk_list = [wk for k in range(numk)]
        nelectron += delta_charge

            # numk = the total number of k-points    
            # nbmax = the number of bands
            # ne_list is the band indices associated with wannier orbital construction.
            # wk_list is a list of the k-point weights
            # kpts is a list of all the k-points
            # nelectron is an electron count: sum over momenta of the
            #   Fermi-Dirac distribution, multiplied by the k-point weight
            #   and the spin factor spin_deg [1 if spin-orbit, 2 otherwise]
            #   Then it is rounded: nelectron = int(nelectron+0.5)
            # h1e_list = hmat / k-point weight, then truncated down to the 
            #   wannier part of the basis.
            # hmat_list  = U^\dagger W^\dagger E W U, where:
            #   - W = wannier_list is a list of overlap between band wavefunctions and 
            #           wannier orbitals.
            #   - U = u_wan2csh transforms spherical harmonics and re-orders orbitals
            #   - E = bnd_es is a list of band energies     
            # evals_list=evals,psi0b_list=evecs are the  eigenvalues and eigenvectors of hmat
            # u_wan2csh  is a  matrix that sorts orbitals into first wannier orbitals and then non-wannier orbitals.
            #       In the wannier portion of u_wan2csh, info is included about the basis
            #       transformation from real spherical harmonics to complex spherical harmonics.
            # iso # if control['spin_orbit'] then iso=2, otherwise iso=1
            # ispin # incoming value is 1
            # ismear # control['ismear']
            # delta = inverse temperature,  control['delta']
            
        h5wrt_gbareh(numk,
                nbmax,
                ne_list,
                wk_list,
                kpts,
                nelectron,
                h1e_list,
                hmat_list,
                evals_list,
                psi0b_list,
                u_wan2csh,
                iso=iso,
                ispin=ispin,
                ismear=ismear,
                delta=delta,
                )

# get_total_valence_elec gets the number of electrons by grepping through
#   file fname, broadcasts it to all nodes, and returns it
def get_total_valence_elec(fname):
    comm = MPI.COMM_WORLD
    myrank = comm.Get_rank()
    if myrank == 0:
        with open(fname, "r") as f:
            for line in f.readlines():
                if "valence charge in whole" in line:
                    n_elec = float(line.split()[-1])
                    break
    else:
        n_elec = None
    n_elec = comm.bcast(n_elec, root=0)
    return n_elec

# never used
def get_kvec_para_index(numk):
    comm = MPI.COMM_WORLD
    ncpu = comm.Get_size()
    nk_per_cpu = numk//ncpu
    if nk_per_cpu*ncpu < numk:
        nk_per_cpu += 1
    myrank = comm.Get_rank()
    k_start = nk_per_cpu*myrank
    nk_loc = min(nk_per_cpu, numk-k_start)
    return k_start, nk_loc

# never used
def get_h1e_list_wannier(gwannier, corbs_list):
    # List of one-body parts of Hamiltonian.
    h1e_list = [[]]
    for i, corbs in enumerate(corbs_list):
        norbs = len(corbs)
        h1e_list[0].append(np.zeros((norbs,norbs), dtype=np.complex))
        for j1 in range(norbs):
            _j1 = corbs[j1]
            for j2 in range(norbs):
                _j2 = corbs[j2]
                h1e_list[0][-1][j1,j2] = gwannier.ham_r[(0,0,0)]["h"][_j1,_j2]\
                        /float(gwannier.ham_r[(0,0,0)]["deg"])
    # Unitary transformation from complex Harmonics to wannier.
    u_csh2wan_list = get_u_csh2wan_all([len(corbs) for corbs in corbs_list])

    for i, u_csh2wan in enumerate(u_csh2wan_list):
        h1e_list[0][i] = u_csh2wan.dot(h1e_list[0][i]).dot(u_csh2wan.T.conj())
    return h1e_list

# get_wann2csh creates u_wan2csh, a matrix.
# u_wan2csh sorts orbitals into first wannier orbitals and then non-wannier orbitals.
# In the wannier portion of u_wan2csh, info is included about the basis
#   transformation from creal spherical harmonics to omplex spherical harmonics.
# corbs_list comes from control['impurity_wan'], and is a list of the 
#   orbitals that are wanniers
# nbmax is the number of bands
def get_wann2csh(nbmax, corbs_list):
    '''
    get transformation from wannier basis to complex spherical harmonics basis.
    '''
    
    # ubasis is a basis transformation matrix which moves the corelated orbitals
    # to the first places, followed by the uncorrelated orbitals.
    ubasis = np.zeros((nbmax,nbmax), dtype=np.complex)
    
    # the spin-orbital index remapping accordingly.   
    # this puts the wannier orbitals together
    orbs_map = []
    for i, corbs in enumerate(corbs_list):
        orbs_map.extend(corbs)
    # add in the orbitals which are not orbitals       
    # appending the uncorrelated orbitals
    for i in range(nbmax):
        if i not in orbs_map:
            orbs_map.append(i)
            
    # create ubasis: <original basis | correlated orbs-first basis>
    for i in range(nbmax):
        ubasis[orbs_map[i],i] = 1.

    # Unitary transformation from complex Harmonics to wannier.
    
    
    #    Transformations from complex spherical harmonics to real
    #    spherical harmonics transformation.
    #    Condon–Shortley phase convention for complex spherical harmonics.
    #    See https://en.wikipedia.org/wiki/Spherical_harmonics.
    u_csh2wan_list = get_u_csh2wan_all([len(corbs) for corbs in corbs_list])
    
    # get the transformation from wannier basis to correlated
    # orbital-ordered complex spherical Harmonics basis.
    
    # block_diag creates a block diagonal matrix from u_csh2wan_list.
    u_csh2wan = block_diag(*u_csh2wan_list)
    
    ncorbs = u_csh2wan.shape[0]
    
    u_wan2csh = ubasis.copy()
    # u_wan2csh sorts orbitals into first wannier orbitals and then non-wannier orbitals.
    # In the wannier portion of u_wan2csh, info is included about the basis
    #   transformation from creal spherical harmonics to omplex spherical harmonics.
    u_wan2csh[:,:ncorbs] = ubasis[:,:ncorbs].dot(u_csh2wan.T.conj())
    
    return u_wan2csh


def h5wrt_gbareh(numk,
        nbmax,
        ne_list,
        wk_list,
        kpoints,
        nelectron,
        h1e_list,
        hk_list,
        evals_list,
        psi0b_list,
        u_wan2csh,
        iso=1,
        ispin=1,
        ismear=0,
        delta=0.02585,  # ev unit, room temperature.
        ):
    # single file for the dft band structure information.
    with h5py.File('GBareH.h5', 'w') as f:
        # spin-orbit coupling
        f['/'].attrs['iso'] = iso
        # spin
        f['/'].attrs['ispin'] = ispin
        # k-points dimension
        f['/'].attrs['kptdim'] = numk
        # maximal number of bands
        f['/'].attrs['nbmax'] = nbmax
        # k-points weight
        f['/'].attrs['kptwt'] = wk_list
        # k-points
        f['/'].attrs['kptx'] = kpoints[:,0]
        f['/'].attrs['kpty'] = kpoints[:,1]
        f['/'].attrs['kptz'] = kpoints[:,2]
        # brillouin zone integration method: fermi or gaussian smearing
        f['/'].attrs['ismear'] = ismear
        # smearing factor
        f['/'].attrs['delta'] = delta
        # number of valence electrons in the wannier manifold
        f['/'].attrs['nelectron'] = nelectron
        # number symmetry operations
        f['/'].attrs['symnop'] = 1
        # the index of identity symmetry operation
        f['/'].attrs['SYMIE'] = 1
        f['/'].attrs['mode_hk'] = 0
        f['/'].attrs['u_wan2csh'] = u_wan2csh
        for isp, h1es in enumerate(h1e_list):
            # merge h1es to array, patched by 0 if needed.
            nmax = np.max([h.shape[0] for h in h1es])
            h1e_array = np.zeros((len(h1es), nmax, nmax), \
                    dtype=h1es[0].dtype)
            for i, h in enumerate(h1es):
                n = h.shape[0]
                h1e_array[i][:n,:n] = h.T
            f[f"ispin_{isp}/H1E_LIST"] = h1e_array
            # band indices associated with local orbital construction.
            f[f'ispin_{isp}/ne_list'] = ne_list[isp]
            f[f'ispin_{isp}/e_list'] = evals_list[isp]
            f[f'ispin_{isp}/U_PSIK0_TO_HK0_BASIS_LIST'] = psi0b_list[isp].\
                    swapaxes(1, 2)
            f[f'ispin_{isp}/HK0_LIST'] = hk_list[isp]

# wannier_den_matrix  reads from wannier.dat, GLog.h5, GParam.h5, GBareH.h5.
# wannier_den_matrix writes out wannier_den_matrix.dat
# If file fwannier is present on disk and is an hdf5 file, then it writes out an hdf5 file.
# Otherwise it writes out a binary file.
# wannier_den_matrix contains the following info:
#   'temperature' = 300.0
#   'nqdiv' = the number of k-points, obtained from the shape of wan_den
#   'num_band_in' = the number of bands, obtained from the shape of 
#       wfwannier_list, which is a list of overlap between band wavefunctions and 
#           wannier orbitals, and is obtained from wannier.dat.
#   'n_dmft_wan' = the number of wanniers, obtained from the shape of wfwannier_list
#   'include_band_new' = include_bands, which is the low energy band indices 
#           included in the wannier construction, and is from wannier.dat.
#   'vw_matrix_new' = wf_wannier_list[0]
#   'n_matrix_new' = wan_den[0], with the axes swapped.
#       wan_den is the density matrix resulting from the CyGUTZ calculation.
#       For more information about how wan_den is calculated, see the comments below.
def wannier_den_matrix(wannier_path="./"):
    '''produce the file `wannier_den_matrix.dat` for the feedback from
    g-risb to dft.
    '''
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    # mpiget_wannier_data calls get_wannier_data_binary, which reads data
    #    out of wannier.dat . For more info look at the comments on get_wannier_data_binary.
    #   It broadcasts reals_lat, recip_lat, kpts, bnd_es, include_bands to all mpi processes.
    #   It also splits kpts and wfwannier_list between mpi processes.
    #       After the split, the local process' k-points are stored in k_range,
    #       and its share of wfwannier_list is still called wfwannier_list.
    #   -----------
    #   - reals_lat is the real space lattice vectors
    #   - recip_lat is the reciprocal space lattice vectors
    #   - kpts is a list of all the k-points
    #   - include_bands is low energy band indices included in the wannier construction
    #   - wfwannier_list is a list of overlap between band wavefunctions and 
    #           wannier orbitals [This is split over mpi processes.]
    #   - bnd_es list of band energies
    #   - k_range is a list of the k-points that are assigned to the current process.
    # Things read but not returned directly:
    #   - num_bands is the maximal number of bands (controls the size of wfwannier_list,
    #            bnd_es)
    #   - num_wann is the number of wannier orbitals ( controls the size of wfwannier_list)
    #   - ndiv is the k-point mesh (not returned at all, but the number of k-points
    #           controls  the size of kpts, wfwannier_list, bnd_es)
    _, _, kpts, include_bands, wfwannier_list, bnd_es_in, k_range = \
            mpiget_wannier_data(path=wannier_path)
 
    # get_bands returns the eigenvalues and eigenvectors of hmat, the Gutzwiller
    #   hamiltonian. The eigenvalues first have the Fermi level (from GLog.h5) subtracted from them.      
    # hmat = R (W^\dagger E W - H1) R^T - Lambda
    #   W = wfwannier_list is a list of overlap between band wavefunctions and 
    #           wannier orbitals.
    #   E = bnd_es_in is a list of band energies.
    #   R = "LAMAT" = renormalized local one-body part of the quasiparticle hamiltonian, from GLog.h5
    #   Lambda=  "RMAT" = quasiparticle renormalization matrix, from GLog.h5
    #   H1 = "H1E" = h1e for band structure calculations, from GLog.h5              
    #   R, Lambda, H1 are transformed from the cygutz symmetry adapted basis
    #       to the wannier basis, using information pulled from GBareH.h5 and GParam.h5.
    # ----------
    # k_list is a slice of the list of k-points
    # gmodel._gen_ham can substitute for the wannier-based hamiltonian
    # wfwannier_list is a list of overlap between band wavefunctions and 
    #           wannier orbitals [This is split over mpi processes.]
    # bnd_es_in is a list of eigenvalues
    # mode can be "tb" or "risb" 
    # evmode, if not zero, transforms the returned eigenvectors from  the wannier basis
    #   to the cygutz symmetry adapted basis                
    bnd_es, bnd_vs = mpiget_bndev(kpts, wfwannier_list=wfwannier_list,
            bnd_es_in=bnd_es_in, mode="risb")
    
    nktot = len(kpts)

    with h5py.File("GBareH.h5", "r") as f:
        delta = f["/"].attrs["delta"]
        ismear = f["/"].attrs["ismear"]
        iso = f["/"].attrs["iso"]
        num_elec = f["/"].attrs["nelectron"]

    # correction due to mott solution.
    # ne_mott = sum_impurities GParam.h5/mott/impurity_{i}/num_mott_electrons
    # no_mott = sum_impurities GParam.h5/mott/impurity_{i}/num_mott_orbitals
    # num_elec = GBareH.h5/nelectron - ne_mott
    ne_mott = 0.
    no_mott = 0
    with h5py.File("GParam.h5", "r") as f:
        if "/mott" in f:
            num_imp = f["/"].attrs["num_imp"]
            for i in range(num_imp):
                ne_mott += f[f'/mott/impurity_{i}'].attrs['num_mott_electrons']
                no_mott += f[f'/mott/impurity_{i}'].attrs['num_mott_orbitals']
                
    # num_elec = the number of electrons in the wannier space, minus a mott correction
    # num_elec = GBareH.h5/nelectron - ne_mott
    num_elec -= ne_mott

    
    # set wk_list, a list of k-point weights
    wklist = [1./nktot for i in range(bnd_es.shape[1])]
    
    if rank == 0:
        # get_fermi_level computes and returns the Fermi level, which is defined
        #   as the energy which solves num_e = \sum get_fermi_weight.  The fermi weight
        #   may be governed by a Fermi-Dirac distribution if ismear=-1, or a 
        #   Gaussian distribution if ismear=0. 
        # bnd_es = list of band energies
        # wklist = list of k-point weights
        # num_elec = the number of electrons in the wannier space, minus a mott correction
        # delta = the inverse of the temperature
        # 'ismear' = ismear # control['ismear'] This specifies the 
        #       integration method: fermi or gaussian smearing.
        #       gaussian smearing (ismear=0) or Fermi smearing (-1) or 
        #           tetra-hedron method (-5) will be used for the brillouin zone integration.
        #   brillouin zone integration method: fermi or gaussian smearing
        # 'iso' = iso # if control['spin_orbit'] then iso=2, otherwise iso=1
        efermi = get_fermi_level(bnd_es, wklist, num_elec, delta=delta, \
                ismear=ismear, iso=iso)
    else:
        efermi = None
    efermi = comm.bcast(efermi, root=0)
    
    # reduce bnd_es to local part only for root rank,
    # other ranks already hold only the local part.
    if rank == 0:
        print(f"fermi level: {efermi:.5f} expect 0.")
        ncpu = comm.Get_size()
        nk_cpu = nktot//ncpu
        if nk_cpu*ncpu < nktot:
            nk_cpu += 1
        bnd_es = bnd_es[:nk_cpu]
        wklist = wklist[:nk_cpu]
        
    # get_fermi_weight computes and returns the Fermi weight.
    #   The Fermi weight is defined individually as a function of orbital (including
    #       Mott orbitals and also bands), spin, and k-point.
    #  The fermi weight
    #   may be governed by a Fermi-Dirac distribution if ismear=-1, or a 
    #   Gaussian distribution if ismear=0. 
    #   For Mott orbitals, the Fermi weight is defined as ne_mott/no_mott.
    # efermi = the Fermi energy
    # bnd_es = list of band energies
    # wklist = list of k-point weights
    # delta = the inverse of the temperature
    # 'ismear' = ismear # control['ismear'] This specifies the 
    #       integration method: fermi or gaussian smearing.
    #       gaussian smearing (ismear=0) or Fermi smearing (-1) or 
    #           tetra-hedron method (-5) will be used for the brillouin zone integration.
    #   brillouin zone integration method: fermi or gaussian smearing
    # 'iso' = iso # if control['spin_orbit'] then iso=2, otherwise iso=1   
    # no_mott = sum_impurities GParam.h5/mott/impurity_{i}/num_mott_orbitals
    # ne_mott = sum_impurities GParam.h5/mott/impurity_{i}/num_mott_electrons
    ferwes = get_fermi_weight(efermi, bnd_es, wklist, delta=delta,
            ismear=ismear, iso=iso, no_mott=no_mott, ne_mott=ne_mott)
    
    # get_wannier_den_matrix_risb calculates and returns the density matrix
    # There is special logic if ispin_risb > ispin_dft.
    #   There is also a check on whether the total number of electrons in the density
    #       matrix matches the total number calculated from ferwes.
    #  density matrix = \sum_{spin,k,bands} (1/wk * f_ispin) * (r_mat^\dagger * bnd_vs * ferwek * bndvs^\dagger * r_mat)
    #       + (n_phy_mat-nr_mat)T  
    # \rho_{A,B} =  R^\dagger_{A,a} * <a|psi>f<psi|b> * R_{b,B} 
    #            + n_phys.^{A,B} - n_{sub.}^{A,B}
    #   wk is the k-point weights
    #   f_ispin = 1 if doing spin-orbit, 2 otherwise
    #   bnd_vs is the eigenvectors of the Gutzwiller hamiltonian
    #   ferwek is the Fermi occupation factors, computed using the results of the Gutzwiller calculation
    #   nktot is the number of k-points
    #   r_mat is the quasiparticle renormalization matrix, from GLog.h5
    #   n_phys is the physical density matrix to be added in density calculations, from GLog.h5
    #   nr_mat is the local density part to be subtracted in density calculations, from GLog.h5
    wan_den = get_wannier_den_matrix_risb(bnd_vs, ferwes, wklist, nktot)
    
    # double check total charge
    # chk_charge(wan_den, wklist)
    
    wfwannier_list = comm.gather(wfwannier_list, root=0)
    
    if rank == 0:
        wfwannier_list = np.concatenate(wfwannier_list, axis=1)
        
        # wrt_wan_den writes out wannier_den_matrix.dat
        # If file fwannier is present on disk and is an hdf5 file, then it writes out an hdf5 file.
        # Otherwise it writes out a binary file.
        # wannier_den_matrix contains the following info:
        #   'temperature' = 300.0
        #   'nqdiv' = the number of k-points, obtained from the shape of wan_den
        #   'num_band_in' = the number of bands, obtained from the shape of wfwannier_list
        #   'n_dmft_wan' = the number of wanniers, obtained from the shape of wfwannier_list
        #   'include_band_new' = include_bands
        #   'vw_matrix_new' = wf_wannier_list[0]
        #   'n_matrix_new' = wan_den[0], with the axes swapped.
        wrt_wan_den(wan_den, wfwannier_list, include_bands,
                fwannier="{}/wannier.dat".format(wannier_path))

# never used
def chk_charge(wan_den, wklist):
    nelec = 0.
    for isp, den in enumerate(wan_den):
        for den1, wk1 in zip(den, wklist):
            nelec += np.trace(den1)*wk1
    print(f"total electron in chk_charge: {nelec:.6f}")

# never used
def wannier_den_matrix_lda_chk(wannier_path="./"):
    '''produce the file `wannier_den_matrix.dat` for the feedback from
    g-risb to dft.
    '''
    _, _, kpts, include_bands, _, bnd_es = mpiget_wannier_data(
            path=wannier_path)
    nktot = len(kpts)
    with h5py.File("GBareH.h5", "r") as f:
        num_elec = f["/"].attrs["nelectron"]

    # chop bnd_es
    nbnd = bnd_es.shape[2]
    nwan = int(num_elec/2+3)
    bnd_es = bnd_es[:,:,:nwan]
    # set wk_list
    wklist = [1./nktot for i in range(bnd_es.shape[1])]
    efermi = get_fermi_level(bnd_es, wklist, num_elec, ismear=-1)
    # set fermi weight
    ferwes = get_fermi_weight(efermi, bnd_es, wklist, ismear=-1)
    # setup trivial wannier_den data.
    wan_den = [[]]
    wfwannier_list = [[]]
    vmat = np.zeros((nbnd, nwan), dtype=np.complex)
    np.fill_diagonal(vmat, 1.0)
    for ik in range(nktot):
        wfwannier_list[-1].append(vmat)
        wan_den[0].append(np.diag(ferwes[0][ik]*nktot/(2+0.j)))
    wan_den = np.asarray(wan_den)
    wfwannier_list = np.asarray(wfwannier_list)
    wrt_wan_den(wan_den, wfwannier_list, include_bands,
            fwannier="{}/wannier.dat".format(wannier_path))

# never used
def wannier_den_matrix_lda_chk2(wannier_path="./"):
    '''produce the file `wannier_den_matrix.dat` for the feedback from
    g-risb to dft.
    '''
    _, _, kpts, include_bands, wfwannier_list, bnd_es = mpiget_wannier_data(
            path=wannier_path)
    nktot = len(kpts)
    with h5py.File("GBareH.h5", "r") as f:
        num_elec = f["/"].attrs["nelectron"]
    # set wk_list
    wklist = [1./nktot for i in range(bnd_es.shape[1])]
    efermi = get_fermi_level(bnd_es, wklist, num_elec, ismear=-1)
    # set fermi weight
    ferwes = get_fermi_weight(efermi, bnd_es, wklist, ismear=-1)
    # setup trivial wannier_den data.
    wan_den = [[]]
    for ik in range(nktot):
        dm = np.diag(ferwes[0][ik]*nktot/(2+0.j))
        wan_den[0].append(wfwannier_list[0][ik].T.conj().dot(dm).dot(
            wfwannier_list[0][ik]))
    wan_den = np.asarray(wan_den)
    wrt_wan_den(wan_den, wfwannier_list, include_bands,
            fwannier="{}/wannier.dat".format(wannier_path))

# never used
def wannier_den_matrix_lda_chk3(wannier_path="./"):
    '''produce the file `wannier_den_matrix.dat` for the feedback from
    g-risb to dft.
    '''
    _, _, kpts, include_bands, wfwannier_list, bnd_es = mpiget_wannier_data(
            path=wannier_path)
    nktot = len(kpts)
    with h5py.File("GBareH.h5", "r") as f:
        num_elec = f["/"].attrs["nelectron"]
    # set wk_list
    wklist = [1./nktot for i in range(bnd_es.shape[1])]
    # get eigen-vector from interpolation
    bnd_es2 = [[]]
    bnd_ev2 = [[]]
    for ik in range(nktot):
        hk = wfwannier_list[0][ik].T.conj().dot(np.diag(bnd_es[0][ik])).dot(
                wfwannier_list[0][ik])
        w, v = np.linalg.eigh(hk)
        bnd_es2[0].append(w)
        bnd_ev2[0].append(v)
    bnd_ev2 = np.asarray(bnd_ev2)
    with h5py.File("ev_lda_ref.h5", "w") as f:
            f["e"] = bnd_es2
            f["v"] = bnd_ev2
    efermi = get_fermi_level(bnd_es2, wklist, num_elec, ismear=-1)
    # set fermi weight
    ferwes = get_fermi_weight(efermi, bnd_es2, wklist, ismear=-1)
    # setup trivial wannier_den data.
    wan_den = [[]]
    for ik in range(nktot):
        dm = np.diag(ferwes[0][ik]*nktot/(2+0.j))
        wan_den[0].append(bnd_ev2[0][ik].dot(dm).dot(
            bnd_ev2[0][ik].T.conj()))
    wan_den = np.asarray(wan_den)
    
    wrt_wan_den(wan_den, wfwannier_list, include_bands,
            fwannier="{}/wannier.dat".format(wannier_path))


# wrt_wan_den writes out wannier_den_matrix.dat
# If file fwannier is present on disk and is an hdf5 file, then it writes out an hdf5 file.
# Otherwise it writes out a binary file.
# wannier_den_matrix contains the following info:
#   'temperature' = 300.0
#   'nqdiv' = the number of k-points, obtained from the shape of wan_den
#   'num_band_in' = the number of bands, obtained from the shape of wfwannier_list
#   'n_dmft_wan' = the number of wanniers, obtained from the shape of wfwannier_list
#   'include_band_new' = include_bands
#   'vw_matrix_new' = wf_wannier_list[0]
#   'n_matrix_new' = wan_den[0], with the axes swapped.
def wrt_wan_den(wan_den, wfwannier_list, include_bands, fwannier=None):
    try:
        # check hdf5 version.
        f = h5py.File(fwannier, "r")
        f.close()
        h5wrt_wan_den(wan_den, wfwannier_list, include_bands)
    except:
        fwrt_wan_den(wan_den, wfwannier_list, include_bands)


def fwrt_wan_den(wan_den, wfwannier_list, include_bands,
        fname='wannier_den_matrix.dat',
        ):
    """write wannier_den_matrix.dat file in fortran binary format.
    for spin-restricted case only as of now.
    note the index order wfwannier_list[isp, ibnd, iwann].
    not reliable, integer can be read wrong!
    """
    nktot = wan_den.shape[1]
    wfwannier_list = wfwannier_list.swapaxes(2, 3)
    nband = wfwannier_list.shape[3]
    nwann = wfwannier_list.shape[2]
    with FortranFile(fname, 'w') as f:
        f.write_record(300.0)
        f.write_record(nktot)
        f.write_record([nband for k in range(nktot)])
        f.write_record([nwann for k in range(nktot)])
        f.write_record([include_bands for k in range(nktot)])
        f.write_record(wfwannier_list[0])
        f.write_record(wan_den[0].swapaxes(1, 2))

# writes out wannier_den_matrix.dat
# If file fwannier is present on disk, then it writes out an hdf5 file.
# Otherwise it writes out a binary file.
# wannier_den_matrix contains the following info:
#   'temperature' = 300.0
#   'nqdiv' = the number of k-points, obtained from the shape of wan_den
#   'num_band_in' = the number of bands, obtained from the shape of wfwannier_list
#   'n_dmft_wan' = the number of wanniers, obtained from the shape of wfwannier_list
#   'include_band_new' = include_bands
#   'vw_matrix_new' = wf_wannier_list[0]
#   'n_matrix_new' = wan_den[0], with the axes swapped.
def h5wrt_wan_den(wan_den, wfwannier_list, include_bands,
        fname='wannier_den_matrix.dat',
        ):
    """write wannier_den_matrix.dat file in hdf5 format.
    for spin-restricted case only as of now.
    note the index order wfwannier_list[isp, ibnd, iwann].
    """
    nktot = wan_den.shape[1]
    wfwannier_list = wfwannier_list.swapaxes(2, 3)
    nband = wfwannier_list.shape[3]
    nwann = wfwannier_list.shape[2]
    with h5py.File(fname, 'w') as f:
        f["/temperature"] = [300.0]
        f["/nqdiv"] = [nktot]
        f["/num_band_in"] = [nband for k in range(nktot)]
        f["/n_dmft_wan"] = [nwann for k in range(nktot)]
        f["/include_band_new"] = [include_bands for k in range(nktot)]
        f["/vw_matrix_new"] = wfwannier_list[0]
        f["/n_matrix_new"] = wan_den[0].swapaxes(1, 2)
