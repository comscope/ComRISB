import numpy, h5py
from mpi4py import MPI
from scipy.io import FortranFile

# get_wannier_data_binary reads info out of wannier.dat
#   - reals_lat is the real space lattice vectors
#   - recip_lat is the reciprocal space lattice vectors
#   - kpts is a list of all the k-points
#   - include_bands is low energy band indices included in the wannier construction
#   - wfwannier_list is a list of overlap between band wavefunctions and wannier orbitals
#   - bnd_es list of band energies
# Things read but not returned directly:
#   - num_bands is the maximal number of bands (controls the size of wfwannier_list,
#            bnd_es)
#   - num_wann is the number of wannier orbitals ( controlls the size of wfwannier_list)
#   - ndiv is the k-point mesh (not returned at all, but the number of k-points
#           controls  the size of kpts, wfwannier_list, bnd_es)
def get_wannier_data_binary(fname):
    with FortranFile(fname, "r") as f:
        # real space lattice vectors
        reals_lat = f.read_reals().reshape((3, 3))
        # reciprocal space lattice vectors
        recip_lat = f.read_reals().reshape((3, 3))
        # maximal number of bands
        num_bands = f.read_ints()[0]
        # number of wannier orbitals
        num_wann = f.read_ints()[0]
        # k-point mesh
        ndiv = f.read_ints()
        # total number of k-points
        nqdiv = ndiv[0]*ndiv[1]*ndiv[2]
        # k-points
        kpts = f.read_reals().reshape((nqdiv, 3))
        # low energy band indices included in the wannier construction
        include_bands = f.read_ints()
        # list of overlap between band wavefunctions and wannier orbitals
        wfwannier_list = f.read_reals().view(numpy.complex).reshape(\
                (1, nqdiv, num_wann, num_bands)).swapaxes(2, 3)
        # list of band energies
        bnd_es = f.read_reals().reshape((1, nqdiv, num_bands))
    return reals_lat, recip_lat, kpts, include_bands, wfwannier_list, bnd_es


def h5get_wannier_data(fname):
    with h5py.File(fname, "r") as f:
        # real space lattice vectors
        reals_lat = f["/rbas"][()]
        # reciprocal space lattice vectors
        recip_lat = f["/gbas"][()]
        # maximal number of bands
        num_bands = f["/num_bands"][0]
        # number of wannier orbitals
        num_wann = f["/num_wann"][0]
        # k-point mesh
        ndiv = f["/ndiv"][()]
        # total number of k-points
        nqdiv = ndiv[0]*ndiv[1]*ndiv[2]
        # k-points
        kpts = f["/kpt_latt"][()]
        # low energy band indices included in the wannier construction
        include_bands = f["/include_bands"][()]
        # list of overlap between band wavefunctions and wannier orbitals
        wfwannier_list = f["/v_matrix"][()].reshape(\
                (1, nqdiv, num_wann, num_bands)).swapaxes(2, 3)
        # list of band energies
        bnd_es = f["/eigenvalues"][()].reshape((1, nqdiv, num_bands))
    return reals_lat, recip_lat, kpts, include_bands, wfwannier_list, bnd_es

# get_wannier_data calls get_wannier_data_binary, which reads data
#    out of wannier.dat . For more info look at the comments on get_wannier_data_binary.
def get_wannier_data(path="./"):
    fname = "{}/wannier.dat".format(path)
    try:
        f = h5py.File(fname, "r")
        f.close()
        h5input = True
    except:
        h5input = False
        pass
    if h5input:
        reals_lat, recip_lat, kpts, include_bands, wfwannier_list, bnd_es = \
                h5get_wannier_data(fname)
    else:
        reals_lat, recip_lat, kpts, include_bands, wfwannier_list, bnd_es = \
                get_wannier_data_binary(fname)
    # check orthonormality
    # chk_orthonormality(wfwannier_list)
    return reals_lat, recip_lat, kpts, include_bands, wfwannier_list, bnd_es


def chk_orthonormality(vw_list):
    imat = numpy.eye(vw_list.shape[3])
    print(f"vw_list shape: {vw_list.shape}")
    for vw in vw_list:
        for vwk in vw:
            a = vwk.T.conj().dot(vwk)
            if not numpy.allclose(a, imat):
                print(a)
                raise ValueError("vw not orthonormal.")


# mpiget_wannier_data calls get_wannier_data_binary, which reads data
#    out of wannier.dat . For more info look at the comments on get_wannier_data_binary.
#   It broadcasts reals_lat, recip_lat, kpts, bnd_es, include_bands to all mpi processes.
#   It also splits kpts and wfwannier_list between mpi processes.
#       After the split, the local process' k-points are stored in k_range,
#       and its share of wfwannier_list is still called wfwannier_list.
def mpiget_wannier_data(path="./"):
    '''get the contents in wannier.dat Fortran binary file.
    '''
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    ncpu = comm.Get_size()
    
    if rank == 0:
        reals_lat, recip_lat, kpts, include_bands, wfwannier_list, bnd_es =\
                get_wannier_data(path=path)
    else:
        reals_lat = recip_lat = kpts = include_bands = wfwannier_list = \
                bnd_es = k_range = None
                
    reals_lat = comm.bcast(reals_lat, root=0)
    recip_lat = comm.bcast(recip_lat, root=0)
    kpts = comm.bcast(kpts, root=0)
    bnd_es = comm.bcast(bnd_es, root=0)
    include_bands = comm.bcast(include_bands, root=0)
    
    if rank == 0:
        from pygrisb.mpi.mpi import get_k_range_list as get_k_range_list
        k_range_list = get_k_range_list(kpts.shape[0])
        # change wfwannier_list to blocks
        wfwannier_bklist = [wfwannier_list[:, \
                k_range_list[i][0]:k_range_list[i][1],:,:] \
                for i in range(ncpu)]
    else:
        wfwannier_bklist = None
        k_range_list = None

    # distributed evenly among cores with respect to k-points
    wfwannier_list = comm.scatter(wfwannier_bklist, root=0)
    k_range = comm.scatter(k_range_list, root=0)
    
    return reals_lat, recip_lat, kpts, include_bands, wfwannier_list, bnd_es, \
            k_range


if __name__ == "__main__":
    get_wannier_data()
