from scipy.sparse import csr_matrix,coo_matrix
from time import sleep
import h5py



def h5auto_read(f, path, default=None):
    if path in f:
        return f[path][()]
    else:
        return default


def h5auto_write(f, path, data):
    if path in f:
        del f[path]
    f[path] = data


def get_csr_matrix(f, path):
    '''
    Read the csr_matrix located at path in the hdf5 file f.
    '''
    nrow = f[path].attrs["nrow"]
    ncol = f[path].attrs["ncol"]
    data = f[path + "/data"][()]
    base = f[path].attrs["base"]
    # possible one-based to zero-based
    indices = f[path + "/indices"][()] - base
    indptr = f[path + "/indptr"][()] - base
    return csr_matrix((data, indices, indptr), shape=(nrow, ncol))


def get_coo_matrix(f, path):
    '''
    Read the coo_matrix located at path in the hdf5 file f.
    '''
    nrow = f[path].attrs["nrow"]
    ncol = f[path].attrs["ncol"]
    data = f[path + "/data"][()]
    base = f[path].attrs["base"]
    # possible one-based to zero-based
    indi = f[path + "/i"][()] - base
    indj = f[path + "/j"][()] - base
    return coo_matrix((data, (indi, indj)), shape=(nrow, ncol))


def write_csr_matrix(f, path, a):
    '''
    Read the csr_matrix located at path in the hdf5 file f.
    '''
    if path in f:
        del path
    f[path].attrs["nrow"] = a.shape[0]
    f[path].attrs["ncol"] = a.shape[1]
    f[path + "/data"] = a.data
    f[path].attrs["base"] = 0
    f[path + "/indices"] = a.indices
    f[path + "/indptr"] = a.indptr


def write_coo_matrix(f, path, a):
    '''
    Write the coo_matrix located at path in the hdf5 file f.
    '''
    f[path].attrs["nrow"] = a.shape[0]
    f[path].attrs["ncol"] = a.shape[1]
    f[path].attrs["/nnz"] = a.nnz
    f[path + "/data"] = a.data
    f[path].attrs["base"] = 0
    f[path + "/i"] = a.row
    f[path + "/j"] = a.col


def get_hs_rotations(f, imp, valences):
    '''
    Get rotation representations in Hilbert space.
    '''
    Rpr_list = []
    for val in valences:
        Rpr_list.append([])
        dim_rot = f["Impurity_{}/val_block={}/dim_rotations".format( \
                imp, val)][()]
        for i in range(dim_rot):
            Rpr_list[-1].append(get_csr_matrix(f, \
                    "/Impurity_{}/val_block={}/rotation_{}".format( \
                    imp, val, i)))
    return Rpr_list


def h5open(filename, *args, **kwargs):
    while True:
        try:
            hdf5_file = h5py.File(filename, *args, **kwargs)
            break  # success!
        except OSError:
            print(f"{filename} already opened, wait 5s!")
            sleep(5)  # wait a bit

    return hdf5_file
