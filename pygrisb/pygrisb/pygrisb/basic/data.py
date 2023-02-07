import h5py, numpy


def compare_array(path,
        fname1='GLOG.h5',
        fname2='GLOG_REF.h5',
        chk_sum=False,
        ):
    with h5py.File(fname1, 'r') as f:
        data1 = f[path][()]
    with h5py.File(fname2, 'r') as f:
        data2 = f[path][()]
    if chk_sum:
        return abs(numpy.sum(data1) - numpy.sum(data2)) < 1.e-5
    else:
        return numpy.allclose(data1, data2)

