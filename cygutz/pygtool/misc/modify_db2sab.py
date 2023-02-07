import h5py, numpy
from pygrisb.symm.unitary import complx_sph_harm_to_real_harm
'''
Recipe to switch to real-embedding Hamiltonian calculation.
1. modify db2sab
2. init_ga.py --fixsab --realhemb -u ev
'''


csh2rsh = complx_sph_harm_to_real_harm(l=2).u_trans
csh2t2geg = csh2rsh.copy()
csh2t2geg[:, 2] = csh2rsh[:, 3]
csh2t2geg[:, 3] = csh2rsh[:, 2]
csh2t2geg2 = numpy.zeros([10, 10], dtype=numpy.complex)
csh2t2geg2[:5, ::2] = csh2t2geg2[5:, 1::2] = csh2t2geg


with h5py.File("GParam.h5", "a") as f:
    f["/impurity_0/db2sab"][()] = csh2t2geg2
