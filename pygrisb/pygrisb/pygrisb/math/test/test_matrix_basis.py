import numpy
import pygrisb.math.matrix_basis as mb


m_basis = mb.dense_matrix_basis(symbol_matrix=numpy.arange(1,5).reshape(2,2))
for a in m_basis.basis:
    print(a)
