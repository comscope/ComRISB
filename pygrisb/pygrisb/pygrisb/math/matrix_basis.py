import numpy as np
from abc import ABC, abstractmethod


class matrix_basis(ABC):
    """
    base class of complete matrix basis set.
    """
    def __init__(self, symbol_matrix=None, basis=None, dtype=np.complex, \
            btype="hermi"):
        self._symbol_matrix = symbol_matrix
        self._basis = basis
        self.dtype = dtype
        self.btype = btype # general or hermitian matrix basis
        self.set_basis_set()

    @property
    def basis(self):
        return self._basis

    def expansion_coefficient(self, a):
        pass

    @abstractmethod
    def set_basis_set(self):
        pass


class dense_matrix_basis(matrix_basis):
    def set_basis_set(self):
        if self._basis is None:
            matrix_basis = []
            symbol_matrix = np.asarray(self._symbol_matrix, dtype=np.int)
            for element in range(1,symbol_matrix.max()+1):
                spots = np.argwhere(symbol_matrix == element)
                num_spots = len(spots)
                if num_spots == 0:
                    continue
                # for easier access
                spots = spots.T
                zero_matrix = np.zeros_like(symbol_matrix, dtype=self.dtype)
                # hermitian or symmetric matrix
                if "h" == self.btype[0].lower():
                    # Skip if located at lower trigonal block
                    if spots[0][0] > spots[1][0]:
                        continue
                    if spots[0][0] == spots[1][0]:
                        value = 1./np.sqrt(float(num_spots))
                        matrix = zero_matrix.copy()
                        matrix[spots[0], spots[1]] = value
                        matrix_basis.append(matrix)
                    else:
                        # non-zero element
                        value = 1./np.sqrt(float(num_spots*2))
                        matrix = zero_matrix.copy()
                        matrix[spots[0], spots[1]] = value
                        matrix[spots[1], spots[0]] = value
                        matrix_basis.append(matrix)

                        if self.dtype in [np.complex, np.complex_]:
                            value = value * 1.j
                            matrix = np.zeros_like(symbol_matrix,
                                    dtype=np.complex)
                            matrix[spots[0], spots[1]] = value
                            matrix[spots[1], spots[0]] = -value
                            matrix_basis.append(matrix)
                # general matrix
                else:
                    value = 1/np.sqrt(float(num_spots))
                    matrix = zero_matrix.copy()
                    matrix[spots[0], spots[1]] = value
                    matrix_basis.append(matrix)
            self._basis = matrix_basis
