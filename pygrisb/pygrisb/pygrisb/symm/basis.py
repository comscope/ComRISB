import numpy as np
from scipy.linalg import block_diag
import pygrisb.symm.unitary as un
from pygrisb.symm.group import group_decomp


class symm_basis:
    '''base class for symmetry-adapted basis set.
    '''
    def __init__(self):
        self.evaluate()

    @property
    def transform_matrix(self):
        '''return the symmetry adapated unitary transformation.
        '''
        return self.u_trans

    @property
    def symbolic_matrix(self):
        '''return a symbolic matrix representing a general matrix
        commuting with the representation of the group.
        '''
        return self.symbol_matrix

    def display_symbolic_matrix(self):
        print("symbol matrix:")
        emax = np.max(self.symbol_matrix)
        n = len(str(emax))
        for row in self.symbol_matrix:
            print('  '+' '.join([f"{e:{n}d}" for e in row]))

    def spin_blk_swap(self):
        u = un.spin_blk_swap(self.u_trans.shape[0]).u
        self.u_trans = u.dot(self.u_trans)

    def evaluate(self):
        raise NotImplementedError("function evaluate not implemented!")


class symm_basis_pg(symm_basis):
    '''class to get symmetry-adapted basis
    and the general matrix structure commuting with point group
    symmetry operations.
    '''
    def __init__(self, u_list):
        self.u_list = u_list
        # basis dimension
        self.dim = u_list[0].shape[0]
        super().__init__()

    def evaluate(self):
        self.gd = group_decomp(self.u_list)
        self.beautify()
        self.set_symbolic_matrix()

    def beautify(self):
        # get first piece
        self.u_trans = np.empty((self.dim, 0))
        for equ_ireps, equ_evecs in zip(self.gd.equ_ireps_list,
                self.gd.equ_evecs_list):
            for i, irep in enumerate(equ_ireps[1:]):
                ireps_u = un.ireps_unitary_trans(irep, equ_ireps[0])
                u = ireps_u.u
                equ_evecs[i] = equ_evecs[i].dot(u.T.conj())
            self.u_trans = np.hstack((self.u_trans, \
                    np.moveaxis(equ_evecs, 0, -1).reshape((self.dim, -1))))

    def set_symbolic_matrix(self):
        sm = []
        ibase = 1
        for equ_ireps in self.gd.equ_ireps_list:
            # number of equivalent irred. reps.
            n = len(equ_ireps)
            # ired space dimension
            d = equ_ireps[0][0].shape[0]
            sm1 = np.arange(n*n).reshape((n,n)) + ibase
            sm += [sm1]*d
            ibase += n*n
        self.symbol_matrix = block_diag(*sm)


class symm_basis_uniform(symm_basis):
    def __init__(self, n, dtype=np.complex):
        self.n = n
        self.dtype = dtype
        super().__init__()

    def evaluate(self):
        self.symbol_matrix = np.eye(self.n, dtype=np.int)
        self.u_trans = np.eye(self.n, dtype=self.dtype)


class symm_basis_random(symm_basis):
    def __init__(self, n, dtype=np.complex):
        self.n = n
        self.dtype = dtype
        super().__init__()

    def evaluate(self):
        self.symbol_matrix = np.arange(1, self.n**2+1).reshape(self.n, self.n)
        self.u_trans = np.eye(self.n, dtype=self.dtype)


class symm_basis_soc(symm_basis):
    def __init__(self, l, ispin, dtype=np.complex):
        self.l = l
        self.ispin = ispin
        self.dtype = dtype
        self.dim = (l*2+1)*2
        super().__init__()

    def evaluate(self):
        self.set_symbolic_matrix()
        self.u_trans = np.asarray(un.comp_sph_harm_to_relat_harm(self.l).u,
                dtype=self.dtype)

    def set_symbolic_matrix(self):
        '''setup symbolic matrix fir the case with spin-orbit interaction
        and no crystal field effect.
        '''
        elem_base = 1
        l = self.l
        if l == 0:
            if self.ispin == 1:
                elem = [elem_base, elem_base]
            else:
                elem = [elem_base, elem_base + 1]
        else:
            elem = []
            if self.ispin == 1:
                elem += [elem_base for i in \
                        range(int(2 * (l - 0.49)) + 1)]
                elem_base = max(elem) + 1
                elem += [elem_base for i in \
                        range(int(2 * (l + 0.51)) + 1)]
            else:
                elem += [elem_base + i for i in \
                        range(int(2 * (l - 0.49)) + 1)]
                elem_base = max(elem) + 1
                elem += [elem_base + i for i in \
                        range(int(2 * (l + 0.51)) + 1)]
        self.symbol_matrix = np.diag(elem)

