import sys
import numpy as np
import scipy.sparse as sps
import scipy.linalg as slinalg
from itertools import product


class matrix:
    '''base matrix class to represent various matrix forms.
    '''
    def __init__(self):
        self.is_sm = False
        self.is_bd = False
        self.evaluate()

    @property
    def a(self):
        return self._a

    @property
    def is_sparse_matrix(self):
        return self.is_sm

    @property
    def is_block_diagonal(self):
        return self.is_bd

    def evaluate(self):
        raise NotImplementedError("function evaluate not implemented!")


class general_matrix(matrix):
    '''general class for dense and sparse matrix
    with possible block_diagonal form.
    '''
    def __init__(self, a, bkdim_list=None):
        self._a = a
        self.bkdim_list = bkdim_list
        super().__init__()

    def evaluate(self):
        self.is_bd = True if self.bkdim_list is not None else False
        if sps.issparse(self._a):
            self.is_sm = True
        elif type(self._a) is np.ndarray:
            self.is_sm = False
        else:
            raise TypeError("input matrix not sparse or dense!")


class eigen_system:
    '''base class to handle eigen system of a general matrix.
    '''
    def __init__(self, mat, **kwargs):
        self.matrix = mat
        self.kwargs = kwargs
        self.evaluate()

    @property
    def eigen_values(self):
        return self.evals

    @property
    def eigen_vectors(self):
        return self.evecs

    def evaluate(self):
        if self.matrix.is_bd:
            val_list = []
            vec_list = []
            nbase = 0
            for n in self.matrix.bkdim_list:
                submat = matrix(self.matrix.a[nbase:nbase+n, nbase:nbase+n])
                esys = eigen_system(submat, **self.kwargs)
                vals = esys.eigen_values
                vecs = esys.eigen_vectors
                nbase += n
                val_list += vals.tolist()
                vec_list.append(vecs)
            self.evals = np.array(val_list)
            self.vecs = slinalg.block_diag(*vec_list)
        elif self.matrix.is_sm:
            self.evals, self.evecs = sps.linalg.eigsh(self.matrix.a,
                    **self.kwargs)
        else:
            self.evals, self.evecs = np.linalg.eigh(self.matrix.a,
                    **self.kwargs)


class eigen_space:
    '''class to calculate eigen-spaces.
    '''
    def __init__(self, evals, evecs=None, rtol=1.e-6):
        self.evals = evals
        self.evecs = evecs
        self.rtol = rtol
        self.eval_list = None
        self.evecs_list = []
        self.evaluate()

    @property
    def eigen_space_values(self):
        return self.eval_list

    @property
    def eigen_space_vectors(self):
        return self.evecs_list

    @property
    def idx(self):
        return self.es_idx

    @classmethod
    def from_matrix(cls, mat, rtol=1.e-6, **kwargs):
        esys = eigen_system(mat, **kwargs)
        return cls(esys.eigen_values, esys.eigen_vectors, rtol=rtol)

    def print_result(self, log=sys.stdout):
        print(f'Number of eigen spaces = {len(self.eval_list)}', file=log)
        print('Dimension of each eigen space:', file=log)
        for i, evecs in enumerate(self.evecs_list):
            print(i, evecs.shape[1], file=log)

    def evaluate(self):
        self.set_eigen_space_indices()
        es_idx = self.es_idx
        self.eval_list = self.evals[es_idx[:-1]]
        if self.evecs is not None:
            for i,j in zip(es_idx[:-1], es_idx[1:]):
                self.evecs_list.append(self.evecs[:,i:j])

    def set_eigen_space_indices(self):
        # evals should have been ordered.
        # get the pairwise difference.
        ediff = self.evals[1:] - self.evals[:-1]
        # take care of the beginning and end
        beg_end = np.ones(2)+self.rtol
        ediff = np.insert(beg_end, -1, ediff)
        es_idx = np.argwhere(np.abs(ediff) > self.rtol)
        self.es_idx = es_idx.reshape(-1)


    def get_reduce_labels(self, labels, method="average"):
        '''reduce labels for eigen-space. method = 'average' ot 'sum'.
        '''
        if labels is None:
            return None
        s_labels = []
        for i,j in zip(self.es_idx[:-1], self.es_idx[1:]):
            if j<=i:
                continue
            tmp = np.sum(labels[i:j])
            if method == "average":
                tmp /= j-i
            s_labels.append(tmp)
        return np.array(s_labels)


class random_hermitian_matrix(matrix):
    '''class of random hermitian matrix(n, n).
    '''
    def __init__(self, n):
        self.n = n
        super().__init__()

    def evaluate(self):
        x = np.random.rand(self.n, self.n) + \
                1.j*np.random.rand(self.n, self.n)
        self._a = x + x.conj().T


class sym_matrix(matrix):
    '''base class of symmetrized matrix
    '''
    def __init__(self, a, u_list):
        self._a = a
        self.u_list = u_list
        self.order = len(u_list)
        super().__init__()


class sym_dense_matrix(sym_matrix):
    '''class of symmetrized dense matrix
    '''
    def evaluate(self):
        a_sym = np.zeros_like(self._a)
        for u in self.u_list:
            a_sym += u.conj().T.dot(self._a).dot(u)
        self._a = a_sym/self.order


class sym_rnd_hermi_matrix(sym_dense_matrix):
    '''class of random hermitian (dense) matrix, symmetrized
    with respect to a group represented by u_list.
    '''
    def __init__(self, u_list):
        dim = u_list[0].shape[0]
        a = random_hermitian_matrix(dim).a
        super().__init__(a, u_list)


class sym_sparse_matrix(sym_matrix):
    '''class of symmetrized dense matrix
    '''
    def evaluate(self):
        a_sym = self._a.__class__(self._a.shape, dtype=self._a.dtype)
        for u in self.u_list:
            a_sym += u.getH()*self._a*u
        self._a = a_sym/self.order


class proj_expan_dense_matrix(sym_matrix):
    def evaluate(self):
        a_sym = np.zeros_like(self._a)
        for u in self.u_list:
            a_sym += np.vdot(u, self._a) * u
        self._a = a_sym


def get_loewner_matrix_ev(e_vals, f, fp):
    '''
    Loewner matrix, ref: positive definite matrices, Bhatia, p.154.
    f: function; fp: first derivative.
    '''
    loewner = []
    for i, ei in enumerate(e_vals):
        _loewner = []
        for j, ej in enumerate(e_vals):
            if np.abs(ei - ej) < 1.e-10:
                res = fp(ei)
            else:
                res = (f(ei) - f(ej))/(ei - ej)
            _loewner.append(res)
        loewner.append(_loewner)
    return np.array(loewner)


def yield_derivative_f_matrix(a, h_list, f, fp):
    '''
    a = \sum_{c} {c*h_list}. fp(x) = \par f(x) / \par x
    Yields \par f(a) / \par c.
    '''
    w, v = np.linalg.eigh(a)
    vherm = np.conj(v.T)
    loewner = get_loewner_matrix_ev(w, f, fp)
    for h in h_list:
        _h = np.dot(vherm, np.dot(h, v))
        yield np.dot(v, np.dot(loewner*_h, vherm))


def unitary_transform_coulomb_matrix(a, u):
    '''Perform a unitary transformation (u) on the Coulomb matrix (a).
    '''
    a_ = np.asarray(a).copy()
    m_range = range(a.shape[0])
    for i,j in product(m_range, repeat=2):
        a_[i,j,:,:] = u.T.conj().dot(a_[i,j,:,:].dot(u))
    a_ = a_.swapaxes(0,2).swapaxes(1,3)
    for i,j in product(m_range, repeat=2):
        a_[i,j,:,:] = u.T.conj().dot(a_[i,j,:,:].dot(u))
    return a_


def get_hamilt_matrix_from_ev(w, v):
    '''given v_{alpha, n} and w_{n}, get back the Hamiltonian matrix.
    Here alpha is the basis orbital index and n the band index.
    '''
    vh = v.T.conj()
    for i in range(v.shape[0]):
        v[i] *= w
    return v.dot(vh)


def solve_ab_plus_ba_equal_c(b, c):
    '''solve for a given a*b + b*a = c and b is hermitian.
    '''
    err = np.max(np.abs(b-b.T.conj()))
    if err > 1.e-8:
        raise AssertionError(f"b not hermitian with error = {err:.2e}!")
    w,v = np.linalg.eigh(b)
    cp = v.T.conj().dot(c).dot(v)
    ap = [[cp[i,j]/(w[i]+w[j]) for j in range(b.shape[0])] \
            for i in range(b.shape[0])]
    a = v.dot(ap).dot(v.T.conj())
    return a


def get_1stderivative_sroot_entropyn_at_h(a, h):
    '''
    evaluate p_{\squareroot(a(1-a))} / p_h.
    '''
    k0 = a.dot(np.eye(a.shape[0])-a)
    x0 = slinalg.sqrtm(k0)
    k1 = h - a.dot(h) - h.dot(a)
    res = solve_ab_plus_ba_equal_c(x0, k1)
    return res


def get_2ndderivative_sroot_entropyn_at_h(a, h1, h2):
    '''
    evaluate p2_{\squareroot(a(1-a))} / p_h1 / p_h2.
    '''
    p1 = get_1stderivative_sroot_entropyn_at_h(a, h1)
    p2 = get_1stderivative_sroot_entropyn_at_h(a, h2)
    k2 = -h1.dot(h2) - h2.dot(h1)
    c = k2 - (p1.dot(p2) + p2.dot(p1))
    k0 = a.dot(np.eye(a.shape[0])-a)
    x0 = slinalg.sqrtm(k0)
    res = solve_ab_plus_ba_equal_c(x0, c)
    return res


def list_mat_to_array(a_list):
    '''convert a list of matrix to 3d array.
    '''
    dim_list = np.array([len(a) for a in a_list])
    if np.all(np.abs(dim_list-dim_list[0])<.1):
        a_array = np.asarray(a_list)
    else:
        nmax = np.max(dim_list)
        a_array = np.zeros((len(a_list), nmax, nmax),
                dtype=a_list[0][0][0].dtype)
        for i, a in enumerate(a_list):
            n = len(a)
            a_array[i,:n,:n] = a
    return a_array


def get_block_diag_indices(a, tol=1.e-8):
    '''get the block diagonal indices.
    '''
    ind = [0]
    for i in range(a.shape[-1]-1):
        if np.all(np.abs(a[:i+1, i+1:]) < tol) and \
                np.all(np.abs(a[i+1:, :i+1]) < tol):
            ind.append(i+1)
    ind.append(a.shape[-1])
    return ind
