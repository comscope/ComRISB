# with parts from U_matrix.py in TRIQS package.
from itertools import product
import numpy as np
import sys
from scipy.special import factorial as fact


class coulomb_matrix:
    r'''base class for coulomb matrix. chemist convention adopted, namely:

    .. math:: H = \frac{1}{2} \sum_{ijkl,\sigma \sigma'} V_{ikjl}
            a_{i \sigma}^\dagger a_{j \sigma'}^\dagger
            a_{l \sigma'} a_{k \sigma}.
    '''
    def __init__(self):
        self.evaluate_all()

    @property
    def u_orb(self):
        return self._v

    @property
    def u_spin_orb(self):
        return self._v2

    @property
    def u_avg(self):
        return self._u_avg

    @property
    def j_avg(self):
        return self._j_avg

    def evaluate_all(self):
        self.evaluate()
        self.set_uj_avg()
        self.set_full_u_matrix()

    def evaluate(self):
        raise NotImplementedError("evaluate function not implemented.")

    def set_uj_avg(self):
        n = self._v.shape[0]
        u_sum = np.einsum('iijj->', self._v)
        u_avg = u_sum/(n*n)
        j_sum = np.einsum('ijji->', self._v)
        j_avg = (u_sum-j_sum)/max(1.0, n*(n-1))
        if n > 1:
            j_avg = u_avg - j_avg
        self._u_avg = u_avg
        self._j_avg = j_avg

    def unitary_transform(self, u):
        self._v2 = apply_a4_unitary_trans(self._v2, u)

    def print_nnz(self, log=sys.stdout):
        nnz = np.count_nonzero(np.abs(self._v2)>1.e-10)
        print(f"number of nonzero elements in u-matrix = {nnz}.")

    def set_full_u_matrix(self):
        # add spin-components
        self._v2 = add_spin_comp_to_coul_mat(self._v)


class coulomb_matrix_slater(coulomb_matrix):
    '''coulomb matrix for orbitals with quantum number l
    with slater-condo parametrization.
    '''
    def __init__(self, l, slater_integrals):
        self.l = l
        self.slater_integrals = slater_integrals
        super().__init__()

    @classmethod
    def from_uj(cls, l, u_hubb, j_hund):
        fs = cls.get_slater_integrals(u_hubb, j_hund, l=l)
        return cls(l, fs)

    def evaluate(self):
        # Basis of complex spherical harmonics
        # Y_{-2}, Y_{-1}, Y_{0}, Y_{1}, Y_{2}
        # U^{spherical}_{m1 m4 m2 m3} = \sum_{k=0}^{2l} F_k
        # angular_matrix_element(l, k, m1, m2, m3, m4)
        l = self.l
        nl = l*2+1
        self._v = np.zeros((nl,nl,nl,nl),dtype=np.float)
        m_range = range(-self.l, self.l+1)
        for n, f in enumerate(self.slater_integrals):
            k = 2*n
            for m1, m2, m3, m4 in product(m_range, repeat=4):
                self._v[m1+l,m3+l,m2+l,m4+l] += \
                        f*angular_matrix_element(l,k,m1,m2,m3,m4)

    @staticmethod
    def get_slater_integrals(u_hubb, j_hund, l=None, orb=None):
        """get slater radial integrals [f0,f2,f4,f6] given u,j parameters
        for shell of orbital angular momentum l.
        """
        if l is None:
            # orb must be provided.
            l = "spdf".index(orb)
        fs = np.zeros(4,dtype=np.float)
        fs[0] = u_hubb
        if l == 0:
            pass
        elif l == 1:
            fs[1] = j_hund*5.0
        elif l == 2:
            fs[1] = j_hund*14.0/(1.0 + 0.625)
            fs[2] = 0.625*fs[1]
        elif l == 3:
            fs[1] = 6435.0*j_hund/(286.0 + 195.0*0.668 + 250.0*0.494)
            fs[2] = 0.668*fs[1]
            fs[3] = 0.494*fs[1]
        else:
            raise ValueError("l outside of range of [0,1,2,3]!")
        return fs

    @staticmethod
    def get_uj_params(fs, l=None, orb=None):
        '''get u,j parameters given slater radial integrals f0,f1...
        for shell of angular momentum l.
        '''
        if l is None:
            # orb must be provided.
            l = "spdf".index(orb)
        u_hubb = fs[0]
        if l == 0:
            j_hund = 0.
        elif l == 1:
            j_hund = fs[1]/5.0
        elif l == 2:
            j_hund = fs[1]*(1.0 + 0.625)/14.0
        elif l == 3:
            j_hund = fs[1]*(286.0 + 195.0*0.668 + 250.0*0.494)/6435.0
        else:
            raise ValueError("l outside of range of [0,1,2,3]!")
        return u_hubb, j_hund


class coulomb_matrix_kanamori(coulomb_matrix):
    '''coulomb matrix with kanamori parametrization.
    '''
    def __init__(self, n_orb, u_hubb, j_hund):
        self.n = n_orb
        self.u_hubb = u_hubb
        self.j_hund = j_hund
        super().__init__()

    @classmethod
    def from_luj(cls, l, u_hubb, j_hund):
        n_orb = 2*l+1
        return cls(n_orb, u_hubb, j_hund)

    def evaluate(self):
        self._v = np.zeros([self.n]*4, dtype=np.float)
        m_range = range(self.n)
        for m, mp in product(m_range, m_range):
            if m == mp:
                self._v[m, m, mp, mp] = self.u_hubb
            else:
                self._v[m, m, mp, mp] = self.u_hubb - 2.0*self.j_hund
                self._v[m, mp, mp, m] = self.j_hund
                self._v[m, mp, m, mp] = self.j_hund


# Angular matrix elements of particle-particle interaction
# (2l+1)^2 ((l 0) (k 0) (l 0))^2 \sum_{q=-k}^{k} (-1)^{m1+m2+q}
# ((l -m1) (k q) (l m3)) ((l -m2) (k -q) (l m4))
def angular_matrix_element(l, k, m1, m2, m3, m4):
    r"""
    Calculate the angular matrix element

    .. math::
       (2l+1)^2
       \begin{pmatrix}
            l & k & l \\
            0 & 0 & 0
       \end{pmatrix}^2
       \sum_{q=-k}^k (-1)^{m_1+m_2+q}
       \begin{pmatrix}
            l & k & l \\
         -m_1 & q & m_3
       \end{pmatrix}
       \begin{pmatrix}
            l & k  & l \\
         -m_2 & -q & m_4
       \end{pmatrix}.

    Parameters
    ----------
    l : integer
    k : integer
    m1 : integer
    m2 : integer
    m3 : integer
    m4 : integer

    Returns
    -------
    ang_mat_ele : scalar
                  Angular matrix element.

    """
    ang_mat_ele = 0
    for q in range(-k,k+1):
        ang_mat_ele += three_j_symbol((l,-m1),(k,q),(l,m3))* \
                three_j_symbol((l,-m2),(k,-q),(l,m4))* \
                (-1.0 if (m1+q+m2) % 2 else 1.0)
    ang_mat_ele *= (2*l+1)**2 * (three_j_symbol((l,0),(k,0),(l,0))**2)
    return ang_mat_ele


# Wigner 3-j symbols
# ((j1 m1) (j2 m2) (j3 m3))
def three_j_symbol(jm1, jm2, jm3):
    r"""
    Calculate the three-j symbol
    .. math::
       \begin{pmatrix}
        l_1 & l_2 & l_3\\
        m_1 & m_2 & m_3
       \end{pmatrix}.
    Parameters
    ----------
    jm1 : tuple of integers
          (j_1 m_1)
    jm2 : tuple of integers
          (j_2 m_2)
    jm3 : tuple of integers
          (j_3 m_3)
    Returns
    -------
    three_j_sym : scalar
                  Three-j symbol.
    """
    j1, m1 = jm1
    j2, m2 = jm2
    j3, m3 = jm3

    if (m1+m2+m3 != 0 or
        m1 < -j1 or m1 > j1 or
        m2 < -j2 or m2 > j2 or
        m3 < -j3 or m3 > j3 or
        j3 > j1 + j2 or
        j3 < abs(j1-j2)):
        return .0

    three_j_sym = -1.0 if (j1-j2-m3) % 2 else 1.0
    three_j_sym *= np.sqrt(fact(j1+j2-j3)*fact(j1-j2+j3)* \
            fact(-j1+j2+j3)/fact(j1+j2+j3+1))
    three_j_sym *= np.sqrt(fact(j1-m1)*fact(j1+m1)*fact(j2-m2)* \
            fact(j2+m2)*fact(j3-m3)*fact(j3+m3))

    t_min = max(j2-j3-m1,j1-j3+m2,0)
    t_max = min(j1-m1,j2+m2,j1+j2-j3)

    t_sum = 0
    for t in range(t_min,t_max+1):
        t_sum += (-1.0 if t % 2 else 1.0)/(fact(t)*fact(j3-j2+m1+t)* \
                fact(j3-j1-m2+t)*fact(j1+j2-j3-t)*fact(j1-m1-t)*fact(j2+m2-t))

    three_j_sym *= t_sum
    return three_j_sym


# coulomb matrix unitary transformation
def apply_a4_unitary_trans(v2e, u):
    m_range = range(v2e.shape[0])
    v2 = np.asarray(v2e, dtype=u.dtype)
    for i,j in product(m_range, repeat=2):
        v2[i,j,:,:] = u.T.conj().dot(v2[i,j,:,:]).dot(u)
    v2 = v2.swapaxes(0,2).swapaxes(1,3)
    for i,j in product(m_range, repeat=2):
        v2[i,j,:,:] = u.T.conj().dot(v2[i,j,:,:]).dot(u)
    return v2


def add_spin_comp_to_coul_mat(v1):
    # add spin-components with faster spin index
    n = v1.shape[0]
    n2 = n*2
    v2 = np.zeros((n2,n2,n2,n2), dtype=v1.dtype)
    # spin block
    v2[:n,:n,:n,:n] = v2[n:,n:,n:,n:] = v2[:n,:n,n:,n:] = v2[n:,n:,:n,:n] = v1
    return v2
