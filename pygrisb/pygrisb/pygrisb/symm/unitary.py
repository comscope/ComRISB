import numpy as np
from sympy.physics.wigner import clebsch_gordan
from pygrisb.math.matrix_util import random_hermitian_matrix


class basis_trans:
    '''base class to handle basis transformation.
    '''
    def __init__(self):
        self.evaluate()

    @property
    def u(self):
        return self.u_trans

    def evaluate(self):
        raise NotImplementedError("function evaluate not implemented!")


class basis_trans_l(basis_trans):
    '''class to handle various (spin) orbital basis transformation.
    z-component of the angular momentum operator
    is assumed in ascending order unless specified explicitly.
    '''
    def __init__(self, l):
        self.l = l
        super().__init__()

    @classmethod
    def from_dim_m(cls, dim_m):
        l = dim_m//2
        return cls(l)


class relat_harm_to_cubic_relat_harm(basis_trans_l):
    '''class to generate relativistic spherical harmonics to cubic
    relativistic harmonics transformation.
    '''
    def __init__(self, l):
        self.l = l
        super().__init__(self)

    def evaluate(self):
        if self.l == 3:
            jj_to_cubic = np.zeros((14,14))
            jj_to_cubic[0,8] = -np.sqrt(1./6.) # |5/2, -5/2>
            jj_to_cubic[4,8] =  np.sqrt(5./6.) # |5/2, +3/2> G7, 5/2, +
            jj_to_cubic[5,10] = -np.sqrt(1./6.) # |5/2, +5/2>
            jj_to_cubic[1,10] =  np.sqrt(5./6.) # |5/2, -3/2> G7, 5/2, -

            jj_to_cubic[4,0] =  np.sqrt(1./6.) # |5/2, +3/2>
            jj_to_cubic[0,0] =  np.sqrt(5./6.) # |5/2, -5/2> G81, 5/2, +
            jj_to_cubic[1,2] =  np.sqrt(1./6.) # |5/2, -3/2>
            jj_to_cubic[5,2] =  np.sqrt(5./6.) # |5/2, +5/2> G81, 5/2, -

            jj_to_cubic[3,4] =  1. # |5/2, +1/2> G82, 5/2, +
            jj_to_cubic[2,6] =  1. # |5/2, -1/2> G82, 5/2, -

            jj_to_cubic[13,12] =  np.sqrt(5./12.) # |7/2, +7/2>
            jj_to_cubic[ 9,12] =  np.sqrt(7./12.) # |7/2, -1/2> G6, 7/2, +
            jj_to_cubic[ 6,13] =  np.sqrt(5./12.) # |7/2, -7/2>
            jj_to_cubic[10,13] =  np.sqrt(7./12.) # |7/2, +1/2> G6, 7/2, -

            jj_to_cubic[12,11] = -np.sqrt(3./4.) # |7/2, +5/2>
            jj_to_cubic[ 8,11] =  np.sqrt(1./4.) # |7/2, -3/2> G7, 7/2, +
            jj_to_cubic[ 7,9] =  np.sqrt(3./4.) # |7/2, -5/2>
            jj_to_cubic[11,9] = -np.sqrt(1./4.) # |7/2, +3/2> G7, 7/2, -

            jj_to_cubic[13,7] =  np.sqrt(7./12.) # |7/2, +7/2>
            jj_to_cubic[ 9,7] = -np.sqrt(5./12.) # |7/2, -1/2> G81, 7/2, +
            jj_to_cubic[ 6,5] = -np.sqrt(7./12.) # |7/2, -7/2>
            jj_to_cubic[10,5] =  np.sqrt(5./12.) # |7/2, +1/2> G81, 7/2, -

            jj_to_cubic[12,3] = -np.sqrt(1./4.)  # |7/2, +5/2>
            jj_to_cubic[ 8,3] = -np.sqrt(3./4.)  # |7/2, -3/2> G82, 7/2, +
            jj_to_cubic[ 7,1] =  np.sqrt(1./4.)  # |7/2, -5/2>
            jj_to_cubic[11,1] =  np.sqrt(3./4.)  # |7/2, +3/2> G82, 7/2, -
            self.u_trans = jj_to_cubic
        else:
            msg = f'Undefined transformation for l = {self.l}'
            raise ValueError(msg)


class complx_sph_harm_to_real_harm(basis_trans_l):
    '''class to generate complex spherical harmonics to real
    spherical harmonics transformation.
    Condon–Shortley phase convention for complex spherical harmonics.
    See https://en.wikipedia.org/wiki/Spherical_harmonics.
    (VASP uses the same convention.)
    '''
    def evaluate(self):
        dim_m = 2*self.l+1
        csh2rh = np.zeros((dim_m, dim_m), dtype=np.complex)
        # set the orbital indices for the convenience with negative index.
        iy = np.roll(np.arange(dim_m), -self.l)
        for m in range(-self.l,self.l+1):
            if m < 0:
                csh2rh[iy[m], iy[m]] = 1.j/np.sqrt(2.)
                csh2rh[iy[-m], iy[m]] = -1.j/np.sqrt(2.)*(-1)**m
            elif m == 0:
                csh2rh[iy[m], iy[m]] = 1.
            else:
                csh2rh[iy[-m], iy[m]] = 1./np.sqrt(2)
                csh2rh[iy[m], iy[m]] = 1./np.sqrt(2)*(-1)**m
        self.u_trans = csh2rh


def get_u_csh2rh_all(ncorbs_list):
    u_csh2rh_list = [complx_sph_harm_to_real_harm((ncorbs-1)//2).u for \
            ncorbs in ncorbs_list]
    return u_csh2rh_list


class comp_sph_harm_to_relat_harm(basis_trans_l):
    '''class to generate orbital-fast spin-complex spherical harmonics to
    relativistic harmonics transformation.
    '''
    def evaluate(self):
        dim_m = 2*self.l+1
        dim_ms = dim_m*2
        csh2relh = np.zeros((dim_ms, dim_ms), dtype=np.complex)
        iy = np.roll(np.arange(dim_m), -self.l)
        # add slow spin index
        iys = {-0.5:iy, 0.5:[iy[i]+dim_m for i in range(dim_m)]}
        # relativistic_harmonics index
        i_jm = -1
        for i in [-0.5, 0.5]:
            j = self.l+i
            for mj in np.arange(-j, j+1):
                i_jm += 1
                for s in [-0.5,0.5]:
                    csh2relh[iys[s][int(round(mj-s))], i_jm] = \
                            clebsch_gordan(self.l,0.5,j,mj-s,s,mj)
        self.u_trans = csh2relh


def get_u_csh2relh_all(ncorbs_list):
    u_csh2relh_list = [comp_sph_harm_to_relat_harm((ncorbs//2-1)//2).u for \
            ncorbs in ncorbs_list]
    return u_csh2relh_list


class relat_harm_to_spinupdn_comp_sph_harm(basis_trans_l):
    '''class to generate relativistic harmonics to complex spherical harmonics
    with spin block of up and down.
    H. Watanabe 'Operator Methods in Ligand Field Theory'
    Prentice Hall, 1966, Table 1.8-1.
    ordering because of the convention used in WIEN is:
                       mS=1/2        mS=-1/2
                     -L .......L  -L ...... L     (2*(2L+1) columns)
            -(L-1/2)
               .
    J=L-1/2    .
               .
             (L-1/2)
             -L-1/2
               .
    J=L+1/2    .
               .
              L+1/2
    '''
    def evaluate(self):
        L = self.l
        cf = np.zeros((2*(2*L+1), 2*(2*L+1)))
        if L == 0:
            cf[0, 1] = 1.0
            cf[1, 0] = 1.0
        else:
            k1 = -1
            for ms in range(-1, 2, 2):
                ams = -ms/2.
                for ml in range(-L, L+1):
                    k1 = k1+1
                    k2 = -1
                    for mj in range(-2*L+1, 2*L, 2):  # L-1/2 states
                        amj = mj/2.
                        k2 = k2+1
                        d = amj-ml-ams
                        if abs(d) < 0.0001:
                            if ms == 1:
                                cf[k2, k1] = -np.sqrt((L+0.5+amj)/(2*L+1))
                            else:
                                cf[k2, k1] = np.sqrt((L+0.5-amj)/(2*L+1))
                    for mj in range(-2*L-1, 2*L+2, 2):  # L+1/2 states
                        amj = mj/2.
                        k2 = k2+1
                        d = amj-ml-ams
                        if abs(d) < 0.0001:
                            if ms == 1:
                                cf[k2, k1] = np.sqrt((L+0.5-amj)/(2*L+1))
                            else:
                                cf[k2, k1] = np.sqrt((L+0.5+amj)/(2*L+1))
        self.u_trans = cf


class basis_trans_ms(basis_trans):
    '''class to handle various (spin) orbital basis transformation.
    z-component of the angular momentum operator
    is assumed in ascending order unless specified explicitly.
    '''
    def __init__(self, dim_ms):
        self.dim_ms = dim_ms
        super().__init__()

    @classmethod
    def from_l(cls, l):
        dim_ms = 2*(2*l+1)
        return cls(dim_ms)


class orbital_fast_to_spin_fast(basis_trans_ms):
    '''class to generate orbital-fast basis to spin-fast basis transformation.
    '''
    def evaluate(self):
        U = np.zeros((self.dim_ms, self.dim_ms))  # <orbital-fast|spin-fast>
        i = 0
        for j in range(0, self.dim_ms, 2):
            U[i,j] = 1
            i += 1
        for j in range(1, self.dim_ms, 2):
            U[i,j] = 1
            i += 1
        self.u_trans = U


class spin_blk_swap(basis_trans_ms):
    '''class to generate transformation to swap the spin-up, spin-down block.
    '''
    def evaluate(self):
        U = np.zeros((self.dim_ms, self.dim_ms))
        dim_m = self.dim_ms//2
        i = 0
        for j in range(dim_m, self.dim_ms):
            U[i,j] = 1
            i += 1
        for j in range(dim_m):
            U[i,j] = 1
            i += 1
        self.u_trans = U


class random_unitary_trans(basis_trans):
    '''class to generate random unitary transformation.
    '''
    def __init__(self, n):
        self.n = n
        super().__init__()

    def evaluate(self):
        x = random_hermitian_matrix(self.n).a
        w, v = np.linalg.eigh(x)
        self.u_trans = v


class loewdin_unitary_trans(basis_trans):
    '''class to generate random unitary transformation.
    '''
    def __init__(self, a):
        self._a = a
        super().__init__()

    def evaluate(self):
        u, s, v = np.linalg.svd(self._a)
        self.u_trans = u.dot(v)


class ireps_unitary_trans(basis_trans):
    '''class to generate unitary transformation u
    which transforms equivalent irreducible representation a
    to equivalent irreducible representation b, i.e., b=uau^h.
    '''
    def __init__(self, ireps_a, ireps_b):
        self.ireps_a = ireps_a
        self.ireps_b = ireps_b
        super().__init__()

    def evaluate(self):
        m = np.zeros(self.ireps_a[0].shape)
        for irep_a, irep_b in zip(self.ireps_a, self.ireps_b):
            m = m + irep_b.dot(irep_a.T.conj())
        mm = m.dot(m.T.conj())
        vals, vecs = np.linalg.eigh(mm)
        mm = vecs.dot(np.diag(np.sqrt(vals)**(-1))).dot(vecs.T.conj())
        # order is important!
        self.u_trans = mm.dot(m)

    def check_consistency(self):
        for irep_a, irep_b in zip(self.ireps_a, self.ireps_b):
            uauh = self.u_trans.dot(irep_a).dot(self.u_trans.T.conj())
            assert(np.allclose(uauh, irep_b)), "error in ireps_unitary_trans!"


def get_u_csh2wan_all(ncorbs_list):
    ncorbs = ncorbs_list[0]
    if ncorbs%2 == 0: #with spin-orbit interaction
        u_csh2wan_list = get_u_csh2relh_all(ncorbs_list)
    else:
        u_csh2wan_list = get_u_csh2rh_all(ncorbs_list)
    return u_csh2wan_list

