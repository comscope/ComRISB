import numpy, logging
from scipy.linalg import block_diag
import pygrisb.symm.unitary as un


class amop:
    '''base class of angular momentum operater.
    '''
    def __init__(self):
        self.evaluate()

    @property
    def am_op(self):
        return numpy.asarray(self.amop)

    def evaluate(self):
        raise NotImplementedError("function evaluate not implemented!")

    @staticmethod
    def commutation_check(jvec, tol=1.e-5):
        for i in range(0,-3,-1):
            comm = jvec[i].dot(jvec[i+1]) - \
                    jvec[i+1].dot(jvec[i]) - \
                    1.j*jvec[i+2]
            error = numpy.max(numpy.abs(comm))
            assert(error < tol), \
                    "generators not consistent with commutations!"
        logging.info("commutation_check passed.")

    def spin_blk_swap(self):
        u = un.spin_blk_swap(self.amop[0].shape[0]).u
        # u is a symmetric unitary matrix
        self.amop = [u.dot(op).dot(u) for op in self.amop]


class amop_1b(amop):
    '''base class of angular momentum operater in one-particle representation.
    basis orbital assumes a convention of z-component of angular momentum
    operator in ascending order.
    '''
    def __init__(self, l):
        self.l = l
        self.dim_m = int(l*2+1.1)
        super().__init__()


class lamop_csh1(amop_1b):
    '''class of orbital angular momentum operater representation
    in the complex spherical harmonics basis with one spin-component.
    '''
    def set_lopz(self):
        '''
        set matrix l_z.
        '''
        # ascending order.
        m_l = [-self.l+m for m in range(self.dim_m)]
        self.lopz = numpy.diag(m_l)

    def set_lopp(self):
        '''
        set matrix l_{\dagger}.
        '''
        self.lopp = numpy.zeros([self.dim_m, self.dim_m])
        for i in range(self.dim_m-1):
            m = -self.l+i
            self.lopp[i+1, i] = numpy.sqrt(self.l*(self.l+1) - m*(m + 1))

    def evaluate(self):
        self.set_lopz()
        self.set_lopp()
        lopn = self.lopp.T.conj()
        lopx = (self.lopp+lopn)/2.
        lopy = (self.lopp-lopn)/2.j
        self.amop = [lopx, lopy, self.lopz]


class lamop_csh2(amop_1b):
    '''class of orbital angular momentum operater representation
    in the complex spherical harmonics basis with both spin-components.
    '''
    def evaluate(self):
        lamop = lamop_csh1(self.l)
        lvec = lamop.am_op
        self.amop = [block_diag(lop, lop) for lop in lvec]
        self.commutation_check(self.amop)


class samop_csh2(amop_1b):
    '''class of spin angular momentum operater representation
    in the orbital-fast spin complex spherical harmonics basis.
    '''
    def set_sopz(self):
        '''
        set matrix s_z.
        '''
        # ascending order.
        m_s = [-0.5 for m in range(self.dim_m)] + \
                [0.5 for m in range(self.dim_m)]
        self.sopz = numpy.diag(m_s)

    def set_sopp(self):
        '''
        set matrix s_{\dagger}.
        '''
        self.sopp = numpy.zeros([self.dim_m*2, self.dim_m*2])
        s = 0.5
        ms = -0.5
        res = numpy.sqrt(s*(s+1) - ms*(ms+1))
        for i in range(self.dim_m):
            self.sopp[i+self.dim_m, i] = res

    def evaluate(self):
        self.set_sopz()
        self.set_sopp()
        sopn = self.sopp.T.conj()
        sopx = (self.sopp+sopn)/2.
        sopy = (self.sopp-sopn)/2.j
        self.amop = [sopx, sopy, self.sopz]
        self.commutation_check(self.amop)


class jamop_csh2(amop_1b):
    '''class of total angular momentum operater representation
    in the orbital-fast spin complex spherical harmonics basis.
    '''
    def evaluate(self):
        samop = samop_csh2(self.l)
        svec = samop.am_op
        lamop = lamop_csh2(self.l)
        lvec = lamop.am_op
        self.amop = [svec[i]+lvec[i] for i in range(3)]


def get_l_list_from_string(code):
    '''
    Gt l list from letters, e.g., "spd" -> [0, 1, 2].
    '''
    l_list = []
    for c in code:
        l = get_l_from_char(c)
        if l >= 0:
            l_list.append(l)
    return l_list


def get_l_from_char(c):
    return {'s': 0, 'p': 1, 'd': 2, 'f': 3}.get(c, -1)

