from __future__ import print_function

import unittest
import numpy as np
from scipy.linalg import expm,sqrtm
from pygrisb.math.matrix_util import yield_derivative_f_matrix, \
        get_1stderivative_sroot_entropyn_at_h, \
        get_2ndderivative_sroot_entropyn_at_h


class KnowValues(unittest.TestCase):
    def test_derivative_exp_ix(self):
        # random Hermitian matrix of dimention 5.
        n = 5
        x = np.random.rand(n, n) + np.random.rand(n, n)*1.j
        x += np.conj(x.T)
        # Numerically calculate the derivative w.r.t. h.
        h = np.array([ \
            [0,              0, 1./np.sqrt(2.), 0, 0],
            [0,              0, 0,              0, 0],
            [1./np.sqrt(2.), 0, 0,              0, 0],
            [0,              0, 0,              0, 0],
            [0,              0, 0,              0, 0]])

        expix0 = expm(-1.j*x)
        delta = 1.e-7
        expix1 = expm(-1.j*(x + delta*h))
        partial0 = (expix1 - expix0)/delta
        f = lambda x: np.exp(-1.j*x)
        fp = lambda x: -1.j*np.exp(-1.j*x)
        partials = yield_derivative_f_matrix(x, [h], f, fp)
        for partial in partials:
            err = np.abs(np.max(partial - partial0))
            self.assertAlmostEqual(err, 0.)

    def test_derivative_srootn(self):
        # random Hermitian matrix of dimention 5.
        n = 5
        x = np.random.rand(n, n) + np.random.rand(n, n)*1.j
        x += np.conj(x.T)
        _,v = np.linalg.eigh(x)
        w = np.random.rand(n)
        # get a valid density matrix
        x = v.dot(np.diag(w)).dot(v.T.conj())

        # Numerically calculate the derivative w.r.t. h.
        h = np.array([ \
            [0,              0, 1./np.sqrt(2.), 0, 0],
            [0,              0, 0,              0, 0],
            [1./np.sqrt(2.), 0, 0,              0, 0],
            [0,              0, 0,              0, 0],
            [0,              0, 0,              0, 0]])
        delta = 1.e-7
        xp = x + delta*h
        sroot0 = sqrtm((np.eye(x.shape[0])-x).dot(x))
        sroot1 = sqrtm((np.eye(x.shape[0])-xp).dot(xp))
        partial0 = (sroot1 - sroot0)/delta
        partial1 = get_1stderivative_sroot_entropyn_at_h(x, h)

        f = lambda x: np.sqrt((1-x)*x)
        fp = lambda x: 0.5/np.sqrt((1-x)*x)*(1-2*x)
        partials = yield_derivative_f_matrix(x, [h], f, fp)

        for partial2 in partials:
            err = np.abs(np.max(partial2 - partial0))
            self.assertAlmostEqual(err, 0., places=6)
        err = np.abs(np.max(partial1 - partial0))
        self.assertAlmostEqual(err, 0., places=6)


    def test_2ndderivative_srootn(self):
        # random Hermitian matrix of dimention 5.
        n = 5
        x = np.random.rand(n, n) + np.random.rand(n, n)*1.j
        x += np.conj(x.T)
        _,v = np.linalg.eigh(x)
        w = np.random.rand(n)
        # get a valid density matrix
        x = v.dot(np.diag(w)).dot(v.T.conj())

        # Numerically calculate the derivative w.r.t. h.
        h1 = np.array([ \
            [0,              0, 1./np.sqrt(2.), 0, 0],
            [0,              0, 0,              0, 0],
            [1./np.sqrt(2.), 0, 0,              0, 0],
            [0,              0, 0,              0, 0],
            [0,              0, 0,              0, 0]])
        h2 = np.array([ \
            [0,              0, 0, 1./np.sqrt(2.), 0],
            [0,              0, 0,              0, 0],
            [0,              0, 0,              0, 0],
            [1./np.sqrt(2.), 0, 0,              0, 0],
            [0,              0, 0,              0, 0]])
        delta = 1.e-7
        xp = x + delta*h2
        p0 = get_1stderivative_sroot_entropyn_at_h(x, h1)
        p1 = get_1stderivative_sroot_entropyn_at_h(xp, h1)
        partial0 = (p1 - p0)/delta

        partial1 = get_2ndderivative_sroot_entropyn_at_h(x, h1, h2)
        err = np.abs(np.max(partial1 - partial0))
        self.assertAlmostEqual(err, 0., places=6)


    def test_derivative_tr_exp_ix_a(self):
        # random Hermitian matrix of dimention 5.
        n = 5
        a = np.random.rand(n, n)
        x = np.random.rand(n, n) + np.random.rand(n, n)*1.j
        x += np.conj(x.T)
        # Numerically calculate the derivative w.r.t. h.
        h = np.array([ \
            [0,              0, 1./np.sqrt(2.), 0, 0],
            [0,              0, 0,              0, 0],
            [1./np.sqrt(2.), 0, 0,              0, 0],
            [0,              0, 0,              0, 0],
            [0,              0, 0,              0, 0]])

        expix0 = expm(-1.j*x)
        delta = 1.e-7
        expix1 = expm(-1.j*(x + delta*h))
        partial0 = (np.trace(np.dot(expix1, a)) -
                np.trace(np.dot(expix0, a)))/delta
        f = lambda x: np.exp(-1.j*x)
        fp = lambda x: -1.j*np.exp(-1.j*x)
        partials = yield_derivative_f_matrix(x, [h], f, fp)
        for partial in partials:
            # trace
            _partial = np.sum(partial.T*a)
            err = np.abs(_partial - partial0)
            self.assertAlmostEqual(err, 0.)


if __name__=="__main__":
    print("Tests for matrix_util.")
    unittest.main()
