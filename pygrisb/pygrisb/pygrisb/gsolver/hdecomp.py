import numpy, h5py
from openfermion.linalg.givens_rotations import givens_decomposition_square
from pygrisb.gsolver.forest.util import givens_rotation


class hamildecomp:
    def __init__(self,
            t,  # [spin][[]]
            v,
            ):
        self._t = numpy.asarray(t)
        self._v = v  # chemist's convention v_{pqrs}c^\dag_p c^\dag_r c_s c_q
        # spin dimension
        self._nspin = len(t)
        # spin-orbit dimension
        self._iso = self._t.shape[1]//v.shape[0]
        self.set_h()

    @classmethod
    def from_hdf5(cls, fname="h.h5"):
        with h5py.File(fname, "r") as f:
            t = f["/t"][()]
            v = f["/v"][()]
        return cls(t, v)

    def set_h(self):
        # merge one-body part from coulomb v-tensor
        # due to operator reordering from c+c+cc -> c+c c+c
        s = numpy.einsum("ijkl->il", self._v)
        self._h = [t - 0.5*s for t in self._t]

    def set_eigendecomp_h(self):
        self._lambda_h = []
        self._givens_h = []
        for h in self._h:
            w, v = numpy.linalg.eigh(self._h)
            self._lambda_h.append(w)
            givens, _diagonal = givens_decomposition_square(v)
            # _diagonal (+/-1) not matter since only n_i is needed.
            self._givens_h.append(givens)

    def set_eigendecomp_u(self, tol=1.e-4):
        n = self._v.shape[0]
        vs = self._v.reshape((n*n, n*n))
        w, v = numpy.linalg.eigh(vs)
        self._lambda_v = []
        self._givens_v = []
        for w1, v1 in zip(w, v.T):  # column of v is the eigen-vector
            if w1 < -tol:
                raise ValueError(f"eigen-value of v = {w1:.2e} < 0.")
            elif w1 < tol:
                continue
            else:
                v1 *= numpy.sqrt(w)
                v1 = v1.resahpe((n, n))
                w2, v2 = numpy.linalg.eigh(v1)
                self._lambda_v.append(w2)
                givens, _diagonal = givens_decomposition_square(v2)
                self._givens_v.append(givens)


class hamildecompforest(hamildecomp):
    '''
    rigetti forest.
    '''
    def set_eigendecomp_h(self):
        super().set_eigendecomp_h()
        nq = self._v.shape[0]
        givens_prog = []
        for isp, givens in enumerate(self._givens_h):
            circuit = None
            prog = givens_rotation(givens, nq, isp)
            if circuit == None:
                circuit = prog
            else:
                circuit += prog
            givens_prog.append(circuit)
        self._givens_h = givens_prog

    def set_eigendecomp_u(self, tol=1.e-4):
        super().set_eigendecomp_u(tol=tol)
        nq = self._v.shape[0]
        givens_prog = []
        for givens in self._givens_v:
            circuit = None
            for isp in range(2):
                prog = givens_rotation(givens, nq, isp)
                if circuit == None:
                    circuit = prog
                else:
                    circuit += prog
            givens_prog.append(circuit)
        self._givens_v = givens_prog


