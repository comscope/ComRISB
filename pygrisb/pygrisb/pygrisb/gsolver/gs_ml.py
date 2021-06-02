'''machine learning solver originally developed by John Rogers,
with contributions from Nicola Lanata, Tsung-Han Lee, and Yongxin Yao.
For any questions or concerns please email jsr12e@my.fsu.edu
'''
import numpy, h5py, argparse, warnings
from sklearn.metrics.pairwise import rbf_kernel
from pygrisb.gsolver.gsolver import gsolver_h5ml_fsoc


class gsolver_h5ml_fsoc_krr(gsolver_h5ml_fsoc):
    def evaluate(self, tol=1.e-10):
        self.set_ref_data_path()
        with h5py.File(f"{self.ref_path}/kernelridge/fsoc_training_data.h5",
                "r") as f:
            weights = f[f"/val_{self.nval}/weights_laplace"][()]
            vtrain = f[f"/val_{self.nval}/v_training"][()]
            gamma = f["/gamma"][()]
        gaussians = numpy.asarray([numpy.exp(-gamma*(numpy.sum(
                [numpy.sqrt(v1**2+tol) for v1 in self.v-v]))) for v in vtrain])
        self.res = gaussians.dot(weights).reshape(-1)
        self.get_density_matrix()
        # self.check_unique_elements()
        if self._jac:
            # derivative
            derivatives = []
            for i, v1 in enumerate(self.v):
                dcoefs = numpy.asarray([-gamma*(v1-rv1)/numpy.sqrt(
                        (v1-rv1)**2+tol)
                        for rv1 in vtrain[:,i]])
                gaussiansd = gaussians*dcoefs
                derivative = gaussiansd.dot(weights).reshape(-1)
                derivatives.append(derivative)
            self.derivatives = numpy.asarray(derivatives)
            self.get_density_matrix_derivatives()

    def sanity_check(self):
        assert(1<=self.nval<=6), "valence out of range [1,6]!"


class gsolver_h5ml_fsoc_nm1d(gsolver_h5ml_fsoc):
    '''n-mode expansion at 1st order for multivariate interpolation.
    '''
    def evaluate(self, mode=1):
        self.set_ref_data_path()
        with h5py.File(f"{self.ref_path}/normalmode/n{self.nval}/"+
                "fsoc_training_data_1d.h5", "r") as f:
            f0 = f["/f0"][()].reshape(-1)
            ene0 = f["/ene0"][()]
            f1b = numpy.zeros((5,6))
            ene1b = numpy.zeros(5)
            if self._jac:
                f1b_der = numpy.zeros((5,6))
            for i, vi in enumerate(self.v):
                # load array [n_samples, n_features]
                deltas = f[f"/axis_{i}/deltas"][()].reshape(-1)
                if vi < deltas.min() or vi > deltas.max():
                    obj = ["D1", "D2", "Delta", "Delta1", "Delta2"][i]
                    warnings.warn((f"ML solver: {obj} = {vi:.4f} "
                            f"out of range [{deltas.min()}, {deltas.max()}]!"))
                # load energy prediction weights and gamma
                weights = f[f"/axis_{i}/weights_ene"][()].reshape(-1)
                gamma = f[f"/axis_{i}/gammas_ene"][()]
                # gaussians for prediction
                # gaussians = rbf_kernel(vi, deltas, gamma=gamma)
                gaussians = numpy.exp(-gamma*(vi-deltas)**2)
                # energy bar function
                ene1b[i] = gaussians.dot(weights) - ene0

                # density matrix
                for j in range(6):
                    # load dm prediction weights and gamma
                    weights = f[f"/axis_{i}/weights_{j}"][()].reshape(-1)
                    gamma = f[f"/axis_{i}/gammas_{j}"][0,0]
                    # gaussians for prediction
                    # gaussians = rbf_kernel(vi, deltas, gamma=gamma)
                    gaussians = numpy.exp(-gamma*(vi - deltas)**2)
                    # dm bar function
                    f1b[i, j] = gaussians.dot(weights) - f0[j]
                    if self._jac:
                        dcoeffs = -2*gamma*(vi - deltas)
                        dgaussians = dcoeffs*gaussians
                        f1b_der[i, j] = dgaussians.dot(weights)

        self.res = (f0 + numpy.sum(f1b, axis=0)).reshape(-1)
        self._emol = ene0 + numpy.sum(ene1b)
        # convention due to the learning database.
        # energy correction due to the constant potential shift.
        self._emol += (self.e1 + self.e2 - (self.l1 + self.l2))/4* \
                self.h1e.shape[0]
        if mode > 0:
            self.get_density_matrix()
            if self._jac:
                self.derivatives = self._dm1b_der
                self.get_density_matrix_derivatives()
            # self.check_unique_elements()
        else:
            self._ene1b = ene1b
            self._dm1b = f1b
            if self._jac:
                self._dm1b_der = f1b_der

    def sanity_check(self):
        assert(1<=self.nval<=6), "valence out of range [1,6]!"


class gsolver_h5ml_fsoc_nm2d(gsolver_h5ml_fsoc_nm1d):
    '''n-mode expansion at 2nd order for multivariate interpolation.
    '''
    def evaluate(self, mode=1):
        super().evaluate(mode=0)
        with h5py.File( \
                f"{self.ref_path}/normalmode/n{self.nval}/"+ \
                "fsoc_training_data_1d.h5", "r") as f1, \
                h5py.File( \
                f"{self.ref_path}//normalmode/n{self.nval}/"+ \
                "fsoc_training_data_2d.h5", "r") as f2:
            f0 = f1["/f0"][()].reshape(-1)
            ene0 = f1["/ene0"][()]
            # setup 2d predictions
            dm2b = numpy.zeros((10,6))
            ene2b = numpy.zeros(10)
            if self._jac:
                dm2b_der = numpy.zeros((5,6))
            #loop over each axis to make 2d axis predictions
            ij = 0
            for i in range(5):
                for j in range(i+1,5):
                    # 2d predicting
                    deltas = f2[f"/axis_{i}_{j}/deltas"][()].T
                    gammas = f2[f"/axis_{i}_{j}/gammas_ene"][()]
                    vij = numpy.array([self.v[i], self.v[j]]).reshape(-1,1).T
                    gaussians = rbf_kernel(vij, deltas, gamma=gammas)
                    weights = f2[f"/axis_{i}_{j}/weights_ene"][()]
                    ene2b[ij] = gaussians.dot(weights) - ene0 \
                            - self._ene1b[i] - self._ene1b[j]
                    # density matrix
                    for k in range(6):
                        # load dm prediction weights and gamma
                        gammas = f2[f"/axis_{i}_{j}/gammas_{k}"][()]
                        # gaussians for prediction
                        gaussians = rbf_kernel(vij, deltas, gamma=gammas)
                        weights = f2[f"/axis_{i}_{j}/weights_{k}"][()]
                        # dm bar function
                        dm2b[ij,k] = gaussians.dot(weights) - f0[k]
                        if self._jac:
                            dcoeffs = -2*gammas*(vij[:,0]-deltas[:,0])
                            dgaussians = dcoeffs*gaussians
                            derivative = dgaussians.dot(weights)
                            dm2b_der[i, k] += derivative
                            dcoeffs = -2*gammas*(vij[:,1]-deltas[:,1])
                            dgaussians = dcoeffs*gaussians
                            derivative = dgaussians.dot(weights)
                            dm2b_der[j, k] += derivative
                    # calculate 2d contribution to nm expansion
                    dm2b[ij,:] -= self._dm1b[i,:] + self._dm1b[j,:]
                    if self._jac:
                        dm2b_der[i, :] -= self._dm1b_der[i,:]
                        dm2b_der[j, :] -= self._dm1b_der[j,:]
                    ij += 1
        # add 2d correction to 1d nm prediction
        self.res += numpy.sum(dm2b, axis=0)
        self._emol += numpy.sum(ene2b)
        if mode > 0:
            self.get_density_matrix()
            if self._jac:
                self.derivatives = self._dm1b_der + dm2b_der
                self.get_density_matrix_derivatives()
            # self.check_unique_elements()
        else:
            self._ene2b = ene2b
            self._dm2b = dm2b
            if self._jac:
                self._dm2b_der = dm2b_der


class gsolver_h5ml_fsoc_nm3d(gsolver_h5ml_fsoc_nm2d):
    '''n-mode expansion at 3rd order for multivariate interpolation.
    '''
    def evaluate(self):
        super().evaluate(mode=0)
        with h5py.File( \
                f"{self.ref_path}/normalmode/n{self.nval}/"+ \
                "fsoc_training_data_1d.h5", "r") as f1, \
                h5py.File( \
                f"{self.ref_path}//normalmode/n{self.nval}/"+ \
                "fsoc_training_data_3d.h5", "r") as f3:
            f0 = f1["/f0"][()].reshape(-1)
            ene0 = f1["/ene0"][()]
            # setup 3d predictions
            dm3b = numpy.zeros((10,6))
            ene3b = numpy.zeros(10)
            if self._jac:
                dm3b_der = numpy.zeros((5,6))
            #loop over each axis to make 2d axis predictions
            ijk = 0
            for i in range(5):
                for j in range(i+1,5):
                    ij = get_p12_index(5, i, j)
                    for k in range(j+1,5):
                        ik = get_p12_index(5, i, k)
                        jk = get_p12_index(5, j, k)
                        # 3d predicting
                        deltas = f3[f"/axis_{i}_{j}_{k}/deltas"][()].T
                        gammas = f3[f"/axis_{i}_{j}_{k}/gammas_ene"][()]
                        vijk = numpy.array([self.v[i], self.v[j],  \
                                self.v[k]]).reshape(-1,1).T
                        gaussians = rbf_kernel(vijk, deltas, gamma=gammas)
                        weights = f3[f"/axis_{i}_{j}_{k}/weights_ene"][()]
                        ene3b[ijk] = gaussians.dot(weights) - ene0 \
                                - self._ene1b[i] \
                                - self._ene1b[j] \
                                - self._ene1b[k] \
                                - self._ene2b[ij] \
                                - self._ene2b[ik] \
                                - self._ene2b[jk]
                        # density matrix
                        for l in range(6):
                            # load dm prediction weights and gamma
                            gammas = f3[f"/axis_{i}_{j}_{k}/gammas_{l}"][()]
                            # gaussians for prediction
                            gaussians = rbf_kernel(vijk, deltas, gamma=gammas)
                            weights = f3[f"/axis_{i}_{j}_{k}/weights_{l}"][()]
                            # dm bar function
                            dm3b[ijk,l] = gaussians.dot(weights) - f0[l]
                            if self._jac:
                                for m, n in enumerate([i, j, k]):
                                    dcoeffs = -2*gammas*(vijk[:,m]-deltas[:,m])
                                    dgaussians = dcoeffs*gaussians
                                    derivative = dgaussians.dot(weights)
                                    dm3b_der[n, l] += derivative
                        # calculate 2d contribution to nm expansion
                        dm3b[ijk,:] -= self._dm1b[i,:]  \
                                + self._dm1b[j,:] \
                                + self._dm1b[k,:] \
                                + self._dm2b[ij,:] \
                                + self._dm2b[ik,:] \
                                + self._dm2b[jk,:]
                        if self._jac:
                            for n in [i, j, k]:
                                dm3b_der[n, :] -= self._dm1b_der[n,:] + \
                                        self._dm2b_der[n,:]
                        ijk += 1
        # add 3d correction to 1d nm prediction
        self.res += numpy.sum(dm3b, axis=0)
        self._emol += numpy.sum(ene3b)
        self.get_density_matrix()
        if self._jac:
            self.derivatives = self._dm1b_der + self._dm2b_der + dm3b_der
            self.get_density_matrix_derivatives()
        # self.check_unique_elements()


# convenience functions
def get_p12_index(n, i, j):
    return n*i - (1+i)*i//2 + j -i - 1


def get_inline_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--impurity", type=int, default=0,
            help="impurity index (int)")
    parser.add_argument("-n", "--valence", type=int, default=5,
            help="valence of the atom (int)")
    parser.add_argument("-m", "--mode", type=int, default=0,
            help="mode: 0 (ml); 1 (n-mode 1d); else (n-mode 2d)")
    parser.add_argument("-d", "--delta", type=float, default=1.e-4,
            help="finite difference step (float)")
    args = parser.parse_args()
    return args


def gs_h5ml_fsoc():
    args = get_inline_args()
    if args.mode == 0:
        gs_h5ml_fsoc_krr(nval=args.valence, imp=args.impurity)
    elif args.mode == 1:
        gs_h5ml_fsoc_nm1d(nval=args.valence, imp=args.impurity)
    else:
        gs_h5ml_fsoc_nm2d(nval=args.valence, imp=args.impurity)


def check_jacobian():
    args = get_inline_args()
    if args.mode == 0:
        gs = gsolver_h5ml_fsoc_krr(nval=args.valence, imp=args.impurity)
    elif args.mode == 1:
        gs = gsolver_h5ml_fsoc_nm1d(nval=args.valence, imp=args.impurity)
    else:
        gs = gsolver_h5ml_fsoc_nm2d(nval=args.valence, imp=args.impurity)
    gs.check_jacobian(delta=args.delta)


def gs_h5ml_fsoc_krr(nval=5, imp=1, jac=None):
    gs = gsolver_h5ml_fsoc_krr(nval, imp=imp, jac=jac)
    gs.driver()


def gs_h5ml_fsoc_nm1d(nval=5, imp=1, jac=None):
    gs = gsolver_h5ml_fsoc_nm1d(nval, imp=imp, jac=jac)
    gs.driver()


def gs_h5ml_fsoc_nm2d(nval=5, imp=1, jac=None):
    gs = gsolver_h5ml_fsoc_nm2d(nval, imp=imp, jac=jac)
    gs.driver()


def gs_h5ml_fsoc_nm3d(nval=5, imp=1, jac=None):
    gs = gsolver_h5ml_fsoc_nm3d(nval, imp=imp, jac=jac)
    gs.driver()
