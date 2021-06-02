import itertools, sys, logging
import numpy as np
from scipy.sparse import isspmatrix
from pygrisb.math.matrix_util import sym_rnd_hermi_matrix, eigen_space
import pygrisb.symm.angular_momentum as am


class group:
    '''base class for symmetry group.
    '''
    def __init__(self):
        self.evaluate()

    @property
    def elements(self):
        return self.u_list

    def evaluate(self):
        raise NotImplementedError("function evaluate not implemented!")

    def set_product_table(self, tol=1.e-5):
        u_list = np.asarray(self.u_list)
        h = u_list.shape[0]
        table = np.zeros((h, h), dtype=np.int)
        # there may be duplicated representations.
        # but the duplications should be the same.
        deg = 0
        for i,u1 in enumerate(self.u_list):
            for j,u2 in enumerate(self.u_list):
                u12 = u1.dot(u2)
                diff = u_list - u12
                inds = np.where(np.all(np.abs(diff) < tol, axis=(1,2)))[0]
                if deg == 0:
                    if len(inds) == 0:
                        raise ValueError("cannot locate product element. "+ \
                                "check local rotation matrix.")
                    deg = len(inds)
                    if deg > 1:
                        logging.warning(f"representation degeneracy = {deg}!")
                assert(len(inds) == deg), "degeneracies not the same."
                table[i,j] = inds[0]
        self.prod_table = table

    def check_product_table(self):
        self.set_product_table()


class rotation_group(group):
    '''base class for rotation group.
    '''
    @property
    def lie_params(self):
        return self._lie_params

    @property
    def generators(self):
        return self.jvec

    def generator_consistency_check(self, tol=1.e-5):
        am.amop.commutation_check(self.jvec)

    def get_rotation_matrix(self, theta):
        '''
        get rotation matrix givem lie parameters.
        '''
        exponent = self.jvec[0]*theta[0] + self.jvec[1]*theta[1] + \
                self.jvec[2]*theta[2]
        if type(exponent) is np.ndarray:
            from scipy.linalg import expm
        elif isspmatrix(exponent):
            from scipy.sparse.linalg import expm
        else:
            msg = "invalid matrix format: {}".format(type(exponent))
            raise Exception(msg)
        rotation = expm(1.j*exponent)
        return rotation

    def check_double_group(self, tol=1.e-5):
        # check the eigen-values of the jz component.
        w = np.linalg.eigvalsh(self.jvec[2])
        wmin = np.abs(w).min()
        if np.abs(wmin - 0.5) < tol:
            self.idx_dblg = 2
        elif wmin < tol:
            self.idx_dblg = 1
        else:
            msg = f"smallest absolute eigen-value of jz is {wmin}!"
            raise ValueError(msg)


class rotation_group_3d(rotation_group):
    '''class of 3d rotation group.
    '''
    def __init__(self, u_list):
        self.u_list = u_list
        super().__init__()

    @property
    def axis_list(self):
        return self._axis_list

    @property
    def angle_list(self):
        return self.alpha_list

    def evaluate(self):
        self.set_generator()
        self.set_axis_angle_list()
        self.set_lie_params()

    def set_generator(self):
        self.jvec = -1.j*np.array([[[int((i-j)*(j-k)*(k-i)/2)
                for k in range(3)] for j in range(3)] for i in range(3)])

    def set_axis_angle_list(self):
        self._axis_list = []
        self.alpha_list = []
        for u in self.u_list:
            axis, alpha = self.get_axis_angle(u)
            self._axis_list.append(axis)
            self.alpha_list.append(alpha)

    def get_axis_angle(self, u, tol=1.e-7):
        '''given a 3d rotation matrix, calculate the rotation axis and angle.
        '''
        u = np.asarray(u)
        # u is generally nonhermitian,
        # but diagonalizable over the complex field
        vals, vecs = np.linalg.eig(u)
        # pick the label of the eigenvalue 1 (rotation axis)
        vals = np.abs(vals-1)
        axis = vecs[:, vals.argmin()].real
        fix = (u.trace()-1)/2.
        if fix <= -1:
            fix = -1.
        elif fix >= 1:
            fix = 1.
        alpha = np.arccos(fix)
        theta = alpha*axis
        # get the additional sign
        rotation = self.get_rotation_matrix(theta)
        if not np.allclose(rotation, u, atol=tol):
            alpha *= -1
            theta = alpha*axis
            rotation = self.get_rotation_matrix(theta)
            assert(np.allclose(rotation, u, atol=tol)), \
                    "lie params not consistent!"
        return axis, alpha

    def set_lie_params(self):
        self._lie_params = [[alpha*axis, (alpha+2*np.pi)*axis]
                for alpha, axis in zip(self.alpha_list, self._axis_list)]
        self._lie_params = np.array(self._lie_params).swapaxes(0,1)

    def check_lie_params(self):
        '''
        check the representation generated by lie parameters and rotation
        generator jvec against the given rotation matrices.
        '''
        for theta, u in zip(self.lie_params[0], self.u_list):
            u_ = self.get_rotation_matrix(theta)
            assert(np.allclose(u_, u)), "error in_lie parameters!"

    def print_rotations(self, log=sys.stdout):
        print(" 3d rotation list with axis and angle:", file=log)
        for i, u, axis, alpha in zip(itertools.count(), self.u_list, \
                self._axis_list, self.alpha_list):
            print((
                    "{:3d}:\n"
                    "    {:5.2f} {:5.2f} {:5.2f}\n"
                    "    {:5.2f} {:5.2f} {:5.2f}\n"
                    "    {:5.2f} {:5.2f} {:5.2f}\n"
                    "    {:5.2f} {:5.2f} {:5.2f} {:6.1f}").format(
                    i, *u.ravel(), *axis, alpha), file=log)


class rotation_group_md(rotation_group):
    '''class of multi-dimensional rotation group.
    '''
    def __init__(self, jvec, lie_params):
        self.jvec = jvec
        self._lie_params = lie_params
        super().__init__()

    def evaluate(self):
        self.check_double_group()
        self.set_rotation_matrices()

    def set_rotation_matrices(self, rtol=1.e-7):
        self.u_list = [self.get_rotation_matrix(theta) for theta in \
                np.array(self._lie_params[:self.idx_dblg]).reshape(-1,3)]


class point_group_3d(rotation_group_3d):
    '''class of 3d point group.
    '''
    @property
    def non_rotation_generator_type(self):
        '''flag for non_rotation_generator.
        0: pure rotation group; 1: inversion; 2: rotoreflection
        '''
        return self._non_rotation_generator_type

    @property
    def non_rotation_generator_lie_params(self):
        '''for rotoreflection, we need the lie parameter
        for the rotation c. (rotoreflection=c*i)
        '''
        return self._non_rotation_generator_lie_params

    @property
    def pure_rotation_list(self):
        return self.rotation_list

    def evaluate(self):
        self.set_rotation_subgroup_and_nonrot_generator()
        super().evaluate()
        # set self.u_list to self.full_u_list to avoid confusion
        # in the class property of elements.
        self.u_list = self.full_u_list
        self.set_non_rotation_generator_lie_params()

    def set_non_rotation_generator_lie_params(self, tol=1.e-5):
        if self._non_rotation_generator_type == 2:
            axis, alpha = self.get_axis_angle(self.rotoreflection_rot)
            self._non_rotation_generator_lie_params = [\
                    alpha*axis, (alpha+2*np.pi)*axis]
            # make sure the lie parameter is different
            # from the lie parameters of pure rotation elements
            # of the point group.
            lie_params = np.asarray(self._lie_params)[0,:,:]
            diff = lie_params - \
                    np.asarray(self._non_rotation_generator_lie_params[0])
            assert(not np.any(np.all(np.abs(diff) < tol, axis=1))), \
                    "non_rotation_generator_lie_param present duplicated!"
        else:
            # for compatibility reason, dummy values, not to be used.
            self._non_rotation_generator_lie_params = [[0.]*3]*2

    def set_rotation_subgroup_and_nonrot_generator(self, tol=1.e-5):
        '''get the rotation subgroup and a non-rotation generator.
        for 32 point groups in crystalline, only one non-rotation
        (inversion or rotoreflection) generator is needed.
        '''
        inversion_op = -np.eye(3, dtype=np.int)
        u_list = []
        non_rotation_generator_type = 0
        # sort out rotations
        for u in self.u_list:
            if np.abs(np.linalg.det(u)-1) < tol:
                u_list.append(u)
            elif non_rotation_generator_type != 1 and  \
                    np.allclose(u, inversion_op):
                non_rotation_generator_type = 1
            elif non_rotation_generator_type == 0:
                rotoreflection_rot = -u
                non_rotation_generator_type = 2
        self._non_rotation_generator_type = non_rotation_generator_type
        if non_rotation_generator_type == 2:
            self.rotoreflection_rot = rotoreflection_rot
        # point group elements, reordered, copy the list first.
        self.full_u_list = u_list[:]
        if non_rotation_generator_type > 0:
            assert(len(self.u_list)//len(u_list) == 2), \
                    "elements in rotation group and point group not match."
            if non_rotation_generator_type == 2:
                non_rot_generator = -rotoreflection_rot
            else:
                non_rot_generator = inversion_op
            self.full_u_list += [u.dot(non_rot_generator) for u in u_list]
        # pure rotation group elements
        self.rotation_list = self.u_list = u_list


class point_group_md(rotation_group_md):
    '''class for point symmetry group, including rotations, inversion,
    and mirror and rotoreflection,
    which are product of rotation and inversion.
    '''
    def __init__(self, jvec, lie_params, l, ne=1,
            non_rotation_generator_type=0,
            non_rotation_generator_lie_params=None):
        self.l = l
        self.ne = ne
        self._non_rotation_generator_type = non_rotation_generator_type
        self._non_rotation_generator_lie_params = \
                non_rotation_generator_lie_params
        super().__init__(jvec, lie_params)

    def evaluate(self):
        self.check_double_group()
        self.set_rotation_matrices()
        self.generator_consistency_check()
        self.check_product_table()
        # make a copy
        self.rotation_list = self.u_list[:]
        self.add_non_rotations()

    def add_non_rotations(self):
        if self._non_rotation_generator_type == 0:
            return

        self.u_list = []
        n = len(self.rotation_list)//self.idx_dblg
        sign = 1 if self.l%2 == 0 or self.ne%2 == 0 else -1
        inversion = np.eye(self.jvec[0].shape[0])*sign
        for i in range(self.idx_dblg):
            self.u_list.extend(self.rotation_list[i*n:(i+1)*n])
            # add inversion first
            non_rot_generator = inversion.copy()
            if self._non_rotation_generator_type == 2:
                # get rotoreflection
                rot = self.get_rotation_matrix(
                        self._non_rotation_generator_lie_params[i])
                non_rot_generator = non_rot_generator.dot(rot)
            self.u_list.extend([e.dot(non_rot_generator) for e in \
                    self.rotation_list[i*n:(i+1)*n]])

class group_decomp:
    '''class for irreducible representation decomposition of a group
    representation specified u_list.
    '''
    def __init__(self, u_list):
        self.u_list = u_list
        self.order = len(u_list)
        self.dim = u_list[0].shape[0]
        self.evaluate()

    @property
    def unq_chi_list(self):
        return self.unique_chi_list

    @property
    def equ_ireps_list(self):
        return self.equiv_ireps_list

    @property
    def equ_evecs_list(self):
        return self.equiv_evecs_list

    def evaluate(self):
        self.construct_ireps()
        self.group_ireps()
        self.print_chi_space_info()

    def set_ireps(self):
        '''ireps_list dimensions: n_ireps * n_group_elements *
        ired space dim * ired space dim.
        '''
        self.ireps_list = [[vs.T.conj().dot(u).dot(vs) for u in self.u_list] \
                for vs in self.evecs_list]
        self.n_ireps = len(self.ireps_list)

    def set_characters(self):
        self.chi_list = np.array([[irep.trace() for irep in ireps] \
                for ireps in self.ireps_list])

    def check_sum_chi2(self):
        '''
        Check sum of character**2.
        '''
        for chi in self.chi_list:
            sum_chi2 = np.vdot(chi, chi).real
            diff = np.abs(sum_chi2-self.order)
            msg = f"error in sum_chi2 = {diff}!"
            assert(diff < 1.e-5), msg

    def construct_ireps(self):
        while True:
            # symmetrized random hermitian matrix
            smat = sym_rnd_hermi_matrix(self.u_list)
            # get eigen space
            espace = eigen_space.from_matrix(smat)
            self.evecs_list = espace.eigen_space_vectors
            self.set_ireps()
            self.set_characters()
            try:
                self.check_sum_chi2()
                break
            except AssertionError:
                pass

    def group_ireps(self, rtol=1.e-6):
        '''
        Group the equivalent irreducible representations.
        '''
        done_list = []
        self.equiv_ireps_list = []
        self.equiv_evecs_list = []
        self.unique_chi_list = []
        for i1, chi1, ireps1, evecs1 in zip(itertools.count(), self.chi_list,
                self.ireps_list, self.evecs_list):
            if i1 in done_list:
                continue
            equiv_ireps = [ireps1]
            equiv_evecs = [evecs1]
            self.unique_chi_list.append(chi1)
            if i1 < self.n_ireps-1:
                for i2, chi2, ireps2, evecs2 in zip(itertools.count(i1+1),
                        self.chi_list[i1+1:], self.ireps_list[i1+1:],
                        self.evecs_list[i1+1:]):
                    res = chi1.dot(chi2.conj())/self.order
                    if np.abs(res - 1) < rtol:
                        equiv_ireps.append(ireps2)
                        equiv_evecs.append(evecs2)
                        done_list.append(i2)
            self.equiv_ireps_list.append(np.asarray(equiv_ireps))
            self.equiv_evecs_list.append(np.asarray(equiv_evecs))

    def print_chi_space_info(self, log=sys.stdout):
        for i, equiv_evecs in enumerate(self.equiv_evecs_list):
            print(("chi_space {}: {} equivalent ireps \n"
                    "              {} basis vectors.").format(i, \
                    len(equiv_evecs), equiv_evecs[0].shape))
