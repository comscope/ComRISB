import h5py, numpy, os, pprint, subprocess, warnings
from scipy.linalg import block_diag


class gsolver:
    '''base class to solve Gutzwiller embedding Hamiltonian.
    '''
    def __init__(self, h1e, h2e):
        # one-body and two-body part of the e,bedding Hamiltonian
        self.h1e = numpy.asarray(h1e)
        self.h2e = numpy.asarray(h2e)
        # initialize the density matrix and energy
        self._dm = None
        self._emol = 0.
        self._derivatives_d = None
        self._derivatives_l = None

    @property
    def dm(self):
        # one-particle density matrix
        return self._dm

    @property
    def emol(self):
        # expectation value of the Hamiltonian (molecular energy).
        return self._emol

    def show_results(self):
        print(f"energy = {self._emol}")
        print("density matrix: ")
        pprint.pprint(self._dm)

    def driver(self):
        raise NotImplementedError("driver not implemented!")


class gsolver_h5(gsolver):
    '''base class to solve Gutzwiller embedding Hamiltonian
    with hdf5 interface.
    '''
    def __init__(self,
            imp=0,
            jac=None,
            read_v2e=False,
            ):
        # impurity index
        self.imp = imp
        self._emol = 0.
        self._jac = jac
        self.read_v2e = read_v2e

    def driver(self):
        self.evaluate()
        self.save_results()

    def evaluate(self):
        raise NotImplementedError("function not implemented!")

    def set_parameters(self):
        with h5py.File('HEmbed.h5', 'r') as f:
            self.d = f[f'impurity_{self.imp}/D'][()].T
            self.h1e = f[f'impurity_{self.imp}/H1E'][()].T
            self.lam = f[f'impurity_{self.imp}/LAMBDA'][()].T
            if self.read_v2e:
                self.v2e = f[f'impurity_{self.imp}/V2E'][()].T
            else:
                self.v2e = None

    def set_full_v1e(self):
        norb = self.h1e.shape[0]
        self.v1e = block_diag(self.h1e, -self.lam)
        self.v1e[:norb, norb:] = self.d.T
        self.v1e[norb:, :norb] = self.d.conj()


    def set_full_v12e(self):
        self.set_full_v1e()
        norb = self.h1e.shape[0]
        norb2 = norb*2
        v2e = numpy.zeros((norb2, norb2, norb2, norb2), dtype=self.v2e.dtype)
        v2e[:norb, :norb, :norb, :norb] = self.v2e
        self.v2e = v2e

    def save_results(self):
        with h5py.File(f'HEmbedRes_{self.imp}.h5', 'w') as f:
            f['/DM'] = self._dm.T
            f['/'].attrs["emol"] = self._emol
            if self._jac:
                na2 = self._dm.shape[0]//2
                if self._derivatives_d is not None:
                    # bath density matrix (nb_{ab} = <f_b f_a^\dagger>
                    # (nc_var) = \delta_{ab}-<f_a^\dagger f_b>)
                    f['/NCV_DERI_D'] = -self._derivatives_d[:,na2:,na2:]. \
                            swapaxes(1,2)
                    # r0_{aA} = <c_A^\dagger f_a>, a natural transposition
                    f['/R0_DERI_D'] = self._derivatives_d[:,:na2,na2:]
                if self._derivatives_l is not None:
                    f['/NCV_DERI_LC'] = -self._derivatives_l[:,na2:,na2:]. \
                            swapaxes(1,2)
                    f['/R0_DERI_LC'] = self._derivatives_l[:,:na2,na2:]


class gsolver_h5ed(gsolver_h5):
    '''parallel exact diagonalization solver.
    '''
    def __init__(self, imp=0, mpiexec=["mpirun", "-np", "2"], path="./", \
            nval_bot=0, nval_top=None):
        super().__init__(imp)
        self.mpiexec = mpiexec
        self.path = path
        self._nval_bot = nval_bot
        self._nval_top = nval_top
        if not os.path.isfile(f"GEDInp_{self.imp}.h5"):
            self.read_v2e = True

    def evaluate(self):
        self.set_parameters()
        self.set_full_v1e()
        self.set_ged_inp()
        cmd = self.mpiexec + [f"{self.path}/exe_ed", "-i", str(self.imp)]
        print(f" running {' '.join(cmd)}")
        subprocess.run(cmd)
        self.load_ged_out()

    def set_ged_inp(self, rtol=1.e-10):
        with h5py.File(f"GEDInp_{self.imp}.h5", "a") as f:
            if self.v2e is not None:
                ind_list = numpy.where(numpy.abs(self.v2e) > rtol)
                # absorb the 1/2 factor.
                val_list = self.v2e[ind_list]/2.
                # to fortran convention
                ind_list = numpy.array(ind_list)+1
                f["/v2e/INDS"] = ind_list.T
                f["/v2e/vals"] = val_list
                f["/v2e/nnz"] = [len(val_list)]
            if "/fock_basis" not in f:
                from pygrisb.mbody.basis import double_fock_basis \
                        as dfock_basis
                fb = dfock_basis(self.h1e.shape[0], nval_bot=self._nval_bot, \
                        nval_top=self._nval_top)
                fb.calc()
                fb.h5save(f=f, path="/fock_basis")

            ind_list = numpy.where(numpy.abs(self.v1e) > rtol)
            val_list = self.v1e[ind_list]
            if "/v1e" in f:
                del f["/v1e"]
            # to fortran convention.
            ind_list = numpy.array(ind_list)+1
            f["/v1e/INDS"] = ind_list.T
            f["/v1e/vals"] = val_list
            f["/v1e/norb"] = self.v1e.shape[0]
            f["/v1e/nnz"] = [len(val_list)]

    def load_ged_out(self):
        # get the calculation results in cmp_sph_harm basis with
        # faster spin index.
        with h5py.File(f"GEDOut_{self.imp}.h5", "r") as f:
            self._dm = f["/DM"][()].T
            self._emol = f["/emol"][0]


class gsolver_h5trans_ed(gsolver_h5ed):
    '''parallel exact diagonalization solver with hamiltonian transformed
    into complex spherical hamrmonics basis with spin-faster index.
    '''
    def evaluate(self):
        self.set_parameters()
        self.get_csh2sab()
        self.trans_v12_tocsh()
        super().evaluate()
        self.trans_dm_tosab()

    # get complex harmonics (spin-faster) basis
    # to symmetry adapted basis transformation
    def get_csh2sab(self):
        with h5py.File('GParam.h5', 'r') as f:
            self._u = f[f'/impurity_{self.imp}/db2sab'][()]
        from pygrisb.symm.unitary import orbital_fast_to_spin_fast as of2sf
        u_of2sp = of2sf(self._u.shape[0]).u
        self._u = u_of2sp.T.dot(self._u)

    def trans_v12_tocsh(self):
        self.d = self._u.conj().dot(self.d).dot(self._u.T)
        self.lam = self._u.dot(self.lam).dot(self._u.T.conj())
        self.h1e = self._u.dot(self.h1e).dot(self._u.T.conj())
        if self.v2e is not None:
            from pygrisb.mbody.coulomb_matrix import apply_a4_unitary_trans
            self.v2e = apply_a4_unitary_trans(self.v2e, self._u.T.conj())

    def trans_dm_tosab(self):
        # append the transformation for the bath
        u2_sf2sab = block_diag(self._u, self._u)
        # unitary transformation to symmetry-adapted basis
        self._dm = u2_sf2sab.T.dot(self._dm).dot(u2_sf2sab.conj())


class gsolver_h5ml(gsolver_h5):
    '''machine learning eigen solver class.
    '''
    def __init__(self, imp, jac=None, db_path=None):
        if db_path is None:
            self.set_ref_data_path()
        else:
            self.ref_path = db_path
        super().__init__(imp, jac=jac)

    def driver(self):
        self.set_parameters()
        self.sanity_check()
        self.reduce_parameters()
        self.evaluate()
        self.save_results()

    def sanity_check(self):
        raise NotImplementedError("function not implemented!")

    def set_ref_data_path(self):
        if "GMLDB_ROOT" not in os.environ:
            raise FileNotFoundError("Please set environ var GMLDB_ROOT!")
        else:
            self.ref_path = os.environ["GMLDB_ROOT"]


class gsolver_h5ml_fsoc(gsolver_h5ml):
    '''machine learning eigen solver class for f-eletron system gutzwiller
    embedding hamiltonian with spin-orbit interactin only.
    '''
    def __init__(self, nval, imp=1, jac=None, db_path=None):
        self.nval = nval
        self.l = 3
        super().__init__(imp, jac=jac, db_path=db_path)

    def reduce_parameters(self, tol=1.e-6):
        # unique parameters for h_embed with spin-orbit interaction only.
        self.e1 = self.h1e[0,0].real
        self.e2 = self.h1e[6,6].real
        self.l1 = self.lam[0,0].real
        self.l2 = self.lam[6,6].real
        # convert d to negative real convention.
        d1 = self.d[0,0]
        # convention is d1 <= 0
        self.d1_phase = -d1/numpy.abs(d1)
        self.d1 = (d1*numpy.conj(self.d1_phase)).real
        d2 = self.d[6,6]
        self.d2_phase = -d2/numpy.abs(d2)
        self.d2 = (d2*numpy.conj(self.d2_phase)).real
        if abs(self.d1_phase.imag) > tol or abs(self.d2_phase.imag) > tol or \
                self.d1_phase.real < 0 or self.d2_phase.real < 0:
            warnings.warn("Positive or complex D detected!")
        # Reparameterize
        delta = (self.e1+self.e2+self.l1+self.l2)/4.0
        delta1 = (self.e1-self.e2)/2.0
        delta2 = (self.l1-self.l2)/2.0
        # Build vector representing point to predict
        self.v = numpy.array([self.d1,self.d2,delta,delta1,delta2])

    def get_density_matrix(self):
        self._dm = self.recover_full_matrix(self.res)

    def recover_full_matrix(self, vec, mode="dm"):
        na2 = (2*self.l+1)*2
        na4 = na2*2
        dm = numpy.zeros([na4, na4], dtype=self.d.dtype)
        # j = 5/2 block
        idx = numpy.arange(2*self.l)
        dm[idx, idx] = vec[0]
        dm[idx+na2, idx+na2] = vec[2]
        dm[idx, idx+na2] = vec[4]*numpy.conj(self.d1_phase)
        dm[idx+na2, idx] = vec[4]*self.d1_phase
        # j = 7/2 block
        idx = numpy.arange(2*self.l, na2)
        dm[idx, idx] = vec[1]
        dm[idx+na2, idx+na2] = vec[3]
        dm[idx, idx+na2] = vec[5]*numpy.conj(self.d2_phase)
        dm[idx+na2, idx] = vec[5]*self.d2_phase
        # normalize dm to half-filling case.
        dm -= dm.trace()/na4*numpy.eye(na4)
        if mode == "dm":
            dm += 0.5*numpy.eye(na4)
        return dm

    def recover_full_matrix_list(self, vec_list):
        dm_list = []
        for vec in vec_list:
            dm = self.recover_full_matrix(vec, mode="derivative")
            dm_list.append(dm)
        return numpy.asarray(dm_list)

    def get_density_matrix_derivatives(self, tol=1.e-8):
        assert(abs(self.d1_phase.imag)<tol and abs(self.d2_phase.imag)<tol), \
                "error: not real hembed version."
        # wrt d, rescaled to matrix basis.
        self.derivatives[0,:] /= numpy.sqrt(6.)
        self.derivatives[1,:] /= numpy.sqrt(8.)
        derivatives_d = self.derivatives[:2,:]
        self._derivatives_d = self.recover_full_matrix_list(derivatives_d)
        # wrt lambda
        res1 = self.derivatives[2,:]/4. + self.derivatives[4,:]/2.
        res2 = self.derivatives[2,:]/4. - self.derivatives[4,:]/2.
        res1 /= numpy.sqrt(6.)
        res2 /= numpy.sqrt(8.)
        derivatives_l = numpy.asarray([res1, res2])
        self._derivatives_l = self.recover_full_matrix_list(derivatives_l)

    def check_unique_elements(self):
        print(self._dm[[0,6,14,20,0,6], [0,6,14,20,14,20]].real)

    def check_jacobian(self, delta=1.e-4):
        # if gaussian function with variables on exponents,
        # step has to be relatively bigger for finite difference somehow.
        self._jac = True
        # get analytical jacobian
        self.driver()
        # keep a copy of dm
        dm_ref = self._dm
        self._jac = False
        # check derivative_d numerically
        # d1 (j=5/2)
        self.d[0,0] += 1./numpy.sqrt(6.)*delta
        self.check_jacobian1(dm_ref, self._derivatives_d[0], delta)
        self.d[0,0] -= 1./numpy.sqrt(6.)*delta
        # d2 (j=7/2)
        self.d[6,6] += 1./numpy.sqrt(8.)*delta
        self.check_jacobian1(dm_ref, self._derivatives_d[1], delta)
        self.d[6,6] -= 1./numpy.sqrt(8.)*delta
        # lambda_c (j=5/2)
        self.lam[0,0] += 1./numpy.sqrt(6.)*delta
        self.check_jacobian1(dm_ref, self._derivatives_l[0], delta)
        self.lam[0,0] -= 1./numpy.sqrt(6.)*delta
        # lambda_c (j=7/2)
        self.lam[6,6] += 1./numpy.sqrt(8.)*delta
        self.check_jacobian1(dm_ref, self._derivatives_l[1], delta)
        self.lam[6,6] -= 1./numpy.sqrt(8.)*delta

    def check_jacobian1(self, dm_ref, jac_ref, delta):
        self.reduce_parameters()
        self.evaluate()
        jac = (self._dm - dm_ref)/delta
        # print("reference gradients")
        # for row in jac_ref:
        #     print("".join(f" {e:.2f}" for e in row))
        # print("numerical gradients")
        # for row in jac:
        #     print("".join(f" {e:.2f}" for e in row))
        print(f" max jac diff = {numpy.max(numpy.abs(jac_ref-jac)):.2e}")
