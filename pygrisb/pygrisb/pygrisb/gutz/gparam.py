import numpy as np
import h5py
from pygrisb.math.matrix_basis import dense_matrix_basis
from pygrisb.math.matrix_util import list_mat_to_array

class gparam:
    '''class to write/read gutzwiller simulation parameter files.
    '''
    def __init__(self, kwargs: dict):
        self.set_default_params1()
        self.update(kwargs)
        self.set_default_params2()

    def update(self, kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)

    def set_default_params1(self):
        raise NotImplementedError("function evaluate not implemented!")

    def set_default_params2(self):
        raise NotImplementedError("function evaluate not implemented!")

    def h5compat_list_save(self, h5compat_list, h5f):
        for key in h5compat_list:
            try:
                if key.lower() == "symie":
                    h5f["/"].attrs[key.upper()] = getattr(self, key)
                else:
                    h5f["/"].attrs[key] = getattr(self, key)
            except:
                pass

    @staticmethod
    def h5load(fname):
        ps = {}
        with h5py.File(fname, 'r') as f:
            for key, val in f["/"].attrs.items():
                try:
                    ps[key] = val[()]
                except:
                    pass
            if "/impurity_0" in f:
                for i in range(ps["num_imp"]):
                    for key, val in f[f"/impurity_{i}"].items():
                        if key in ps:
                            ps[key].append(val[()])
                        else:
                            ps[key] = [val[()]]
        return ps


class gparam_atms(gparam):
    '''class to write/read GParam_atms.h5 file,
    which contains the essential parameters for correlated atoms (impurities).
    '''
    def __init__(self, na2_list, kwargs: dict):
        self.na2_list = na2_list
        super().__init__(kwargs)

    def set_default_params1(self):
        self.dtype = np.complex
        self.iso = 1
        self.ispin = 1
        self.giembeddiag = 0
        self.gmaxiter = 10000
        self.dc_mode = 0
        self.dc_u_avg_list = None
        self.dc_j_avg_list = None
        self.symbol_matrix_list = None
        self.matrix_basis_list = None
        self.matbs_dim_list = None
        self.v2e_list = None
        self.ityp_list = None
        self.imap_list = None

    def set_default_params2(self):
        self.num_imp = len(self.na2_list)
        if self.dc_u_avg_list is None:
            self.dc_u_avg_list = [0. for _ in self.na2_list]
        if self.dc_j_avg_list is None:
            self.dc_j_avg_list = [0. for _ in self.na2_list]
        if self.symbol_matrix_list is None:
            self.symbol_matrix_list = [np.arange(1,n**2+1).reshape((n,n)) \
                    for n in self.na2_list]
        if self.matrix_basis_list is None:
            self.matrix_basis_list = [\
                    dense_matrix_basis(symbol_matrix).basis \
                    for symbol_matrix in self.symbol_matrix_list]
        if self.matbs_dim_list is None:
            self.matbs_dim_list = [len(matbs) for matbs in
                    self.matrix_basis_list]
        matrix_basis_r_list = getattr(self, "matrix_basis_r_list", None)
        if matrix_basis_r_list is not None:
            self.matbs_r_dim_list = [len(matbs) for matbs in
                    matrix_basis_r_list]
        if self.v2e_list is None:
            self.v2e_list = [np.zeros((n,n,n,n), dtype=self.dtype) \
                    for n in self.na2_list]
        if self.ityp_list is None:
            self.ityp_list = np.arange(self.num_imp)
        if self.imap_list is None:
            self.imap_list = np.arange(self.num_imp)

    def h5save(self):
        with h5py.File('GParam.h5', 'w') as f:
            # root entries
            h5compat_list = [
                    "iso",
                    "unit",
                    "spin_order",
                    "ispin",
                    "giembeddiag",
                    "gmaxiter",
                    "dc_mode",
                    "na2_list",
                    "dc_u_avg_list",
                    "dc_j_avg_list",
                    "ityp_list",
                    "imap_list",
                    "num_imp",
                    "dc_nelf_list",
                    "nval_top_list",
                    "nval_bot_list",
                    "matbs_dim_list",
                    "matbs_r_dim_list",
                    ]
            self.h5compat_list_save(h5compat_list, f)

            # impurity-specific entries.
            h5sub_list = [
                    "sx",
                    "sy",
                    "sz",
                    "lx",
                    "ly",
                    "lz",
                    "db2sab",
                    "symm_operations_3d",
                    "lie_params",
                    "lie_params2",
                    "symm_operations_csh",
                    "symbol_matrix",
                    "matrix_basis",
                    "matrix_basis_r",
                    "v2e",
                    ]
            for key in h5sub_list:
                val_list = getattr(self, f"{key}_list", None)
                if val_list is not None:
                    for i, val in enumerate(val_list):
                        if key in ["v2e"]:
                            # dump in fortran order
                            f[f'/impurity_{i}/{key.upper()}'] = val.T
                        else:
                            f[f'/impurity_{i}/{key}'] = val

    @staticmethod
    def params():
        return super().h5load("GParam.h5")

    def add_vext_list(self, vext_list, iext=1):
        ''' add external local potential.
        iext = 0: initial perturbation;
               1: additional terms added to h1e and init_la1
              -1: additional terms added to h1e, but not init_la1
        '''
        with h5py.File('GParam.h5', 'a') as f:
            if "/vext" in f:
                del f["/vext"]
            for i, vext in enumerate(vext_list):
                f[f"/vext/impurity_{i}/v"] = vext
            f["/vext"].attrs["givext"] = iext

    def add_shft_init_la1(self, delta):
        ''' add initial diagonal constant shift lambda
        in quasiparticle hamiltonian.
        '''
        with h5py.File('GParam.h5', 'a') as f:
            f["/"].attrs["shft_init_la1"] = delta


class gparam_bnds(gparam):
    '''class to write/read GBareH.h5 file,
    which contains the essential parameters for correlated atoms (impurities).
    '''
    def __init__(self, kptwt, nelectron, nbmax, **kwargs):
        self.kptwt = kptwt
        self.nelectron = nelectron
        self.nbmax = nbmax
        super().__init__(kwargs)

    def set_default_params1(self):
        self.dtype = np.complex
        self.iso = 1
        self.ispin = 1
        self.symnop = 1
        self.symie = 1
        self.ensemble = 0
        self.ismear = 0
        self.delta = 0.01
        self.ne_list = None
        self.h1e_list = None

    def set_default_params2(self):
        self.kptdim = len(self.kptwt)
        if self.ne_list is None:
            self.ne_list = [[self.nbmax, 0, self.nbmax] for w in self.kptwt]
        if self.h1e_list is None:
            # assume 1 bnad case
            self.h1e_list = [np.zeros((2, 2), dtype=self.dtype)]
        self.num_imp = len(self.h1e_list)

    def h5save(self):
        h5compat_list = ["iso", "ispin", "kptdim", "nbmax", "kptwt", \
                "ismear", "ensemble", "delta", "nelectron", "symnop", \
                "symie", "ne_list"]
        with h5py.File("GBareH.h5", "w") as f:
            self.h5compat_list_save(h5compat_list, f)
            h1e_list = list_mat_to_array(self.h1e_list)
            f["/ispin_0/H1E_LIST"] = h1e_list.swapaxes(1,2)

    @staticmethod
    def params():
        return super().h5load("GBareH.h5")
