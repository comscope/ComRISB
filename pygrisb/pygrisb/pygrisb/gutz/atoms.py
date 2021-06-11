import os, sys, h5py, itertools, spglib, glob, json
import numpy as np
import pygrisb.mbody.coulomb_matrix as c_matrix
import pygrisb.math.matrix_basis as mb
import pygrisb.symm.basis as sb
import pygrisb.symm.angular_momentum as am
import pygrisb.symm.group as gp
from pygrisb.gutz.molecule import cell_point_groups_analyzer
from pygrisb.basic.units import rydberg_to_ev as rydberg

class structure:
    """atomic structure class.
    """
    def __init__(self, symbols, scaled_positions, cell, case=None,
                 locrot_list=None):
        self.symbols = symbols
        self.scaled_positions = np.asarray(scaled_positions)
        self.cell = np.asarray(cell)
        self.case = case
        self.locrot_list = locrot_list
        self.complete_struct_params()
        self.ps = {}

    @classmethod
    def from_numbers(cls, numbers, scaled_positions, cell, case=None,
            locrot_list=None):
        atom_symbols = [*atom_data.keys()]
        symbols = [atom_symbols[i] for i in numbers]
        return cls(symbols, scaled_positions, cell, case=case,
                locrot_list=locrot_list)

    @property
    def z_numbers(self):
        return self.numbers

    @property
    def volume(self):
        if not hasattr(self, "_volume"):
            self._volume = np.linalg.det(self.cell)
        return self._volume

    def complete_struct_params(self):
        self.set_symbols()
        self.set_positions()
        self.set_default_equivalent_atom_indices()

    def set_symbols(self):
        atom_symbols = [*atom_data.keys()]
        if not hasattr(self, "numbers"):
            self.numbers = np.array([atom_symbols.index(s)
                    for s in self.symbols])

    def set_positions(self):
        if not hasattr(self, "positions"):
            self.positions = self.scaled_positions.dot(self.cell)

    def set_default_equivalent_atom_indices(self):
        # Getting general information about symmetry of the lattice
        if self.case is not None:
            from pygrisb.dft.wien import get_equivalent_atom_indices
            self.default_idx_equivalent_atoms = \
                    get_equivalent_atom_indices(self.case)
        else:
            cell = (self.cell, self.scaled_positions, self.numbers)
            dataset = spglib.get_symmetry_dataset(cell, symprec=1e-5)
            self.default_idx_equivalent_atoms = \
                    dataset['equivalent_atoms'].tolist()

    def update(self, ps):
            self.ps.update(ps)


class cor_structure(structure):
    '''dtructure class with local electron correlations.
    '''
    def set_all_attributes(self, fixsab=False, realhemb=False, valtrunc=False):
        self.set_properties_unique_toall(valtrunc=valtrunc)
        self.set_sl_vector_csh2_list()
        # for symmetry-adapted basis
        if self.ps["orbital_polarization"] == 1:
            self.set_site_symm_operations_3d()
            self.set_lie_params_list()
            self.set_symm_operations_csh_list()

        if fixsab:
            self.load_sab_list()
        else:
            self.set_symm_adaptive_basis()
            self.set_full_symm_adaptive_basis()

        self.set_matrix_basis_list(realhemb=realhemb)
        self.set_coulomb_list(realhemb=realhemb)
        self.set_sl_vector_sab_list()

    def set_properties_unique_toall(self, valtrunc=False):
        cor_atm_list = []
        ityp_list = []
        imap_list = []
        if self.ps["dc_mode"] == 2:
            dc_nelf_list = []
        cor_orb_list = []
        u_list = []
        j_list = []
        fs_list = []
        idx_equivalent_atoms = list(self.ps["idx_equivalent_atoms"])
        for i, s in enumerate(self.symbols):
            if s in self.ps["unique_cor_symbol_list"]:
                cor_atm_list.append(i)
                ityp_list.append(self.ps["unique_cor_symbol_list"].index(s))
                cor_orb_list.append(self.ps["unique_cor_orb_list"]\
                        [ityp_list[-1]])

                if self.ps["dc_mode"] == 2:
                    dc_nelf_list.append(self.ps["unique_nelf_list"] \
                            [ityp_list[-1]])

                u_list.append(self.ps["unique_u_list"][ityp_list[-1]])
                j_list.append(self.ps["unique_j_list"][ityp_list[-1]])
                fs_list.append(self.ps["unique_fs_list"][ityp_list[-1]])
                idx_equ = idx_equivalent_atoms[i]
                imap = idx_equivalent_atoms.index(idx_equ)
                imap_list.append(imap)
        self.cor_atm_list = cor_atm_list
        self.cor_orb_list = cor_orb_list
        self.l_list = ["spdf".index(cor_orb) for cor_orb in cor_orb_list]
        self.u_list = u_list
        self.j_list = j_list
        self.fs_list = fs_list
        na2_list = [l*4+2 for l in self.l_list]
        if valtrunc:
            nval_top_list = [atom_data[self.symbols[i]][1] \
                    for i in cor_atm_list]
            nval_bot_list = [atom_data[self.symbols[i]][0] \
                    for i in cor_atm_list]
        else:
            nval_top_list = na2_list
            nval_bot_list = [0]*len(cor_atm_list)
        imap_list = cleanup_imap_list(imap_list)
        self.ps.update({"ityp_list": ityp_list,
                "imap_list": imap_list,
                "na2_list": na2_list,
                "nval_top_list": nval_top_list,
                "nval_bot_list": nval_bot_list})
        if self.ps["dc_mode"] == 2:
            self.ps["dc_nelf_list"] = dc_nelf_list

    def load_sab_list(self):
        csh2sab_list = []
        symbol_matrix_list = []
        with h5py.File("GParam.h5", "r") as f:
            for i, _ in enumerate(self.ps["imap_list"]):
               csh2sab_list.append(
                       f[f"/impurity_{i}/db2sab"][()])
               symbol_matrix_list.append(
                       f[f"/impurity_{i}/symbol_matrix"][()])
        self.ps["db2sab_list"] = csh2sab_list
        self.ps["symbol_matrix_list"] = symbol_matrix_list

    def set_symm_adaptive_basis(self):
        '''get symmetry-adaptive basis transformation from (orbital-fast)
        complex spherical harmonics basis and the associated
        symbolic matrix for general matrices commutative
        with the symmetry operations (crystal field).
        '''
        csh2sab_list = []
        symbol_matrix_list = []
        iso = self.ps["iso"]
        ispin = self.ps["ispin"]
        orb_pol = self.ps["orbital_polarization"]
        for i, l, na2, imap in zip(itertools.count(), self.l_list, \
                self.ps["na2_list"], self.ps["imap_list"]):
            # check equivalent atoms
            if imap == i:
                naso = na2//(3-iso)
                # if not explicitly consider crystal discrete symmetry
                if orb_pol != 1:
                    if self.ps["orbital_polarization"] == 0:
                        if iso == 1:
                            # uniform case, no splitting.
                            sab = sb.symm_basis_uniform(naso)
                        else:
                            sab = sb.symm_basis_soc(l, ispin)
                    elif self.ps["orbital_polarization"] == 2:
                        sab = sb.symm_basis_random(naso)
                else:
                    # check internal consistency.
                    symm_ops = self.ps["symm_operations_csh_list"][i]
                    assert(symm_ops[0].shape[0]==naso), \
                            "conflict between symm_operations_csh and iso!"
                    sab = sb.symm_basis_pg(symm_ops)
                if iso == 2 and self.ps["spin_order"] == -1:
                    sab.spin_blk_swap()
                csh2sab_list.append(sab.transform_matrix)
                symbol_matrix_list.append(sab.symbolic_matrix)
            else:
                csh2sab_list.append(csh2sab_list[imap])
                symbol_matrix_list.append(symbol_matrix_list[imap])
        self.ps["db2sab_list"] = csh2sab_list
        self.ps["symbol_matrix_list"] = symbol_matrix_list

    def set_full_symm_adaptive_basis(self):
        if self.ps["iso"] == 2:
            # already included spin-component.
            return
        # in the case without spin-orbit interaction,
        # we prefer spin-fast indices for the convenience of
        # adding correlated local quantities of multiple sites.
        csh2sab_list = []
        symbol_matrix_list = []
        for i, csh2sab, symbol_matrix, na2, imap in zip(itertools.count(),
                self.ps["db2sab_list"], self.ps["symbol_matrix_list"],
                self.ps["na2_list"], self.ps["imap_list"]):
            if imap == i:
                na = na2//2
                # two identical spin-block.
                csh2sab2 = np.zeros((na2,na2), dtype=np.complex)
                csh2sab2[:na,::2] = csh2sab2[na:,1::2] = csh2sab
                symbol_matrix2 = np.zeros((na2,na2), dtype=np.int)
                symbol_matrix2[::2,::2] = symbol_matrix
                if self.ps["ispin"] == 1:
                    symbol_matrix2[1::2,1::2] = symbol_matrix
                else:
                    ishift = symbol_matrix.max()
                    symbol_matrix[symbol_matrix>0] += ishift
                    symbol_matrix2[1::2,1::2] = symbol_matrix
            else:
                csh2sab2 = csh2sab_list[imap]
                symbol_matrix2 = symbol_matrix_list[imap]
            csh2sab_list.append(csh2sab2)
            symbol_matrix_list.append(symbol_matrix2)
        self.ps["db2sab_list"] = csh2sab_list
        self.ps["symbol_matrix_list"] = symbol_matrix_list

    def set_matrix_basis_list(self, realhemb=False):
        matrix_basis_list = []
        matrix_basis_r_list = []
        for symbol_matrix, na2 in zip(self.ps["symbol_matrix_list"],
                self.ps["na2_list"]):
            assert(symbol_matrix.shape[0] == na2), \
                    "need to run set_full_symm_adaptive_basis first!"
            if realhemb:
                smb = mb.dense_matrix_basis(symbol_matrix=symbol_matrix,
                        dtype=np.float, btype='general')
                matrix_basis_r_list.append(smb.basis)
                smb = mb.dense_matrix_basis(symbol_matrix=symbol_matrix,
                        dtype=np.float)
            else:
                smb = mb.dense_matrix_basis(symbol_matrix=symbol_matrix)
            matrix_basis_list.append(smb.basis)
        self.ps["matrix_basis_list"] = matrix_basis_list
        if len(matrix_basis_r_list) > 0:
            self.ps["matrix_basis_r_list"] = matrix_basis_r_list

    def set_coulomb_list(self, realhemb=False):
        '''set coulomb matrix in orbital-fast
        complex spherical harmonics basis.
        '''
        v2e_list = []
        u_avg_list = []
        j_avg_list = []
        for l, u, j, fs, db2sab in zip(self.l_list, self.u_list, self.j_list, \
                self.fs_list, self.ps["db2sab_list"]):
            # convert energy unit from rydberg to ev
            if "ryd" in self.ps["unit"].lower():
                u /= rydberg
                j /= rydberg
            if self.ps["u_matrix_type"] == 1:
                cm = c_matrix.coulomb_matrix_slater.from_uj(l, u, j)
            elif self.ps["u_matrix_type"] == 2:
                cm = c_matrix.coulomb_matrix_kanamori.from_luj(l, u, j)
            elif self.ps["u_matrix_type"] == 3:
                cm = c_matrix.coulomb_matrix_slater(l, fs)
            else:
                raise NotImplementedError("u_matrix_type = "+\
                        f"{self.ps['u_matrix_type']} not implemented.")
            u_avg_list.append(cm.u_avg)
            j_avg_list.append(cm.j_avg)
            # transform to symmetry-adapted basis
            cm.unitary_transform(db2sab)
            v2e_list.append(cm.u_spin_orb)
        if realhemb:
            for v2e in v2e_list:
                if np.max(np.abs(v2e.imag)) > 1.e-6:
                    raise AssertionError("real version with complex v2e!")
            v2e_list = [v2e.real for v2e in v2e_list]
        self.ps.update({"dc_u_avg_list": u_avg_list,
                "dc_j_avg_list": j_avg_list,
                "v2e_list": v2e_list})

    def set_sl_vector_csh2_list(self):
        '''set spin and orbital angular momentum operator in the orbital-fast
        complex spherical harmonics basis.
        '''
        s_vec_list = []
        l_vec_list = []
        for l in self.l_list:
            samop = am.samop_csh2(l)
            if self.ps["spin_order"] == -1:
                # spin up-dn convention
                samop.spin_blk_swap()
            s_vec_list.append(samop.am_op)
            lamop = am.lamop_csh2(l)
            if self.ps["spin_order"] == -1:
                # spin up-dn convention
                lamop.spin_blk_swap()
            l_vec_list.append(lamop.am_op)
        self.s_vec_csh2_list = s_vec_list
        self.l_vec_csh2_list = l_vec_list

    def set_s_vec_csh2_list(self, s_vec_list):
        '''for external call.
        '''
        self.s_vec_csh2_list = s_vec_list

    def set_sl_vector_sab_list(self, mode='all'):
        '''set spin and orbital angular momentum operator
        in the symmetry-adaptive basis.
        '''
        self.ps["sz_list"] = [u.T.conj().dot(s[2]).dot(u) for s,u in \
                zip(self.s_vec_csh2_list, self.ps["db2sab_list"])]
        if mode == "sz_only":
            return
        self.ps["sx_list"] = [u.T.conj().dot(s[0]).dot(u) for s,u in \
                zip(self.s_vec_csh2_list, self.ps["db2sab_list"])]
        self.ps["sy_list"] = [u.T.conj().dot(s[1]).dot(u) for s,u in \
                zip(self.s_vec_csh2_list, self.ps["db2sab_list"])]
        if hasattr(self, "l_vec_csh2_list"):
            self.ps["lx_list"] = [u.T.conj().dot(l[0]).dot(u) for l,u in \
                    zip(self.l_vec_csh2_list, self.ps["db2sab_list"])]
            self.ps["ly_list"] = [u.T.conj().dot(l[1]).dot(u) for l,u in \
                    zip(self.l_vec_csh2_list, self.ps["db2sab_list"])]
            self.ps["lz_list"] = [u.T.conj().dot(l[2]).dot(u) for l,u in \
                    zip(self.l_vec_csh2_list, self.ps["db2sab_list"])]

    def set_lie_params_list(self):
        lie_params_list = []
        non_rotation_generator_type_list = []
        non_rotation_generator_lie_params_list = []
        # symm_operations_3d may get reordered.
        symm_operations_3d_list = []
        for sym_ops in self.ps["symm_operations_3d_list"]:
            pg3d = gp.point_group_3d(sym_ops)
            lie_params_list.append(pg3d.lie_params)
            non_rotation_generator_type_list.append(\
                    pg3d.non_rotation_generator_type)
            non_rotation_generator_lie_params_list.append(\
                    pg3d.non_rotation_generator_lie_params)
            symm_operations_3d_list.append(pg3d.elements)
        self.ps["lie_params_list"] = lie_params_list
        self.ps["non_rotation_generator_type_list"] = \
                non_rotation_generator_type_list
        self.ps["non_rotation_generator_lie_params_list"] = \
                non_rotation_generator_lie_params_list
        self.ps["symm_operations_3d_list"] = symm_operations_3d_list

    def set_symm_operations_csh_list(self):
        # if not consider crystal-field effect.
        if self.ps["orbital_polarization"] != 1:
            return
        symm_operations_csh_list = []
        for svec, lvec, lie_params, non_rotation_generator_type, \
                non_rotation_generator_lie_params, l \
                in zip(self.s_vec_csh2_list, \
                self.l_vec_csh2_list, \
                self.ps["lie_params_list"], \
                self.ps["non_rotation_generator_type_list"], \
                self.ps["non_rotation_generator_lie_params_list"], \
                self.l_list):
            if self.ps["iso"] == 2:
                jvec = svec + lvec
            else:
                m = 2*l+1
                jvec = lvec[:, :m, :m]
            pgmd = gp.point_group_md(jvec, lie_params, l,
                    non_rotation_generator_type=\
                    non_rotation_generator_type,\
                    non_rotation_generator_lie_params=\
                    non_rotation_generator_lie_params)
            symm_operations_csh_list.append(pgmd.elements)
        self.ps["symm_operations_csh_list"] = symm_operations_csh_list

    def set_site_symm_operations_3d(self):
        cpgs = cell_point_groups_analyzer(self)
        self.ps["symm_operations_3d_list"] = cpgs.symm_ops_list

    def search_struct(self, log=sys.stdout):
        # check case.struct file first.
        struct_files = glob.glob("*.struct")
        n_struct = len(struct_files)
        cif_files = glob.glob("*.cif")
        n_cif = len(cif_files)
        locrot_list = None
        case = None
        if n_cif > 1:
            raise ValueError("multiple cif files present!")
        elif n_struct > 1:
            raise ValueError("multiple struct files present!")
        elif n_struct == 1:
            material, case = read_wien_structure()
            from pygrisb.dft.wien import get_local_rotations
            locrot_list = get_local_rotations(f'{get_wiencase()}.struct')
        elif n_cif == 1:
            cif_file = cif_files[0]
            material = read_cif(cif_file)
            print(f"structure info read from {cif_file}.", file=log)
        elif os.path.isfile("POSCAR"):
            material = read_vasp_poscar()
        else:
            raise ValueError("missing structure file(cif, struct or POSCAR)!")
        self.symbols = material.get_chemical_symbols()
        self.cell= material.get_cell()
        self.scaled_positions = material.get_scaled_positions()
        self.case = case
        self.locrot_list = locrot_list


class disk_cor_structure(cor_structure):
    '''dtructure class with local electron correlations.
    '''
    def __init__(self, fname):
        self.fname = fname
        try:
            self.load_struct()
        except:
            self.search_struct()
            self.save_struct()
        self.complete_struct_params()
        self.ps = {}

    def load_struct(self):
        raise NotImplementedError("load_struct not implemented.")

    def save_struct(self):
        raise NotImplementedError("save_struct not implemented.")


class h5_cor_structure(disk_cor_structure):
    '''structure class utilizing hdf5 file for storage.
    '''
    def __init__(self, fname="ginit.h5"):
        super().__init__(fname)

    def load_struct(self):
        '''not very good using hdf5, which is designed to store large dataset.
        '''
        with h5py.File(self.fname, "r") as f:
            symbols = f['/struct/symbols'][()]
            self.symbols = [s.decode() for s in symbols]
            self.cell = f['/struct/cell'][()]
            self.scaled_positions = f['/struct/scaled_positions'][()]
            if '/struct/case' in f:
                self.case = f['/struct/case'][()]
            else:
                self.case = None
            if '/struct/locrot_list' in f:
                self.locrot_list = f['/struct/locrot_list'][()]
            else:
                self.locrot_list = None

    def save_struct(self):
        with h5py.File(self.fname, 'w') as f:
            f['/struct/symbols'] = self.symbols
            f['/struct/cell'] = self.cell
            f['/struct/scaled_positions'] = self.scaled_positions
            if self.case is not None:
                f['/struct/case'] = self.case


class json_cor_structure(disk_cor_structure):
    '''structure class utilizing hdf5 file for storage.
    '''
    def __init__(self, fname="ginit.json"):
        super().__init__(fname)

    def load_struct(self):
        data = json.load(open(self.fname, "r"))
        # get the group of struct.
        data = data["struct"]
        self.symbols = data["symbols"]
        self.cell = np.array(data["cell"])
        self.scaled_positions = np.array(data["scaled_positions"])
        self.case = data.get("case", None)
        self.locrot_list = data.get("locrot_list", None)

    def save_struct(self):
        data = {
                "symbols": self.symbols,
                "cell": self.cell.tolist(),
                "scaled_positions": self.scaled_positions.tolist()}
        if self.case is not None:
            data["case"] = self.case
        if self.locrot_list is not None:
            data["locrot_list"] = self.locrot_list
        data = {"struct": data}
        with open(self.fname, 'w') as f:
            json.dump(data, f, indent=4)


def read_cif(cif_file):
    from pymatgen.core import Structure
    struct = Structure.from_file(cif_file, primitive=True)
    import pymatgen.symmetry as symmetry
    struct = symmetry.analyzer.SpacegroupAnalyzer( \
            struct).get_primitive_standard_structure()
    from  pymatgen.io.ase import AseAtomsAdaptor
    material = AseAtomsAdaptor.get_atoms(struct)
    return material


def read_vasp_poscar():
    from ase.io.vasp import read_vasp
    material = read_vasp(filename='POSCAR')
    return material


def read_wien_structure():
    structure_file = get_wiencase() + '.struct'
    from ase.io.wien2k import read_struct
    material = read_struct(structure_file)
    return material, structure_file


def get_wiencase():
    # Determine WIEN2k case name
    case = os.path.basename(os.getcwd())
    # directory name is not case (happens when submitting to cluster)
    if not os.path.isfile(case + '.struct'):
        files = glob.glob('*.struct')
        if len(files) < 1:
            raise IOError('No struct file present!')
        elif len(files) > 1:
            raise IOError('Multiple struct files present!')
        else:
            case = os.path.splitext(os.path.basename(files[0]))[0]
    return case


def cleanup_imap_list(imap_list):
    '''Clean up imap list in case correlated atoms are not ordered first
    in the struct file.
    '''
    imax = np.max(imap_list)+1
    index_map = np.zeros(imax, dtype=np.int)-1
    for i in range(imax):
        if i in imap_list:
            index_map[i] = imap_list.index(i)
    _imap_list = []
    for imap in imap_list:
        _imap_list.append(index_map[imap])
    return _imap_list



atom_data = {
        "X": [0, 0],    "H": [0, 2],
        "He":[0, 2],    "Li":[0, 2],
        "Be":[0, 2],    "B": [0, 6],
        "C": [0, 6],    "N": [0, 6],
        "O": [0, 6],    "F": [0, 6],
        "Ne":[0, 6],    "Na":[0, 2],
        "Mg":[0, 2],    "Al":[0, 6],
        "Si":[0, 6],    "P": [0, 6],
        "S": [0, 6],    "Cl":[0, 6],
        "Ar":[0, 6],    "K": [0, 2],
        "Ca":[0, 2],    "Sc":[0,10],
        "Ti":[0,10],    "V": [0,10],
        "Cr":[0,10],    "Mn":[0,10],
        "Fe":[0,10],    "Co":[0,10],
        "Ni":[0,10],    "Cu":[0,10],
        "Zn":[0,10],    "Ga":[0, 6],
        "Ge":[0, 6],    "As":[0, 6],
        "Se":[0, 6],    "Br":[0, 6],
        "Kr":[0, 6],    "Rb":[0, 2],
        "Sr":[0, 2],    "Y": [0,10],
        "Zr":[0,10],    "Nb":[0,10],
        "Mo":[0,10],    "Tc":[0,10],
        "Ru":[0,10],    "Rh":[0,10],
        "Pd":[0,10],    "Ag":[0,10],
        "Cd":[0,10],    "In":[0, 6],
        "Sn":[0, 6],    "Sb":[0, 6],
        "Te":[0, 6],    "I": [0, 6],
        "Xe":[0, 6],    "Cs":[0, 2],
        "Ba":[0, 2],    "La":[0, 2],
        "Ce":[0, 3],    "Pr":[0, 4],
        "Nd":[0, 5],    "Pm":[0, 8],
        "Sm":[0, 9],    "Eu":[0,10],
        "Gd":[0,11],    "Tb":[0,12],
        "Dy":[0,13],    "Ho":[0,14],
        "Er":[0,14],    "Tm":[0,14],
        "Yb":[0,14],    "Lu":[0,14],
        "Hf":[0,10],    "Ta":[0,10],
        "W": [0,10],    "Re":[0,10],
        "Os":[0,10],    "Ir":[0,10],
        "Pt":[0,10],    "Au":[0,10],
        "Hg":[0,10],    "Tl":[0, 6],
        "Pb":[0, 6],    "Bi":[0, 6],
        "Po":[0, 6],    "At":[0, 6],
        "Rn":[0, 6],    "Fr":[0, 2],
        "Ra":[0, 2],    "Ac":[0,2],
        "Th":[0, 2],    "Pa":[0, 3],
        "U": [0, 4],    "Np":[0,14],
        "Pu":[0,14],    "Am":[0,14],
        "Cm":[0,14],    "Bk":[0,14],
        "Cf":[0,14],    "Es":[0,14],
        "Fm":[0,14],    "Md":[0,14],
        "No":[0,14],    "Lr":[0,14],
        "Rf":[0,10],
        }
