import numpy, logging, itertools, warnings
from pymatgen.core.structure import Molecule
from pymatgen.symmetry.analyzer import PointGroupAnalyzer


class gmolecule(Molecule):
    '''Adding equivalent atomm indices to Molecule.
    The center of mass is also fixed at [0,0,0],
    which is the location of the atom whose local symmetry is to be
    analysized.
    '''
    def __init__(self, species, coords, eq_indices=None, mol_id=0, **kwargs):
        super().__init__(species, coords, **kwargs)
        # equivalent atomm indices
        self.eq_indices = eq_indices
        # molecule id to avoid duplicated file names.
        self._id = mol_id
        # get succinct molecular name
        self.set_name()

    @property
    def center_of_mass(self):
        # fix center at [0,0,0].
        return numpy.zeros(3)

    @property
    def name(self):
        return self._name

    def set_name(self):
        # beautify name by removing psaces and 1.
        self._name = self.formula.replace(" ", "").replace("1","")

    def save_xyz(self):
        # save molecule structure to simple xyz file,
        # for the convenience of 3d view.
        self.to("xyz", f"{self.name}_{self._id}.xyz")


class gpoint_group_analyzer(PointGroupAnalyzer):
    '''extended class based on PointGroupAnalyzer,
    with convenient access to pure rotations and inversion symmetry.
    it further screens the rotations using the additional
    equivalent atom indices.
    '''
    def __init__(self, mol):
        super().__init__(mol)
        # update pure rotations with screening
        # based on equivalent atom indices.
        self.set_symm_ops()

    @property
    def n_sym_ops_post_screening(self):
        return self._n_symm_ops_post_screening

    @property
    def sym_ops_post_screening(self):
        return self._symm_ops_post_screening

    def set_symm_ops(self):
        # get full point group symmetry operations.
        symmops = self.get_symmetry_operations()
        self._symm_ops = [o.rotation_matrix for o in symmops]
        self._n_symm_ops = len(self._symm_ops)
        # perform screening based on equivalent atom indices.
        self._symm_ops_post_screening = [\
                op for op in self._symm_ops \
                if self.is_compatible(self.mol, op)]
        self._n_symm_ops_post_screening = \
                len(self._symm_ops_post_screening)
        # make sure of the integer ratio
        print(f"n_symm_ops = {self._n_symm_ops} from point-group analysis\n" +\
                + 11*" " \
                + f"= {self._n_symm_ops_post_screening} after screening.")
        assert(self._n_symm_ops % self._n_symm_ops_post_screening == 0), \
                "sym op screening error."

    @staticmethod
    def is_compatible(mol, rotation, tol=0.01):
        # screening based on equivalent atom indices.
        for i, coord in enumerate(mol.cart_coords):
            coordp = rotation.dot(coord)
            diff = mol.cart_coords - coordp
            inds = numpy.where(numpy.all(numpy.abs(diff) < tol, axis=1))[0]
            # make sure rotation can map each atom to
            # one and only one of itself or images.
            if len(inds) == 0:
                warnings.warn("cannot locate symmetry-mapped atom.")
                return False
            assert(len(inds) == 1), "overlap atoms identified!"
            if mol.eq_indices is not None and \
                    mol.eq_indices[i] != mol.eq_indices[inds[0]]:
                return False
        return True


class cell_point_groups_analyzer:
    '''get point group operations for selected atoms in a periodic structure.
    '''
    def __init__(self, atoms):
        # atoms object
        self.atoms = atoms
        # selected atom indices for point group analysis.
        self.evaluate()

    @property
    def symm_ops_list(self):
        return self._symm_ops_list

    @staticmethod
    def get_local_atomic_environ(atoms, iat, nmax=10):
        '''get the local environment of atom iat as a molecule
        centered at atom iat.
        '''
        scaled_positions = numpy.copy(atoms.scaled_positions)
        scaled_positions -= scaled_positions[iat]
        n = (2*nmax)**3*len(scaled_positions)
        molecule_positions = numpy.empty((n,3))
        molecule_symbols = numpy.empty(n, dtype='<U2')
        pair_dist = numpy.empty(n)
        eq_indices = numpy.empty(n, dtype=numpy.int)
        ii = 0
        for i, j, k in itertools.product(range(-nmax, nmax), repeat=3):
            for jat, sp in enumerate(scaled_positions):
                n_sp = sp + numpy.array([i, j, k])
                n_p = n_sp.dot(atoms.cell)
                pair_dist[ii] = numpy.linalg.norm(n_p)
                molecule_positions[ii,:] = n_p
                molecule_symbols[ii] = atoms.symbols[jat]
                eq_indices[ii] = atoms.ps["idx_equivalent_atoms"][jat]
                ii += 1
        # Get reasonable dist_cut
        idx_sort = pair_dist.argsort()
        pair_dist = pair_dist[idx_sort]
        # 12 is a good spiritual number
        dist_cut = pair_dist[12] + 2.0
        logging.info("Default dist_cut to extract a centered cluster" +\
                f"\n for symmetry evaluation = {dist_cut}")
        idx_local = pair_dist < dist_cut
        molecule_positions = molecule_positions[idx_sort][idx_local]
        molecule_symbols = molecule_symbols[idx_sort][idx_local]
        eq_indices = eq_indices[idx_sort][idx_local]
        # apply local rotation if required.
        if atoms.locrot_list is not None:
            molecule_positions = molecule_positions.dot(atoms.locrot_list[iat])
        return gmolecule(molecule_symbols, molecule_positions,
                eq_indices=eq_indices, mol_id=iat)

    def evaluate(self):
        symm_ops_list = []
        for i in self.atoms.cor_atm_list:
            mol = self.get_local_atomic_environ(self.atoms, i)
            gpg = gpoint_group_analyzer(mol)
            print(f"correlated atom {i} with point group: {gpg.sch_symbol}.")
            symm_ops_list.append(gpg.sym_ops_post_screening)
        self._symm_ops_list = symm_ops_list
