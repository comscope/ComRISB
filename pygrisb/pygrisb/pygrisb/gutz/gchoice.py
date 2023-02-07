import sys, h5py, logging, inquirer, numpy, pprint, json, warnings
from pygrisb.mbody.coulomb_matrix import coulomb_matrix_slater


class gchoice:
    '''class for aquiring user choices to perform gutzwiller local correlation
    simulations.
    '''
    def __init__(self, symbols, default_idx_equivalent_atoms):
        self.symbols = symbols
        self.default_idx_equivalent_atoms = default_idx_equivalent_atoms
        # core parameters dict
        self.ps = {}
        self.run()

    @property
    def params(self):
        return self.ps

    def set_default_gchoice(self):
        if self.default_idx_equivalent_atoms is None:
            self.default_idx_equivalent_atoms = list(range(len(self.symbols)))
        # unit used in the gutzwiller calculation.
        if '-u' in sys.argv:
            self.ps["unit"] = sys.argv[sys.argv.index('-u') + 1]
            logging.info(f"unit is set to {self.ps['unit']}.")
        else:
            self.ps["unit"] = 'rydberg'
        # spin up/dn order convention
        if '-s' in sys.argv:
            val = int(sys.argv[sys.argv.index('-s') + 1])
            val = 1 if val > 0 else -1
            self.ps["spin_order"] = val
            orders = {1:"dnup", -1:"updn"}
            logging.info(f"spin order convention is set to {orders[val]}.")
        else:
            self.ps["spin_order"] = -1

        self.ps.update({
                "ispin": 1,
                "orbital_polarization": 1,
                "iso": 1,
                "u_matrix_type": 1,
                "dc_mode": 12,
                "idx_equivalent_atoms": None,
                "unique_cor_symbol_list": [],
                "unique_u_list": [],
                "unique_j_list": [],
                "unique_fs_list": [],
                "unique_cor_orb_list": [],
                "unique_nelf_list": [],
                "giembeddiag": -1})

    def save_params(self):
        raise NotImplementedError("save_params not implemented.")

    def load_params(self):
        raise NotImplementedError("load_params not implemented.")

    def print_info(self):
            pprint.pprint(self.ps)

    def get_input(self):
        '''dummy function, to be overwritten.
        '''
        raise NotImplementedError("get_input not implemented.")

    def run(self):
        '''main driver function.
        '''
        try:
            self.load_params()
            self.update_fs()
        except:
            self.set_default_gchoice()
            self.get_input()
            self.save_params()
        # Sanity check
        if self.ps["giembeddiag"] == 10:
            if self.ps["dc_mode"] != 1:
                warnings.warn("hf type calculation: dc_mode reset to 1!")
                self.ps["dc_mode"] = 1

    def update_fs(self):
        if "unique_fs_list" not in self.ps:
            self.ps["unique_fs_list"] = []
            for cor_orb, u, j in zip(self.ps["unique_cor_orb_list"], \
                    self.ps["unique_u_list"], self.ps["unique_j_list"]):
                fs = coulomb_matrix_slater.get_slater_integrals(u, j, \
                        orb=cor_orb).tolist()
                self.ps["unique_fs_list"].append(fs)


class h5_gchoice(gchoice):
    def __init__(self, symbols, default_idx_equivalent_atoms=None,
            fname="ginit.h5"):
        self.fname = fname
        super().__init__(symbols, default_idx_equivalent_atoms)

    def save_params(self):
        '''save user inputs to hdf5 file to be referenced later.
        '''
        with h5py.File(self.fname, 'a') as f:
            if '/gchoice' in f:
                del f['/gchoice']
            for key, val in self.ps.items():
                if val is not None:
                    # str to bytes conversion
                    if type(val) is list and type(val[0]) is str:
                        val = [v.encode() for v in val]
                    f[f'/gchoice/{key}'] = val

    def load_params(self):
        '''load inputs saved previously.
        not very good using hdf5 here, as it is designed to store large arrays.
        '''
        with h5py.File(self.fname, 'r') as f:
            for key, val in f["/gchoice"].items():
                val = val[()]
                # bytes to str conversion
                if type(val) is numpy.ndarray and type(val[0]) is numpy.bytes_:
                    val = [v.decode() for v in val]
                self.ps[key] = val


class json_gchoice(gchoice):
    def __init__(self, symbols, default_idx_equivalent_atoms=None,
            fname="ginit.json"):
        self.fname = fname
        super().__init__(symbols, default_idx_equivalent_atoms)

    def save_params(self):
        '''save user inputs to json file to be referenced later.
        '''
        data = json.load(open(self.fname, "r"))
        data["gchoice"] = self.ps
        with open(self.fname, "w") as f:
            json.dump(data, f, indent = 4)

    def load_params(self):
        '''load inputs saved previously.
        '''
        data = json.load(open(self.fname, "r"))
        self.ps = data["gchoice"]


class gchoice_screen(json_gchoice):
    '''class to handle a list of questions for user input over the screen,
    to guide the set up of gutzwiller calculation.
    final results can also be stored in hdf5 file
    to be referenced later.
    '''
    import inquirer
    def get_input(self):
        '''question and answer session over the screen.
        '''
        print("\nUser inputs to initialize G-RISB simulation.")
        questions = [
                inquirer.List('ispin',
                        message="Break spin-symmetry",
                        choices=[("no", 1), ("yes", 2)],),

                inquirer.List('orbital_polarization',
                        message="Break orbital-symmetry",
                        choices=[
                                ("no", 0),
                                ("crystal field effect", 1), \
                                ("full symmetry breaking", 2)],),

                inquirer.List('iso',
                        message="Include spin-orbit coupling",
                        choices=[("no", 1), ("yes", 2)],),

                inquirer.List('u_matrix_type',
                        message="Parametrize Coulomb-matrix",
                        choices=[
                                ("Slater-Condo with [U,J]", 1),
                                ("Slater-Condo with [F0,F2,...]", 3),
                                ("Kanamori with [U,J]", 2),
                                ("Manual input", 0)],),

                inquirer.List('dc_mode',
                        message="Coulomb double counting",
                        choices=[
                                ("FLL dc (updated at charge iter., Rec.)", 12),
                                ("Fix dc potential", 2),
                                ("FLL dc self-consistent", 1),
                                ("No dc", 0)],),

                inquirer.List('giembeddiag',
                        message="Solution of embedding Hamiltonian",
                        choices=[
                                ("Valence truncation ED (VTED)", -1),
                                ("VTED with Sz symmetry", -2),
                                ("VTED with S=0", -3),
                                ("VTED with Jz symmetry", -4),
                                ("ML (kernel-ridge)", -10),
                                ("ML (normal-mode-1d)", -11),
                                ("DMGR (expt.)", -12),
                                ("ML (normal-mode-2d)", -13),
                                ("ML (normal-mode-3d)", -14),
                                ("HF (debugging only)", 10)],),]

        ans = inquirer.prompt(questions)
        self.ps.update(ans)

        # Equivalent atom indices
        print(("Equivalent atom indices:\n"
                "    [0 0 0 1 1] means 0-2 and 3-4 are two sets of eq. atms."))
        question = [
                inquirer.List('idx_equivalent_atoms',
                        message="Equivalent atom indices",
                        choices=[
                                self.default_idx_equivalent_atoms,
                                "modify"],),]
        ans = inquirer.prompt(question)

        if ans['idx_equivalent_atoms'] == "modify":
            def indices_accepted(_, x):
                try:
                    x = list(map(int, x.split()))
                    return len(x) == len(self.symbols) and \
                            0 <= min(x) <= max(x) < len(self.symbols)
                except:
                    return False

            question = [
                    inquirer.Text('idx_equivalent_atoms',
                            message="Enter indices (number + space only)",
                            validate=indices_accepted),]
            ans = inquirer.prompt(question)
            ans['idx_equivalent_atoms'] = \
                    list(map(int, ans['idx_equivalent_atoms'].split()))
        self.ps.update(ans)

        # atom-specific information:
        for i, s in enumerate(self.symbols):
            if s in self.symbols[:i]:
                continue
            print('\n ' + '-'*12 + f"\n atom {i} {s}")
            question = [
                    inquirer.Confirm('cor',
                            message="Is this atom correlated",
                            default=True),]
            ans = inquirer.prompt(question)
            if not ans['cor']:
                continue

            self.ps["unique_cor_symbol_list"].append(s)

            question = [
                    inquirer.List('cor_orb',
                            message="Correlated orbital",
                            choices=['s', 'p', 'd', 'f'],),]
            ans = inquirer.prompt(question)
            cor_orb = ans["cor_orb"]
            self.ps["unique_cor_orb_list"].append(cor_orb)

            if self.ps["u_matrix_type"] in [1, 2]:
                def uj_accepted(_, x):
                    try:
                        x = list(map(float, x.split()))
                        return len(x) == 2
                    except:
                        return False

                question = [
                        inquirer.Text('uj',
                                message="Enter U J(sep. by space, eV)",
                                validate=uj_accepted,),]
                ans = inquirer.prompt(question)
                uj = list(map(float, ans["uj"].split()))
                fs = coulomb_matrix_slater.get_slater_integrals(*uj,
                        orb=cor_orb).tolist()
            elif self.ps["u_matrix_type"] == 3:
                fentry = {'s':"F0", 'p':"F0 F2", 'd':"F0 F2 F4",
                        'f':"F0 F2 F4 F6"}
                def fs_accepted(_, x):
                    try:
                        x = list(map(float, x.split()))
                        return len(x) == len(fentry[cor_orb].split())
                    except:
                        return False
                question = [
                        inquirer.Text('fs',
                                message=f'Enter {fentry[cor_orb]}'+\
                                        '(sep. by space, eV)',
                                validate=fs_accepted,),]
                ans = inquirer.prompt(question)
                fs = list(map(float, ans["fs"].split()))
                uj = coulomb_matrix_slater.get_uj_params(fs, orb=cor_orb)

            if self.ps["u_matrix_type"] > 0:
                self.ps["unique_u_list"].append(uj[0])
                self.ps["unique_j_list"].append(uj[1])
                self.ps["unique_fs_list"].append(fs)

            if self.ps["dc_mode"] == 2:
                def float_accepted(_, x):
                    try:
                        x = float(x)
                        return x >= 0
                    except:
                        return False
                question = [
                    inquirer.Text('nf',
                            message=f"Enter {cor_orb} electron number",
                            validate=float_accepted,),]
                ans = inquirer.prompt(question)
                nf = float(ans["nf"])
                self.ps["unique_nelf_list"].append([nf/2.,nf/2.])
