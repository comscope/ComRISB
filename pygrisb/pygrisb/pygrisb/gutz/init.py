import argparse
import pygrisb.gutz.atoms as gatms
import pygrisb.gutz.gparam as gp
import pygrisb.gutz.gchoice as gc


class init:
    '''base class to handle gutzwiller initialization
    with multiple interfaces.
    '''
    def run(self):
        self.get_inline_args()
        self.set_atoms()
        self.set_gchoice()
        self.atoms.set_all_attributes(
                fixsab=self.args.fixsab,
                realhemb=self.args.realhemb,
                valtrunc=self.args.valtrunc,
                )
        self.set_gparam_atms()

    def get_inline_args(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("--fixsab", action="store_true",
                help="fix the symm adapted basis during reinitialization.")
        parser.add_argument("--realhemb", action="store_true",
                help="real embedding hamiltonian.")
        parser.add_argument("--valtrunc", action="store_true",
                help="valence truncation for ed.")
        parser.add_argument("-u", type=str, default="rydberg",
                help="energy unit in CyGutz calculation [rydberg(dflt), ev].")
        parser.add_argument("-s", type=int, default=-1,
                help="spin convention [-1:updn(dflt); 1: dnup]")
        self.args, unknown = parser.parse_known_args()

    def set_gchoice(self):
        raise NotImplementedError("function evaluate not implemented!")

    def set_atoms(self):
        self.atoms = gatms.json_cor_structure()

    def set_gparam_atms(self):
        '''save Gparam_atms.h5.
        '''
        gpa = gp.gparam_atms(self.atoms.ps["na2_list"], self.atoms.ps)
        gpa.h5save()


class screen_init(init):
    '''initialization of ga with inputs from screen.
    '''
    def set_gchoice(self):
        g_choice = gc.gchoice_screen(self.atoms.symbols, \
                default_idx_equivalent_atoms= \
                self.atoms.default_idx_equivalent_atoms)
        self.atoms.update(g_choice.params)


if __name__=="__main__":
    s_init = screen_init()
    s_init.run()
