from pygrisb.gutz.molecule import gmolecule, gPGAnalyzer
import pprint, logging

logging.basicConfig(level=logging.DEBUG)

mol = gmolecule(["H", "H", "H", "H"], [[0,0,0], [0,1,0.1], [1,1,0.1], [1,0,0]])
#mol = gmolecule(["H", "H", "H", "H"], [[-1,-1,0], [-1,1,0], [1,-1,0], [1,1,0]])

pprint.pprint(mol.species)
pprint.pprint(mol.name)
pprint.pprint(mol.center_of_mass)
pprint.pprint(mol.distance_matrix)

pnalyzer = gPGAnalyzer(mol)
pprint.pprint(pnalyzer.sch_symbol)
pprint.pprint(pnalyzer.non_pure_rotation_generator)
pprint.pprint(pnalyzer.pure_rotations_post_screening)
