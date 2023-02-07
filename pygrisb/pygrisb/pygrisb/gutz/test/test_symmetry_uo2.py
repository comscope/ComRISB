from phonopy.structure.symmetry import Symmetry
from phonopy.interface.wien2k import parse_wien2k_struct
from phonopy.units import Bohr
from phonopy.interface.vasp import write_vasp


def clean_scaled_positions(cell):
    positions = cell.get_scaled_positions()
    for pos in positions:
        for i in (0, 1, 2):
            # The following %19.16f follows write_vasp
            if float("%19.16f" % pos[i]) >= 1:
                pos[i] -= 1.0
    cell.set_scaled_positions(positions)



cell, npts, r0s, rmts = parse_wien2k_struct("UO2.struct")
sym = Symmetry(cell)
print(sym.get_pointgroup())
print(sym.get_site_symmetry(0))
print(sym._dataset["site_symmetry_symbols"])


positions = cell.get_scaled_positions()
lattice = cell.get_cell() * Bohr
cell.set_cell(lattice)
cell.set_scaled_positions(positions)
clean_scaled_positions(cell)
#write_vasp("POSCAR.wien2k", cell, direct=True)

