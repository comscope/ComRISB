import numpy
from qiskit.chemistry import FermionicOperator
from qiskit.aqua.operators.legacy.op_converter import to_weighted_pauli_operator
from qiskit.chemistry.components.variational_forms import UCCSD
from qiskit.aqua.components.optimizers import L_BFGS_B
from qiskit import Aer
from qiskit.quantum_info import Pauli
from qiskit.aqua.operators import WeightedPauliOperator
from qiskit.aqua.operators.legacy import op_converter
from qiskit.aqua.algorithms import VQE
from qiskit.aqua import QuantumInstance
import itertools
from qiskit import QuantumRegister, QuantumCircuit, execute
from qiskit.aqua.components.initial_states import Custom
from qiskit.chemistry.components.initial_states import HartreeFock
import scipy


#Constructing the connectivity matrices for h1 one body Hamiltonian
with open('v1e.dat','r') as f:
    lines=f.readlines()[1:]
    num_sites=4
    eg_h1=numpy.zeros((2*num_sites,2*num_sites))
    for line in lines:
        elems=line.split()
        eg_h1[int(elems[0])][int(elems[1])]=float(elems[2])
        eg_h1[int(elems[0])+num_sites][int(elems[1])+num_sites]=float(elems[2])

with open('v2e.dat','r') as f:
    num_sites=4
    eg_h2=numpy.zeros((2*num_sites,2*num_sites,2*num_sites,2*num_sites))
    for line in f:
        if "#" in line:
            continue
        line = line.split()
        i,j,k,l = map(int, line[:4])
        val = float(line[4])
        eg_h2[i,j,k,l] \
                = eg_h2[i+num_sites,j+num_sites,k,l] \
                = eg_h2[i,j,k+num_sites,l+num_sites] \
                = eg_h2[i+num_sites,j+num_sites,k+num_sites,l+num_sites] \
                = 0.5*val  # convention with 0.5 factor included.

def qubit_Hamiltonian_eg(h1,h2):
    qubit_op=FermionicOperator(h1=h1,h2=h2).mapping('jordan_wigner')
    return qubit_op

#qubit operator rep for the eg Hamiltonian
qubit_eg = qubit_Hamiltonian_eg(eg_h1,eg_h2)
eg_mat = op_converter.to_matrix_operator(qubit_eg).dense_matrix
w, v = numpy.linalg.eigh(eg_mat)
print(w)
