import numpy as np
import scipy.linalg as la
import os,sys,inspect, re
import argparse
import itertools
from qutip import *

data_path = '/scratch/midway2/yunong/quantum_optimal_control/pulses/outputpulses/'

from quantum_optimal_control.helper_functions.grape_functions import *
from quantum_optimal_control.main_grape.grape import Grape
from quantum_optimal_control.compilation.parser import *
from quantum_optimal_control.compilation.qcircuit import *
from quantum_optimal_control.compilation.qgate import *
from quantum_optimal_control.compilation.custgate import *
from quantum_optimal_control.compilation.simulator import *

argParser = argparse.ArgumentParser(description='Verify the unitary')
argParser.add_argument('--circ', dest='circ', default='')
argParser.add_argument('--time', dest='time', default='')
argParser.add_argument('--qutip', dest='qut', action='store_true')
argParser.add_argument('--no-qutip', dest='qut', action='store_false')
args = argParser.parse_args()

argParser.set_defaults(qut=False)

def is_unitary(m):
    return np.allclose(np.eye(m.shape[0]), m.H * m)

qasm_file = args.circ
time = args.time
qut = args.qut

qp = QasmParser(qasm_file)
qc = QCircuit(qp.qubits)
qubit_num = len(qc.wires)
qc.add_gate(custQGate(qb.gates))

mapping = qc.mapping
final_U = qc.op_tab[0].unitary()
U2 = Simulator.qutip_unitary(qc.op_tab[0])
print "is your result unitary?", is_unitary(np.matrix(final_U))
print U2
print "is qutip result unitary?", is_unitary(np.matrix(U2))

print "equal?",np.allclose(final_U, U2, atol=1E-3)
# Might have a global phase difference

