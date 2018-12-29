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
from quantum_optimal_control.compilation.hamiltonian import *
from quantum_optimal_control.compilation.gen_pulse import *

argParser = argparse.ArgumentParser(description='generate optimal pulses')
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

qp = QasmParser(qasm_file) # Here assume a qasm_file that output by inject_swap.py, so qc.mapping is valid
#print qp.qubits
qc = QCirc(qp.qubits)
qubit_num = len(qc.wires)
g_list = []
for g in qp.gates:
    g_list.append(g)

cg = custQGate(g_list)
qc.add_gate(cg)

mapping = qc.mapping
final_U = cg.unitary()

run_circuit(qasm_file, qubit_num, final_U, mapping, time, if_qutip =qut )
