from quantum_optimal_control.compilation.parser import QasmParser
from quantum_optimal_control.compilation.qgate import *
from quantum_optimal_control.compilation.gcirc import *
from quantum_optimal_control.compilation.qcircuit import *
from quantum_optimal_control.compilation.simulator import *
# When compile for the first time, draw_circuit can be slow
from quantum_optimal_control.compilation.visualization import draw_circuit
from quantum_optimal_control.compilation.custgate import *

import time

qp = QasmParser("../qasm/ising_model.n10.qasm")
# qp = QasmParser("../qasm/simple.qasm")
# qp = QasmParser("../qasm/qaoa20.qasm")
# qp = QasmParser("../qasm/qaoa_line.qasm")
# qp = QasmParser("../qasm/qaoa_two_step.qasm")
qc = QCirc(qp.qubits)

for g in qp.gates:
    qc.add_gate(custQGate([g]))

qc.detect_commutation()
    
draw_circuit(qc, "com_qc")

print "now in cls"
qc.cls(filename="output.qasm")
print "circuit length:", qc.total_length()

qc.simple_optimization()

qc.diagonal_merge()
qc.detect_commutation()

# scheduled = qc.cls(filename="output_merged.qasm")
draw_circuit(qc, "d_merged_qc")
print "circuit length:", qc.total_length()
