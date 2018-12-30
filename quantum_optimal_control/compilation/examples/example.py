from quantum_optimal_control.compilation.parser import QasmParser
from quantum_optimal_control.compilation.qgate import *
from quantum_optimal_control.compilation.gcirc import *
from quantum_optimal_control.compilation.qcircuit import *
from quantum_optimal_control.compilation.simulator import *
# When compile for the first time, draw_circuit can be slow
from quantum_optimal_control.compilation.visualization import draw_circuit
from quantum_optimal_control.compilation.custgate import *

import time

# qp = QasmParser("../qasm/ising_model.n10.qasm")
# qp = QasmParser("../qasm/simple2.qasm")
# qp = QasmParser("../qasm/simple.qasm")
# qp = QasmParser("../qasm/qaoa20.qasm")
# qp = QasmParser("../qasm/qaoa_line.qasm")
# qp = QasmParser("../qasm/qaoa_two_step.qasm")
# qp = QasmParser("../qasm/sr2.qasm")
qp =QasmParser("../qasm/square_root.n02.qasmf")
qc = QCirc(qp.qubits)

for g in qp.gates:
    qc.add_gate(custQGate([g]))

qc.detect_commutation()

for key in qc.commutation_set.keys():
    if isinstance(key, str):
        print "key wire", key
        for c_set in qc.commutation_set[key]:
            print "c_set"
            for c_gate in c_set:
                print "c_gate"
                for gq in c_gate.gate_list:
                    print qg.name, qg.wires
    else:
        print key[1]
        for qg in key[0].gate_list:
            print qg.name, qg.wires

# draw_circuit(qc, "com_qc")
qc.detect_commutation()
qc.simple_optimization()

scheduled = qc.cls(filename="output_merged.qasm")

print "666"
for wire in qc.op_tab[6].wires:
    for ccc in qc.op_tab[6].suc[wire]:
        print "suc"
        for qg in ccc.gate_list:
            print qg.name, qg.wires

    for ccc in qc.op_tab[6].pred[wire]:
        print "pred"
        for qg in ccc.gate_list:
            print qg.name, qg.wires

print "set pred"
for ccc in qc.pred_set[qc.op_tab[6]]:
    print "set pred"
    for qg in ccc.gate_list:
        print qg.name, qg.wires

print "set suc"
for ccc in qc.suc_set[qc.op_tab[6]]:
    print "set suc"
    for qg in ccc.gate_list:
        print qg.name, qg.wires

print "after cls"
for merged_gate in qc.op_tab:
    print "gate", qc.op_tab.index(merged_gate)
    for qg in merged_gate.gate_list:
        print qg.name, qg.wires

    print "his pred"
    for cc in qc.pred_set[merged_gate]:
        print "pred", qc.op_tab.index(cc)
        for qg in cc.gate_list:
            print qg.name, qg.wires

    print "his suc"
    for cc in qc.suc_set[merged_gate]:
        print "suc", qc.op_tab.index(cc)
        for qg in cc.gate_list:
            print qg.name, qg.wires

qc.diagonal_merge()

qc.detect_commutation()
draw_circuit(qc, "d_merge_qc")
print "after D_MERGE"
for merged_gate in qc.op_tab:
    print "gate", qc.op_tab.index(merged_gate)
    for qg in merged_gate.gate_list:
        print qg.name, qg.wires

    print "his pred"
    for cc in qc.pred_set[merged_gate]:
        print "pred", qc.op_tab.index(cc)
        for qg in cc.gate_list:
            print qg.name, qg.wires

    print "his suc"
    for cc in qc.suc_set[merged_gate]:
        print "suc", qc.op_tab.index(cc)
        for qg in cc.gate_list:
            print qg.name, qg.wires

# draw_circuit(qc, "d_merge_qc")
qc.simple_optimization()
scheduled = qc.cls(filename="output_merged.qasm")
draw_circuit(qc, "cls_merge_qc")
qc.detect_commutation()
# draw_circuit(qc, "d_merged_qc")

for c_gate in qc.op_tab:
    print c_gate
qc.block_merge(7)
draw_circuit(qc, "block_merged_qc")
