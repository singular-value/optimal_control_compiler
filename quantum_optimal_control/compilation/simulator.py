import numpy as np 
from errors import do_error
from qutip import QubitCircuit

# A simple quantum simulator, good enough for 10-20 qubits
# Calculate unitary for GCirc
# Accelarated by putting more gates in one layer

class Simulator:

    def __init(self, qubit_state_num = 2): 
        self.qubit_state_num = qubit_state_num 

    @staticmethod 
    def calculate_unitary(circ):

        qubit_num = len(circ.wires)

        final_U = np.identity(2 ** qubit_num, dtype = np.complex)

        for n_layer, layer in enumerate(circ.layers):

            qstate_list = [[_GateMasterDef(name = 'Id')] * qubit_num]
                
            # In the same layer, gates don't share qubits 

            for g in layer:

                op = circ.op_tab[g]

                if op.name == 'CNOT':

                    qstate_list_ext = [None] * len(qstate_list)

                    for i in range(len(qstate_list)):
                        qstate_list_ext[i] = list(qstate_list[i])

                    ctrl = circ.qb2idx[op.wires[0]] - 1
                    tgt = circ.qb2idx[op.wires[1]] - 1
                        
                    for i in range(len(qstate_list)):

                        qstate_list[i][ctrl] = _GateMasterDef(name = 'P0')
                        qstate_list[i][tgt] = _GateMasterDef(name = 'Id')
                        qstate_list_ext[i][ctrl] = _GateMasterDef(name = 'P1')
                        qstate_list_ext[i][tgt] = _GateMasterDef(name = 'X')

                    qstate_list = qstate_list + qstate_list_ext

                elif op.name == 'CZ':

                    qstate_list_ext = [None]*len(qstate_list)

                    for i in range(len(qstate_list)):
                        qstate_list_ext[i] = list(qstate_list[i])

                    ctrl = circ.qb2idx[op.wires[0]] - 1
                    tgt = circ.qb2idx[op.wires[1]] - 1
                        
                    for i in range(len(qstate_list)):

                        qstate_list[i][ctrl] = _GateMasterDef(name = 'P0')
                        qstate_list[i][tgt] = _GateMasterDef(name = 'Id')
                        qstate_list_ext[i][ctrl] = _GateMasterDef(name = 'P1')
                        qstate_list_ext[i][tgt] = _GateMasterDef(name = 'Z')

                    qstate_list = qstate_list + qstate_list_ext

                elif op.name == 'Rx' or op.name == 'Ry' or op.name == 'Rz':
                    rot_mat = _GateMasterDef(name = op.name, para = op.rot)
                    
                    for i in range(len(qstate_list)):
                        qstate_list[i][circ.qb2idx[op.wires[0]]-1] = rot_mat

                elif op.name == 'SWAP':
                    # Magic trick
                        
                    qstate_list_ext = [None]*len(qstate_list)
                    qstate_list_ext1 = [None]*len(qstate_list)
                    qstate_list_ext2 = [None]*len(qstate_list)

                    for i in range(len(qstate_list)):
                        qstate_list_ext[i] = list(qstate_list[i])
                        qstate_list_ext1[i] = list(qstate_list[i])
                        qstate_list_ext2[i] = list(qstate_list[i])


                    ctrl = circ.qb2idx[op.wires[0]] - 1
                    tgt = circ.qb2idx[op.wires[1]] - 1
                        
                    for i in range(len(qstate_list)):

                        qstate_list[i][ctrl] = _GateMasterDef(name = 'P0')
                        qstate_list[i][tgt] = _GateMasterDef(name = 'P0')
                        qstate_list_ext[i][ctrl] = _GateMasterDef(name = 'P10')
                        qstate_list_ext[i][tgt] = _GateMasterDef(name = 'P01')
                        qstate_list_ext1[i][ctrl] = _GateMasterDef(name = 'P01')
                        qstate_list_ext1[i][tgt] = _GateMasterDef(name = 'P10')
                        qstate_list_ext2[i][ctrl] = _GateMasterDef(name = 'P1')
                        qstate_list_ext2[i][tgt] = _GateMasterDef(name = 'P1')

                    qstate_list = qstate_list + qstate_list_ext + qstate_list_ext1 + qstate_list_ext2

                else:

                    mat = _GateMasterDef(name = op.name)
                    for i in range(len(qstate_list)):
                        qstate_list[i][circ.qb2idx[op.wires[0]]-1] = mat

            n_layer += 1 

            crt = np.zeros([2 ** qubit_num, 2 ** qubit_num])

            for state in qstate_list:
                crt = crt + _kron_list(state)

            final_U = np.dot(crt, final_U)  

        return final_U

    @staticmethod
    def qutip_unitary(c_gate):
        
        circ = c_gate.gcirc

        qc = QubitCircuit(len(circ.wires))
        for gate in circ.op_tab: 
            if gate.name == "H":
                tgt = circ.qb2idx[gate.wires[0]]-1
                qc.add_gate("RX", targets=tgt, arg_value =np.pi/2)
                qc.add_gate("RZ", targets=tgt, arg_value =np.pi/2)
                qc.add_gate("RX", targets=tgt, arg_value =np.pi/2)

            if gate.name == "T":
                tgt = circ.qb2idx[gate.wires[0]]-1
                qc.add_gate("RZ", targets=tgt, arg_value =np.pi/8)
            if gate.name == "Tdag":
                tgt = circ.qb2idx[gate.wires[0]]-1
                qc.add_gate("RZ", targets=tgt, arg_value =-np.pi/8)

            if gate.name == "S":
                tgt = circ.qb2idx[gate.wires[0]]-1
                qc.add_gate("RZ", targets=tgt, arg_value =np.pi/4)

            if gate.name == "X":
                tgt = circ.qb2idx[gate.wires[0]]-1
                qc.add_gate("RX", targets=tgt, arg_value =np.pi)
            if gate.name == "Y":
                tgt = circ.qb2idx[gate.wires[0]]-1
                qc.add_gate("RY", targets=tgt, arg_value =np.pi)
            if gate.name == "Z":
                tgt = circ.qb2idx[gate.wires[0]]-1
                qc.add_gate("RZ", targets=tgt, arg_value =np.pi)
            if gate.name == "Rz":
                tgt = circ.qb2idx[gate.wires[0]]-1
                ang = gate.rot
                qc.add_gate("RZ", targets=tgt, arg_value = ang)
            if gate.name == "Rx":
                tgt = circ.qb2idx[gate.wires[0]]-1
                ang = gate.rot
                qc.add_gate("RX", targets=tgt, arg_value = ang)
            if gate.name == "CNOT":
                ctrl = circ.qb2idx[gate.wires[0]]-1
                tgt = circ.qb2idx[gate.wires[1]]-1
                qc.add_gate("CNOT",controls=ctrl, targets=tgt)

        ulist = qc.propagators()
        u = gate_sequence_product(ulist)
        return u.full()
        
def _kron_list(args):
    ret = args[0]
    for item in args[1:]:
        ret = np.kron(ret, item)
    return ret

# Was a constant, keep the naming anyway
def _GateMasterDef(name = '', para = None):

    if name == 'H':
        return 1./np.sqrt(2) * np.array([[1.0,1.0],[1.0,-1.0]], dtype = np.complex)
    if name == 'X':
        return np.array([[0.0, 1.0],[1.0,0.0]], dtype = np.complex)
    if name == 'Y':
        return np.array([[0.0, -1.0j],[1.0j,0.0]], dtype = np.complex)
    if name == 'CNOT': 
        return np.array([[1.0,0.0,0.0, 0.0],[0.0,1.0,0.0, 0.0],[0.0,0.0,0.0, 1.0],[0.0,0.0,1.0, 0.0]], dtype = np.complex)
    if name == 'CZ': 
        return np.array([[1.0,0.0,0.0, 0.0],[0.0,1.0,0.0, 0.0],[0.0,0.0,1.0, 0.0],[0.0,0.0,0.0, -1.0]], dtype = np.complex)
    if name == 'Z':
        return np.array([[1.0, 0.0],[0.0,-1.0]], dtype = np.complex)
    if name == 'T':
        return np.array([[1.0, 0.0],[0.0,np.exp(1j*np.pi/4.0)]], dtype = np.complex)
    if name == 'S':
        return np.array([[1.0, 0.0],[0.0,np.exp(1j*np.pi/2.0)]], dtype = np.complex)
    if name == 'Sdag':
        return np.array([[1.0, 0.0],[0.0,-np.exp(1j*np.pi/2.0)]], dtype = np.complex)
    if name == 'Tdag':
        return np.array([[1.0, 0.0],[0.0,-np.exp(1j*np.pi/4.0)]], dtype = np.complex)
    if name == 'Rz':
        return np.array([[np.exp(-1j * para / 2), 0],[0, np.exp(1j * para / 2)]], dtype = np.complex)
    if name == 'Rx':
        return np.array([[np.cos(para/2), -1j * np.sin(para / 2)], [-1j * np.sin(para / 2), np.cos(para / 2)]], dtype = np.complex)
    if name == 'Ry':
        return np.array([[np.cos(para / 2), - np.sin(para / 2)], [np.sin(para / 2), np.cos(para / 2)]], dtype = np.complex)
    if name == 'P0':
        return np.array([[1.0, 0.0], [0.0,0.0]], dtype = np.complex) 
    if name == 'P1':
        return np.array([[0.0, 0.0], [0.0,1.0]], dtype = np.complex)
    if name == 'P01':
        return np.array([[0.0, 0.0], [1.0,0.0]], dtype = np.complex)
    if name == 'P10':
        return np.array([[0.0, 1.0], [0.0,0.0]], dtype = np.complex)
    if name == 'Id':
        return np.identity(2)

    return None
