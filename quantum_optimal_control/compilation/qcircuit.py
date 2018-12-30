"""
- This module provides the class for the quantum circuit of customized gates
- Not assume to be simulable
- Commutation relations are regarded
"""

from collections import defaultdict
from itertools import compress
import numpy as np
import networkx as nx
import copy
import sys
from quantum_optimal_control.compilation.simulator import Simulator
from quantum_optimal_control.compilation.qgate import QGate
from quantum_optimal_control.compilation.custgate import custQGate
from quantum_optimal_control.compilation.gcirc import GCirc
# from quantum_optimal_control.compilation.errors import do_error
# from quantum_optimal_control.compilation.visualization import draw_circuit
import time

# For large graphs
sys.setrecursionlimit(20000) 

class QCirc(object):
    """
    - This class provides the quantum circuit structure
    """

    def __init__(self, wires, mapped = False):

        self.wires = wires
        self.qb_tab = {}  # keep track of the gates on each wire

        self.commutation_set = {}

        for wire in self.wires:
            self.qb_tab[wire] = []
            self.commutation_set[wire] = []

        self.op_tab = []
        self.execution_time = {} # execute time for each gate
        self.mapping = []
        self.critical_paths = []
        self.paths = defaultdict(lambda : [])

        self.pred_set = defaultdict(lambda : [])
        self.suc_set = defaultdict(lambda : [])
        self.com_set = defaultdict(lambda : [])

        if mapped:
            for wire in self.wires:
                re_bit = re.compile('qubit_(\d)_(\d_').search(wire)
                self.mapping.append((int(re_bit.group(1)), int(re_bit.group(2))))

    def add_gate(self, cgate):
        """
        - Add gate to the circuit, handle the QGate relations
        """

        self.op_tab.append(cgate) # Put the gate into the gate list

        for wire in cgate.wires:
            if len(self.qb_tab[wire]) != 0:
                self.qb_tab[wire][-1].add_suc(wire, cgate)
                cgate.add_pred(wire, self.qb_tab[wire][-1])
            self.qb_tab[wire].append(cgate)

    def simple_optimization(self):
        """
        - Repeatedly cancel self-inverse gates through commutation relations
        """

        pre_len = 0
        cur_len = len(self.op_tab)

        while cur_len != pre_len:

            pre_len = cur_len
            self.detect_commutation()
            self._simple_cancellation()
            cur_len = len(self.op_tab)

        return

    def _simple_cancellation(self):
        """
        - Execute after detecting commutation relations, but before merging.
        """

        # pylint: disable=too-many-locals
        # pylint: disable=too-many-nested-blocks
        # pylint: disable=too-many-branches

        cancel_set = defaultdict(lambda: [])

        for wire in self.wires:
            for n_c, c_set in enumerate(self.commutation_set[wire]):
                for gate in c_set:
                    if len(gate.wires) == 1:

                        cancel_set[(n_c, gate.gate_list[0].name, wire)].append(gate)

                    # - Not worry about cz and cy for now
                    if gate.gate_list[0].name == "CNOT" and gate.wires[0] == wire:

                        cnot_key = (n_c, wire,
                                    self.commutation_set[(gate, gate.wires[1])],
                                    gate.wires[1])
                        cancel_set[cnot_key].append(gate)

        for key in cancel_set:
            if len(cancel_set[key]) > 1:
                if key[1] == 'Rz' or key[1] == 'Rx' or key[1] == 'Ry':

                    op_pos = self.op_tab.index(cancel_set[key][0])
                    qb_pos = self.qb_tab[key[2]].index(cancel_set[key][0])

                    rot = 0.0

                    for gate_to_rm in cancel_set[key]:

                        rot += gate_to_rm.gate_list[0].rot
                        self.op_tab.remove(gate_to_rm)


                        self.qb_tab[key[2]].remove(gate_to_rm)

                        for suc_gate in gate_to_rm.suc[key[2]]:

                            suc_gate.remove_pred(key[2], gate_to_rm)

                        for pre_gate in gate_to_rm.pred[key[2]]:
                            pre_gate.remove_suc(key[2], gate_to_rm)

                        self.commutation_set[key[2]][key[0]].remove(gate_to_rm)

                    qgate = QGate(key[1], [key[2]], rot=rot)
                    cgate_to_add = custQGate([qgate])
                    self.commutation_set[key[2]][key[0]].append(cgate_to_add)
                    if key[0] >= 1:
                        for pred in self.commutation_set[key[2]][key[0]-1]:
                            pred.add_suc(key[0], cgate_to_add)

                    if key[0] < len(self.commutation_set[key[2]]) - 1:
                        for suc in self.commutation_set[key[2]][key[0] + 1]:
                            suc.add_pred(key[0], cgate_to_add)
                    self.op_tab.insert(op_pos, cgate_to_add)
                    self.qb_tab[key[2]].insert(qb_pos, cgate_to_add)

                if key[1] in ('Z', 'X', 'Y', 'H') or len(key) == 4:
                    for gate_to_rm in cancel_set[key][:(len(cancel_set[key])/2*2)]:

                        self.op_tab.remove(gate_to_rm)

                        for wire in gate_to_rm.wires:

                            self.qb_tab[wire].remove(gate_to_rm)

                            for suc_gate in gate_to_rm.suc[wire]:

                                suc_gate.remove_pred(wire, gate_to_rm)

                                if len(suc_gate.pred[wire]) == 0:

                                    for pre in gate_to_rm.pred[wire]:

                                        suc_gate.add_pred(wire, pre)
                                        pre.add_suc(wire, suc_gate)

                            for pred_gate in gate_to_rm.pred[wire]:

                                pred_gate.remove_suc(wire, gate_to_rm)

                        if len(key) != 4:

                            self.commutation_set[key[2]][key[0]].remove(gate_to_rm)

                        else:

                            self.commutation_set[key[1]][key[0]].remove(gate_to_rm)
                            self.commutation_set[key[3]][key[2]].remove(gate_to_rm)

        return

    def detect_commutation(self):
        """
        - Detect commutation relations through matrix multiplication,
        - The DAG will be changed
        - A dictionary set self.commutation_set groups gates for each wire,
        - also provides the index of the commutation of gates on a specific wire.
        """

        # pylint: disable=too-many-nested-blocks
        # pylint: disable=too-many-branches

        for wire in self.wires:
            self.commutation_set[wire] = []

        for wire in self.wires:
            for node in self.qb_tab[wire]:

                if len(self.commutation_set[wire]) == 0:
                    self.commutation_set[wire].append([node])

                elif node not in self.commutation_set[wire][-1]:

                    if commute(node, self.commutation_set[wire][-1][-1]):
                        self.commutation_set[wire][-1].append(node)

                    else:
                        self.commutation_set[wire].append([node])

                self.commutation_set[(node, wire)] = len(self.commutation_set[wire]) - 1

        for wire in self.wires:

            for c_set_ind, c_set in enumerate(self.commutation_set[wire]):

                for node1 in c_set:
                    for node2 in c_set:
                        if node1 != node2:

                            if node2 in node1.suc[wire]:
                                node1.remove_suc(wire, node2)
                                node2.remove_pred(wire, node1)

                    if c_set_ind == len(self.commutation_set[wire]) - 1:
                        continue

                    for next_node in self.commutation_set[wire][c_set_ind + 1]:
                        if next_node not in node1.suc[wire]:
                            node1.add_suc(wire, next_node)
                            next_node.add_pred(wire, node1)
            

        for c_gate in self.op_tab:
            self.suc_set[c_gate] = []
            self.pred_set[c_gate] = []
            for wire in c_gate.wires:
                for c_suc in c_gate.suc[wire]:
                    if c_suc not in self.suc_set[c_gate]:
                        self.suc_set[c_gate].append(c_suc)

                for c_pred in c_gate.pred[wire]:
                    if c_pred not in self.pred_set[c_gate]:
                        self.pred_set[c_gate].append(c_pred)

        self.suc_set[c_gate] = list(set(self.suc_set[c_gate]))
        self.pred_set[c_gate] = list(set(self.pred_set[c_gate]))
    
        self.com_set[c_gate] = self.gate_commutation_set(c_gate)
        return

    def gate_commutation_set(self, gate):
        """
        - Return the gates that commute with a gate
        """

        # pylint: disable=too-many-nested-blocks


        rt_set = []

        if len(gate.wires) == 1:

            g_wire = gate.wires[0]
            gate_key = (gate, gate.wires[0])
            gate_ind = self.commutation_set[gate_key]
            gate_set = self.commutation_set[g_wire][gate_ind]

            for gate_commute in gate_set:

                rt_set.append(gate_commute)

            return rt_set

        for wire in gate.wires:
            
            c_set_ind = self.commutation_set[(gate, wire)]
            c_set = self.commutation_set[wire][c_set_ind]

            for gate_n in c_set:

                if self.op_tab.index(gate_n) != self.op_tab.index(gate):
                    if len(gate_n.wires) == 1:

                        rt_set.append(gate_n)

                    else:
                        comm = True
                        for g_wire in set(gate.wires).intersection(set(gate_n.wires)):

                            gate_ind = (gate, g_wire)
                            gate_n_ind = (gate_n, g_wire)
                            gate_set = self.commutation_set[gate_ind]
                            gate_n_set = self.commutation_set[gate_n_ind]

                            if gate_set != gate_n_set:

                                comm = False

                                break
                        if comm:
                            if gate_n not in rt_set:
                                rt_set.append(gate_n)

        return list(rt_set)

    def _can_merge(self, cgate1, cgate2):

        """
        - This function implements the merge condition
        """

        if cgate1 is cgate2:

            return False

        if set(cgate1.wires) & set(cgate2.wires) == set():

            return False

        rt_val = True

        for wire in set(cgate1.wires) & set(cgate2.wires):

            cgate1_ind = self.commutation_set[(cgate1, wire)]
            cgate2_ind = self.commutation_set[(cgate2, wire)]

            if cgate1_ind != cgate2_ind and (cgate1 not in cgate2.suc[wire]
                                             and
                                             cgate2 not in cgate1.suc[wire]):

                rt_val = False

        return rt_val

    def _can_merge_list(self, cgate_list):

        for c_gate in cgate_list:
            
            c_gate_cannot_merge = True
            for wire in c_gate.wires:
                for c_suc in c_gate.suc[wire]:
                    if c_suc in cgate_list:
                        c_gate_cannot_merge = False

                for c_pred in c_gate.pred[wire]:
                    if c_pred in cgate_list:
                        c_gate_cannot_merge = False

                for n_gate in self.commutation_set[wire][self.commutation_set[(c_gate,wire)]]:

                    if n_gate != c_gate and n_gate in cgate_list:
                        c_gate_cannot_merge = False

            for qg in c_gate.gate_list:
                if qg.name == "H" or qg.name == "Tdag" or qg.name == "T" or qg.name == "S": # Disable for debugging square-root
                    return False
        return not c_gate_cannot_merge
        
    def merge_gate(self, cgate1, cgate2, op_tab=None, qb_tab=None, commutation_set=None):
        """
        - Merge two customized gates.
        - Two gates can be merged only if on every wire, one gate is the pred of the other
          or two commute.
        - Because the circuit is a DAG, it's not possible that on one wire cgate1 is the pred
          and on another wire cgate2 is the pred.
        - Return: the merged gate
        """

        # pylint: disable=too-many-branches
        # pylint: disable=too-many-locals
        # pylint: disable=too-many-statements
        # pylint: disable=too-many-nested-blocks

        if not self._can_merge(cgate1, cgate2):
            return cgate1
        

        c_pred, c_suc = cgate1, cgate2

        wire_set = list(set(cgate1.wires + cgate2.wires))
        wire_interset = list(set(c_pred.wires) & set(c_suc.wires))

        # - Determine which is the pred and suc, if commute, choose cgate1.
        # - The main benefit is not having to write everything twice concerning the two cases.
        for wire in wire_set:

            if cgate1 in cgate2.suc[wire]:

                c_pred, c_suc = cgate2, cgate1
                break

        # Handle the qgate inside between pred and suc
        for wire in wire_interset:

            pred_bfence_ind = c_pred.fence[wire][1]
            suc_ffence_ind = c_suc.fence[wire][0]

            pred_bfence = c_pred.gcirc.op_tab[pred_bfence_ind]
            suc_ffence = c_suc.gcirc.op_tab[suc_ffence_ind]

            if c_suc in c_pred.suc[wire]:

                pred_bfence.suc[wire].remove(c_suc)
                suc_ffence.pred[wire].remove(c_pred)

            pred_bfence.suc[wire].append(suc_ffence)
            suc_ffence.pred[wire].append(pred_bfence)

        rt_gate = custQGate(c_pred.gate_list + c_suc.gate_list)

        # Replace c_pred and c_suc with rt_gate in ob_tab
        # Insert the new cgate, keep the op_tab more or less unchanged

        c_pred_op_tab_ind = self.op_tab.index(c_pred)
        
        self.op_tab.remove(c_pred)
        self.op_tab.remove(c_suc)
        self.op_tab.insert(c_pred_op_tab_ind, rt_gate)
        
        # - Handle the pred-suc relations between the new gate and old preds and sucs.
        # - For each wire, there are 4 (orthogonal) situations need to be handled:
        #   1. c_pred commutes with c_suc / c_pred doesn't commute with c_suc
        #   2. Only one of c_pred and c_suc appears on this wire / Both c_pred and c_suc appear
        #   3. c_pred(c_suc) is in a commutation_set / c_pred(c_suc) is alone in a commutation set
        #   4. c_pred in a commutation set, c_suc alone(or vice versa) / both alone / both in c. set
        # - So there are 2x2x2x3 = 24 cases in total
        # - The following codes try to be as condensed as possible
    
        for wire in wire_set:

            if wire in c_pred.wires:

                c_pred_com_ind = self.commutation_set[(c_pred, wire)]
                c_pred_qb_ind = self.qb_tab[wire].index(c_pred)
                
                self.commutation_set[wire][c_pred_com_ind].remove(c_pred)
                self.qb_tab[wire].remove(c_pred)

                for c_pred_pred in c_pred.pred[wire]:
                    print "c_pred_pred", c_pred_pred.gate_list[0].name
                    c_pred_pred.remove_suc(wire, c_pred)

                    if c_pred in self.suc_set[c_pred_pred]:
                        self.suc_set[c_pred_pred].remove(c_pred)

                for c_pred_suc in c_pred.suc[wire]:

                    c_pred_suc.remove_pred(wire, c_pred)
                    if c_pred in self.pred_set[c_pred_suc]:
                        self.pred_set[c_pred_suc].remove(c_pred)

                rt_gate_pred_set = c_pred.pred[wire]

                rt_gate_suc_set = c_pred.suc[wire]
                rt_gate_commutation_ind = c_pred_com_ind
                rt_gate_qb_ind = c_pred_qb_ind

            if wire in c_suc.wires:

                c_suc_com_ind = self.commutation_set[(c_suc, wire)]
                c_suc_qb_ind = self.qb_tab[wire].index(c_suc)

                self.commutation_set[wire][c_suc_com_ind].remove(c_suc)
                self.qb_tab[wire].remove(c_suc)

                for c_suc_pred in c_suc.pred[wire]:
                    c_suc_pred.remove_suc(wire, c_suc)
                    if c_suc in self.suc_set[c_suc_pred]:
                        self.suc_set[c_suc_pred].remove(c_suc)

                for c_suc_suc in c_suc.suc[wire]:
                    c_suc_suc.remove_pred(wire, c_suc)
                    if c_suc in self.pred_set[c_suc_suc]:
                        self.pred_set[c_suc_suc].remove(c_suc)

                if wire not in c_pred.wires:
                    rt_gate_pred_set = c_suc.pred[wire]
                    rt_gate_commutation_ind = c_suc_com_ind
                rt_gate_suc_set = c_suc.suc[wire]
                rt_gate_qb_ind = c_suc_qb_ind


            if wire in wire_interset and (c_pred_com_ind != c_suc_com_ind):

                if len(self.commutation_set[wire][c_pred_com_ind]) != 0:

                    rt_gate_pred_set = self.commutation_set[wire][c_pred_com_ind]
                    rt_gate_commutation_ind = c_pred_com_ind + 1

                else:
                    del self.commutation_set[wire][c_pred_com_ind]

                    for ind in range(c_pred_com_ind,
                                     len(self.commutation_set[wire])):
                        for cg_comm in self.commutation_set[wire][ind]:
                            self.commutation_set[(cg_comm, wire)] -= 1

                    c_suc_com_ind -= 1
                    rt_gate_commutation_ind = c_pred_com_ind

                if len(self.commutation_set[wire][c_suc_com_ind]) != 0:

                    rt_gate_suc_set = self.commutation_set[wire][c_suc_com_ind]

                else:
                    del self.commutation_set[wire][c_suc_com_ind]
                    for ind in range(c_suc_com_ind,
                                     len(self.commutation_set[wire])):
                        for cg_comm in self.commutation_set[wire][ind]:
                            self.commutation_set[(cg_comm, wire)] -= 1

                self.commutation_set[wire].insert(rt_gate_commutation_ind,
                                                  [])

                for ind in range(rt_gate_commutation_ind + 1,
                                 len(self.commutation_set[wire])):
                    for cg_comm in self.commutation_set[wire][ind]:
                        self.commutation_set[(cg_comm, wire)] += 1
            

            for rt_pred in rt_gate_pred_set:
                rt_pred.add_suc(wire, rt_gate)
                rt_gate.add_pred(wire, rt_pred)
                if rt_gate not in self.suc_set[rt_pred]:
                    self.suc_set[rt_pred].append(rt_gate)
                if rt_pred not in self.pred_set[rt_gate]:
                    self.pred_set[rt_gate].append(rt_pred)

            for rt_suc in rt_gate_suc_set:
                rt_suc.add_pred(wire, rt_gate)
                rt_gate.add_suc(wire, rt_suc)
                if rt_gate not in self.pred_set[rt_suc]:
                    self.pred_set[rt_suc].append(rt_gate)
                if rt_suc not in self.suc_set[rt_gate]:
                    self.suc_set[rt_gate].append(rt_suc)

            self.commutation_set[wire][rt_gate_commutation_ind].append(rt_gate)
            self.commutation_set[(rt_gate, wire)] = rt_gate_commutation_ind
            self.qb_tab[wire].insert(rt_gate_qb_ind, rt_gate)
        
        self.com_set[rt_gate] = self.gate_commutation_set(rt_gate)
        
        for wire in rt_gate.wires:
            rt_com_ind = self.commutation_set[(rt_gate, wire)]
            for cc_gate in self.commutation_set[wire][rt_com_ind]:
                self.com_set[cc_gate] = self.gate_commutation_set(cc_gate)

        rt_gate.et = c_pred.et

        return rt_gate

    def diagonal_merge(self):
        """
        - Find gates to merge into diagonal matrices
        - For commutation-aware scheduling
        """
        # pylint: disable=too-many-branches
        
        print "OP_tab_before_d_merge", len(self.op_tab)
        circ_merge_set = []

        m_mask = defaultdict(lambda: False)

        for gate in self.op_tab:

            if len(gate.wires) == 1:
                continue

            merge_set = []
            commute_merge_set = [[gate]]

            if len(self.gate_commutation_set(gate)) == 0:

                merge_set += commute_merge_set
                suc_set = self._diagonal_merge_suc(gate, commute_merge_set)
                merge_set += suc_set
                merge_set += self._diagonal_merge_pred(gate, commute_merge_set + suc_set)

            else:

                two_qubit_block = True
                for gate_n in self.gate_commutation_set(gate):
                    if set(gate_n.wires).issubset(set(gate.wires)):
                        
                        commute_merge_set_copy = copy.copy(commute_merge_set)
                        for cms in commute_merge_set_copy:
                            commute_merge_set.append(cms+[gate_n])

                    else:

                        two_qubit_block = False

                merge_set += commute_merge_set
                suc_set = self._diagonal_merge_suc(gate, commute_merge_set)
                merge_set += suc_set
                merge_set += self._diagonal_merge_pred(gate, commute_merge_set)

                if two_qubit_block:

                    sort_suc_set = _sort_and_deduplicate(suc_set)
                    sort_merge_set = _sort_and_deduplicate(commute_merge_set)
                    merge_set += self._diagonal_merge_pred(gate,
                                                           sort_suc_set + sort_merge_set)

            max_int = 0
            rt_set = []

            for m_set in merge_set:
                
                m_set = list(set(m_set))
                set_wire = []

                for gate_to_merge in m_set:

                    set_wire = list(set(set_wire) | set(gate_to_merge.wires))

                temp_qc = GCirc(set_wire)

                mset_pass = True
                for gate_to_merge in m_set:
                    for qg_to_merge in gate_to_merge.gate_list:

                        q_gate_to_merge = QGate(qg_to_merge.name,
                                                qg_to_merge.wires,
                                                rot=qg_to_merge.rot)

                        temp_qc.add_gate(q_gate_to_merge)

                    if m_mask[gate_to_merge]:
                        mset_pass = False

                if not mset_pass:
                    continue

                if is_diagonal(temp_qc.unitary()) and len(m_set) >= max_int:

                    if not self._can_merge_list(m_set):
                       continue

                    if len(m_set) > max_int:
                        rt_set = m_set

                        max_int = len(m_set)
                
            if len(rt_set) >= 2:
            
                for rt_set_gate in rt_set:
                    m_mask[rt_set_gate] = True

                circ_merge_set.append(rt_set)
        
        for rt_set in circ_merge_set:

            merged_gate = rt_set[0]
            rt_set.pop(0)

            count = 0 
            
            # while len(rt_set) != 0 and count < len(rt_set) + 1:
            while len(rt_set) != 0:
                if self._can_merge(merged_gate, rt_set[0]):
                    print "Before diag merge"
                    print "merge gate before ind", self.op_tab.index(merged_gate)
                    for qg in merged_gate.gate_list:
                        print qg.name, qg.wires

                    print "his pred and suc"
                    for cc in self.pred_set[merged_gate]:
                        print "pred", self.op_tab.index(cc)
                        for qg in cc.gate_list:
                            print qg.name, qg.wires

                    for cc in self.suc_set[merged_gate]:

                        print "suc", self.op_tab.index(cc)
                        for qg in cc.gate_list:
                            print qg.name, qg.wires

                        
                    print "Before diag merge"
                    print "rt_set before ind", self.op_tab.index(rt_set[0])
                    for qg in rt_set[0].gate_list:
                        print qg.name, qg.wires

                    print "his pred and suc"
                    for cc in self.pred_set[rt_set[0]]:
                        print "pred", self.op_tab.index(cc)
                        for qg in cc.gate_list:
                            print qg.name, qg.wires

                    for cc in self.suc_set[rt_set[0]]:

                        print "suc", self.op_tab.index(cc)
                        for qg in cc.gate_list:
                            print qg.name, qg.wires


                    merged_gate = self.merge_gate(merged_gate, rt_set[0])

                    print "After diag merge"
                    print "merge gate before ind", self.op_tab.index(merged_gate)
                    for qg in merged_gate.gate_list:
                        print qg.name, qg.wires

                    print "his pred and suc"
                    for cc in self.pred_set[merged_gate]:
                        print "pred", self.op_tab.index(cc)
                        for qg in cc.gate_list:
                            print qg.name, qg.wires

                    for cc in self.suc_set[merged_gate]:

                        print "suc", self.op_tab.index(cc)
                        for qg in cc.gate_list:
                            print qg.name, qg.wires

                    rt_set.pop(0)
                    count = 0

                else:
                    count += 1
                    rt_set.append(rt_set.pop(0))

        print "OP_tab_after_d_merge", len(self.op_tab)
        return

    def _old_diagonal_merge_suc(self, gate, current_merge_set):
        
        # Be less flexible and does not consider commutation here.
        
        gate_candidates = []
        wire_max_ind = {}
        for wire in gate.wires:
            gate_wire_ind = self.qb_tab[wire].index(gate)
            wire_max_ind[wire] = gate_wire_ind
            
            start_ind = min(gate_wire_ind + 1, len(self.qb_tab[wire]) - 1)
            for wire_ind in range(start_ind, len(self.qb_tab[wire])):
                gate_to_merge = self.qb_tab[wire][wire_ind]

                if set(gate_to_merge.wires) <= set(gate.wires):
                    gate_candidates.append(gate_to_merge)

                else:
                    wire_max_ind[wire] = wire_ind - 1
                    break
        
        gate_candidates = list(set(gate_candidates))

        for candi_gate in gate_candidates:
            if len(candi_gate.wires) == 2:
                for c_wire in candi_gate.wires:
                    candi_gate_ind = self.qb_tab[wire].index(candi_gate)
                    if candi_gate_ind > wire_max_ind:
                        for comp_wire in gate.wires:
                            comp_candi_gate_ind = self.qb_tab[comp_wire].index(candi_gate)
                            wire_max_ind[wire] = min(comp_candi_gate_ind - 1, wire_max_ind[wire])

        for candi_gate in gate_candidates:
            for wire in candi_gate.wires:
                candi_gate_ind = self.qb_tab[wire].index(candi_gate)
                if candi_gate_ind > wire_max_ind:

                    gate_candidates.remove(candi_gate_ind)

        merge_set = []
        two_qubit_candi = []
        for candi_gate in gate_candidates:
            if len(candi_gate.wires) == 2:
                two_qubit_candi.append(candi_gate)

        two_qubit_candi = sorted(two_qubit_candi,
                                 key = lambda x: self.qb_tab[gate.wires[0]].index(x))
        
        current_whole_set = []
        for tq_ind, tq in enumerate(two_qubit_candi):
            
            tq_set = {}
            for wire in gate.wires:
                
                tq_set[wire] = []
                last_tq_ind = tq_ind - 1
                if last_tq_ind == -1:
                    wire_lower_bound = 0
                else:
                    wire_lower_bound = self.qb_tab[wire].index(two_qubit_candi[last_tq_ind])

                wire_upper_bound = self.qb_tab[wire].index(tq)

                for c_gate in self.qb_tab[wire][wire_lower_bound + 1: wire_upper_bound]:
                    tq_set[wire].append(c_gate)
            
            wire_set_0 = list(map(lambda x:tq_set[gate.wires[0]][:tq_set[gate.wires[0]].index(x)], tq_set[gate.wires[0]]))
            wire_set_1 = list(map(lambda x:tq_set[gate.wires[1]][:tq_set[gate.wires[1]].index(x)], tq_set[gate.wires[1]]))

            for w_set_0 in wire_set_0:
                for w_set_1 in wire_set_1:
                    merge_set.append(current_whole_set + w_set_0 + w_set_1)

            current_whole_set += tq_set[gate.wires[0]] + tq_set[gate.wires[1]] + [tq]
        
        last_set = {}
        for wire in gate.wires:
            last_set[wire] = []
            for candi_gate in gate_candidates:
                if wire in candi_gate.wires:
                    if len(two_qubit_candi) > 0:
                        if self.qb_tab[wire].index(candi_gate) > self.qb_tab[wire].index(two_qubit_candi[-1]):
                            last_set[wire].append(candi_gate)
                    else:
                        last_set[wire].append(candi_gate)


        last_set_0 = list(map(lambda x:last_set[gate.wires[0]][:last_set[gate.wires[0]].index(x)], last_set[gate.wires[0]]))
        last_set_1 = list(map(lambda x:last_set[gate.wires[1]][:last_set[gate.wires[1]].index(x)], last_set[gate.wires[1]]))

        for l_set_0 in last_set_0:
            for l_set_1 in last_set_1:
                merge_set.append(current_whole_set + l_set_0 + l_set_1)
        
        rt_set = merge_set
        merge_set_copy = copy.copy(merge_set)

        for c_set in current_merge_set:
            for m_set in merge_set_copy:

                rt_set.append(m_set + c_set)

        return rt_set

    def _diagonal_merge_suc(self, gate, current_merge_set):
        """
        - Right now doesn't check the combination inside a commutation set,
        which might be helpful in special cases but increase running time.

        - TODO: Having weird bug when printing, temporarily disabling Hadamard
        """

        # pylint: disable=too-many-branches

        rt_set = []
        wire_set = {}
        two_q = []
        for wire in gate.wires:

            wire_set[wire] = []
            gate_set_id = self.commutation_set[(gate, wire)]

            n_max = min(len(self.commutation_set[wire]), gate_set_id + 6) # Only check up to 5 sets

            for set_ind in range(gate_set_id + 1, n_max):

                wire_set[wire].append([])

                two_qubit_block = True

                for gate_n in self.commutation_set[wire][set_ind]:
                    if set(gate_n.wires).issubset(set(gate.wires)):
                        wire_set[wire][-1].append(gate_n)

                        if len(gate_n.wires) == 2 and gate_n not in two_q:

                            two_q.append((gate_n,
                                          self.commutation_set[(gate_n, gate.wires[0])],
                                          self.commutation_set[(gate_n, gate.wires[1])]))
                    else:
                        two_qubit_block = False

                    for qg in gate_n.gate_list:
                        if qg.name == "H":
                            two_qubit_block = False

                if not two_qubit_block:
                    break

        mask = [True] * len(two_q)

        for num, (_, ind_f, ind_s) in enumerate(two_q):

            gate_set_id_0 = self.commutation_set[(gate, gate.wires[0])]
            gate_set_id_1 = self.commutation_set[(gate, gate.wires[1])]

            if ind_f - gate_set_id_0 > len(wire_set[gate.wires[0]]):

                wire_set[gate.wires[1]] = wire_set[gate.wires[1]][:(ind_f - gate_set_id_0)]
                mask[num] = False

            if ind_s - gate_set_id_1 > len(wire_set[gate.wires[1]]):

                wire_set[gate.wires[0]] = wire_set[gate.wires[0]][:(ind_s - gate_set_id_1)]
                mask[num] = False

            two_q = list(compress(two_q, mask))

        for set_ind in range(len(wire_set[gate.wires[0]])+1):

            rt_ele_n = list(reduce(lambda x, y: x + y,
                                   ([[]] + wire_set[gate.wires[0]])[:(set_ind + 1)]))

            for set_ind_s in range(len(wire_set[gate.wires[1]])+1):

                rt_ele = list(set(rt_ele_n +
                                  list(reduce((lambda x, y: x + y),
                                              ([[]] + wire_set[gate.wires[1]])[:(set_ind_s + 1)]))))

                two_qubit_gate_in = False

                for _, two_ind_f, two_ind_s in two_q:
                    if ((two_ind_f >= set_ind and two_ind_s < set_ind_s) or
                            (two_ind_f < set_ind and two_ind_s >= set_ind_s)):

                        two_qubit_gate_in = True

                if two_qubit_gate_in:

                    continue

                for c_set in current_merge_set:

                    rt_set.append(rt_ele + c_set)

        return rt_set

    def _diagonal_merge_pred(self, gate, current_merge_set):

        rt_set = []
        both_pred = [[]]

        # pylint: disable=no-self-use
        # Check 1 1-qubit pred at each wire
        
        for wire in gate.wires:
            if len(gate.pred[wire]) == 1 and len(gate.pred[wire][0].wires) == 1:
                    both_pred.append([gate.pred[wire][0]])

        if len(both_pred) == 2:
            both_pred.append(both_pred[0] + both_pred[1])
        
        for c_set in current_merge_set:
            for p_set in both_pred:
                rt_set.append(c_set + p_set)

        return rt_set

    def cls(self, filename=None):
        """
        - Commutativity aware scheduling
        - Only for 1 qubit and 2 qubit gates
        """
        
        # Initializations
        scheduled = []
        current_com_set = {}
        current_timepoint = [0.0] * len(self.wires)
        com_set = {}
        
        # replaced copy.copy for debugging
        for key in self.commutation_set.keys():
            if isinstance(self.commutation_set[key], int):
                com_set[key] = self.commutation_set[key]
            else: 
                com_set[key] = []
                for term in self.commutation_set[key]:
                    com_set[key].append([])
                    for in_term in term:
                        com_set[key][-1].append(in_term)

        com_set_empty = False
        self.execution_time = {}
        
        for wire in self.wires:
            # Have to use copy here!
            current_com_set[wire] = copy.copy(com_set[wire].pop(0))

        current_gate_set = []
        for wire in self.wires:
            for c_gate in current_com_set[wire]:

                if c_gate in current_gate_set:
                    continue

                elif len(c_gate.wires) == 1:
                    current_gate_set.append(c_gate)

                else: 
                    c_gate_is_current = True
                    for c_wire in c_gate.wires:
                        if c_gate not in current_com_set[c_wire]:
                            c_gate_is_current = False

                    if c_gate_is_current:
                        current_gate_set.append(c_gate)
        counter = 0

        while not com_set_empty or not len(current_gate_set) == 0:

            counter +=1 
            next_gate_ind = 0
            next_gate_set = []    

            while len(next_gate_set) == 0:
                
                next_et = sorted(current_timepoint)[next_gate_ind]
                next_wire_set = [wire for wire_id, wire in enumerate(self.wires) 
                                 if current_timepoint[wire_id] <= next_et]
                
                for c_gate in current_gate_set:
                    if c_gate in next_gate_set:
                        continue 

                    c_gate_next = True
                    for wire in c_gate.wires:
                        if wire not in next_wire_set:
                            c_gate_next = False

                    if c_gate_next:
                        next_gate_set.append(c_gate)

                next_gate_ind += 1
        
            schedule_set = max_matching(next_gate_set)
            
            for c_gate in schedule_set:

                c_gate.et = next_et
                scheduled.append((c_gate, c_gate.et))
                self.execution_time[c_gate] = next_et
                current_gate_set.remove(c_gate)
                next_gate_set.remove(c_gate)
            
                for wire in c_gate.wires:

                    cw_id = self.wires.index(wire)
                    current_com_set[wire].remove(c_gate)
                    current_timepoint[cw_id] = c_gate.duration() + c_gate.et
            
            for wire in self.wires:

                if (len(current_com_set[wire]) == 0
                    and len(com_set[wire]) != 0):
            
                    current_com_set[wire] = copy.copy(com_set[wire].pop(0))

                for c_gate in current_com_set[wire]:

                    if c_gate in current_gate_set:
                        continue

                    if len(c_gate.wires) == 1:
                        current_gate_set.append(c_gate)

                    else: 
                        c_gate_is_current = True
                        for c_wire in c_gate.wires:
                            if c_gate not in current_com_set[c_wire]:
                                c_gate_is_current = False


                        if c_gate_is_current:
                            current_gate_set.append(c_gate)
                
            com_set_empty = True
            for wire in self.wires:
                if len(com_set[wire]) > 0:
                    com_set_empty = False

        # Update the critical path
        sorted_gate_list = sorted(self.execution_time.items(), key = lambda t: t[0].duration() + t[1])

        total_length = sorted_gate_list[-1][1]
        
        last_gates = [u for u, v in sorted_gate_list if u.duration() + v == total_length]

        self.critical_paths += last_gates
        
        for c_gate in last_gates:
            self._find_critical_paths(c_gate)
        
        for c_gate in self.op_tab:
            self._find_critical_paths(c_gate, target_gate=c_gate)

        for c_gate in self.op_tab:
            for wire in c_gate.wires:
                for c_suc in c_gate.suc[wire]:
                    if c_suc not in self.suc_set[c_gate]:
                        self.suc_set[c_gate].append(c_suc)

                for c_pred in c_gate.pred[wire]:
                    if c_pred not in self.pred_set[c_gate]:
                        self.pred_set[c_gate].append(c_pred)

        self.suc_set[c_gate] = list(set(self.suc_set[c_gate]))
        self.pred_set[c_gate] = list(set(self.pred_set[c_gate]))
    
        self.com_set[c_gate] = self.gate_commutation_set(c_gate)

        if filename != None:
            with open(filename, "w") as f_handler:
                for wire in self.wires:
                    f_handler.write("qubit {}\n".format(wire))
                for c_gate, gate_time in scheduled:
                    for q_gate in c_gate.gate_list:
                        # if q_gate.name in ("Rz", "Rx", "Ry"):
                        #     f_handler.write("{} {} {}\n".format(q_gate.name,
                        #                                         q_gate.rot,
                        #                                         ",".join(q_gate.wires)))
                        # else:
                        #     f_handler.write("{} {}\n".format(q_gate.name,
                        #                                      ",".join(q_gate.wires)))

                        if q_gate.name in ("Rz", "Rx", "Ry"):
                            f_handler.write("{} {} {} {}\n".format(gate_time, q_gate.name,
                                                                q_gate.rot,
                                                                ",".join(q_gate.wires)))
                        else:
                            f_handler.write("{} {} {}\n".format(gate_time, q_gate.name,
                                                             ",".join(q_gate.wires)))

        return scheduled

    def _find_critical_paths(self, c_gate, target_gate=None):
        
        if target_gate == None:
            store_list = self.critical_paths
        else:
            store_list = self.paths[target_gate]

        for wire in c_gate.wires:
            if len(c_gate.pred[wire]) == 0:
                continue

            for c_pred in c_gate.pred[wire]:
                if c_pred.et + c_pred.duration() == c_gate.et:
                    if c_pred not in store_list:

                        store_list.append(c_pred)
                        self._find_critical_paths(c_pred, target_gate=target_gate)
            
            c_gate_com_ind = self.commutation_set[(c_gate, wire)]
            for com_gate in self.commutation_set[wire][c_gate_com_ind]:

                if com_gate.et + com_gate.duration() == c_gate.et:
                    if com_gate not in store_list:

                        store_list.append(com_gate)
                        self._find_critical_paths(com_gate, target_gate=target_gate)
        return

    def block_merge(self, wire_width):
        """
        - Generate the final n-qubit blocks
        """

        c_merge_set = {} 

        in_merge_set = defaultdict(lambda: [])
        unmerged_gates = []


        for c_gate in self.op_tab:
            
            c_merge_set[c_gate] = self._max_merge_set([c_gate], wire_width)

            for m_gate in c_merge_set[c_gate]:
                in_merge_set[m_gate].append(c_gate)
         
        for c_gate in self.op_tab:

            if len(c_merge_set[c_gate]) > 1:
                unmerged_gates.append(c_gate)

        while len(unmerged_gates) != 0:

            print "unmerged_gates", len(unmerged_gates)
            for cc_gate in unmerged_gates:
                print "C gate"

                for qg in cc_gate.gate_list:
                    print qg.name, qg.wires
            
            c_gate = max(unmerged_gates,
                         key = lambda cg: len(c_merge_set[cg]))
            
            while c_gate not in self.op_tab:
                unmerged_gates.remove(c_gate)
                c_gate = max(unmerged_gates,
                             key = lambda cg: len(c_merge_set[cg]))

            c_merge_set[c_gate] = self._max_merge_set([c_gate], wire_width)
            gates_to_update = []
            unmerged_gates.remove(c_gate)
            c_merge_set[c_gate].remove(c_gate)

            merged_gate = c_gate
            for g_to_u in in_merge_set[merged_gate]:
                if g_to_u != c_gate:
                    gates_to_update.append(g_to_u)

            print "I am in c_gate", c_gate.gate_list[0].name, c_gate.gate_list[0].wires, self.op_tab.index(c_gate)
            print "c_merge_set len", len(c_merge_set[c_gate])
            for cc_gate in c_merge_set[c_gate]:
                print "C gate", 
                for qg in cc_gate.gate_list:
                    print qg.name, qg.wires
            
            counter = 0

            while len(c_merge_set[c_gate]) != 0:

                print "OP_tab", len(self.op_tab)
                print "##########"
                counter += 1
                if counter > 1000:
                    break
                print "counter", counter

                print "Before merge"
                print "merge gate before ind", self.op_tab.index(merged_gate)
                for qg in merged_gate.gate_list:
                    print qg.name, qg.wires

                print "his pred and suc"
                for cc in self.pred_set[merged_gate]:
                    print "pred", self.op_tab.index(cc)
                    for qg in cc.gate_list:
                        print qg.name, qg.wires

                for cc in self.suc_set[merged_gate]:

                    print "suc", self.op_tab.index(cc)
                    for qg in cc.gate_list:
                        print qg.name, qg.wires
                
                if c_merge_set[c_gate][0] not in self.op_tab:
                    # c_merge_set[c_gate] = self._max_merge_set([c_gate], wire_width)
                    c_merge_set[c_gate].pop(0)
                    break
                print "first gate in cmerge[cgate] ind", self.op_tab.index(c_merge_set[c_gate][0])

                for qg in c_merge_set[c_gate][0].gate_list:
                    print qg.name, qg.wires

                print "his pred"
                for cc in self.pred_set[c_merge_set[c_gate][0]]:
                    print "pred", self.op_tab.index(cc)
                    for qg in cc.gate_list:
                        print qg.name, qg.wires

                print "his suc"
                for cc in self.suc_set[c_merge_set[c_gate][0]]:
                    print "suc", self.op_tab.index(cc)
                    for qg in cc.gate_list:
                        print qg.name, qg.wires

                if self._can_merge(merged_gate, c_merge_set[c_gate][0]):
                    
                    merged_gate = self.merge_gate(merged_gate, c_merge_set[c_gate][0])
                    print "Merged_gate index", self.op_tab.index(merged_gate)

                    gate_removed = c_merge_set[c_gate].pop(0)
                    print "Gate I think I removed"

                    print "___OPTAB len", len(self.op_tab)
                    
                    if c_gate in in_merge_set[gate_removed]:
                        in_merge_set[gate_removed].remove(c_gate)

                    gates_to_update += in_merge_set[gate_removed]
                    gates_to_update = list(set(gates_to_update))
                    
                    if gate_removed in unmerged_gates:
                        print "Now I am removing", gate_removed.gate_list[0].name, gate_removed.wires, "from unmerged"
                        unmerged_gates.remove(gate_removed)
                    

                    print "now check neighbors and sucs and preds"

                    for ccc in self.op_tab:
                        print "ccc", ccc
                        print "sucs"
                        for wire in ccc.wires:
                            print "wire", wire
                            for suc_ccc in ccc.suc[wire]:
                                print "suc", suc_ccc
                                for qg in suc_ccc.gate_list:
                                    print qg.name, qg.wires
                                print "suc_ind", self.op_tab.index(suc_ccc)
                    
                        print "preds"

                        for pred_ccc in self.pred_set[ccc]:
                            print "pred"
                            for qg in pred_ccc.gate_list:
                                print qg.name, qg.wires
                            print "suc_ind", self.op_tab.index(pred_ccc)
                            
                else:
                    print "Nothing happened"
                    c_merge_set[c_gate].append(c_merge_set[c_gate].pop(0))

            c_merge_set[merged_gate] = self._max_merge_set([merged_gate], wire_width)

            if len(c_merge_set[merged_gate]) != 1:
                unmerged_gates.append(merged_gate)

            for g_to_u in gates_to_update:
                c_merge_set[g_to_u] = self._max_merge_set([g_to_u], wire_width)

                if(len(c_merge_set[g_to_u]) == 0):
                    unmerged_gates.remove(g_to_u)

    def _max_merge_set(self, c_gate_list, wire_width):
        """
        - Return the max merge set of a gate, given a wire width limit
        """
        
        c_merge_set = self._find_merge_set(c_gate_list, wire_width)
        if len(c_merge_set) == len(c_gate_list):
            return c_merge_set
        
        wire_set = []
        for c_gate in c_merge_set:
            wire_set += c_gate.wires
        
        wire_set = list(set(wire_set))
        if len(wire_set) > wire_width:
            return c_gate_list

        else:
            return self._find_merge_set(c_merge_set, wire_width)

    def _find_merge_set(self, c_gate_list, wire_width):
        """
        - Find the max merge set of a gate list, 
        - withouting having to merging them
        """
        rt_set = []

        for c_gate in c_gate_list:
            
            rt_set = list(set(self._find_gate_merge_set(c_gate, wire_width))
                       | set(c_gate_list)
                       | set(rt_set))

        return rt_set

    def _find_gate_merge_set(self, c_gate, wire_width):
        """
        - Find the set of gates can be merged with c_gate,
          and not compromising critical path
        """

        c_merge_set = []
        for gate_to_merge in list(set(self.pred_set[c_gate]) | set(self.com_set[c_gate])):

            gate_to_merge_pass =True
            for check_gate in list(set(self.pred_set[c_gate]) |  set(self.com_set[c_gate])):
                if (check_gate.et + check_gate.duration()
                    > gate_to_merge.et
                    and 
                    check_gate != gate_to_merge):

                    gate_to_merge_pass = False
                    break
            
            if not gate_to_merge_pass:
                continue

            for check_gate in list(set(self.suc_set[gate_to_merge]) | set(self.com_set[gate_to_merge])):
                if (check_gate.et < c_gate.duration() + c_gate.et
                    and 
                    check_gate != c_gate):

                    gate_to_merge_pass = False
                    break

            if not gate_to_merge_pass:
                continue

            if gate_to_merge not in c_merge_set:
                c_merge_set.append(gate_to_merge)

        for gate_to_merge in list(set(self.suc_set[c_gate]) | set(self.com_set[c_gate])):
            gate_to_merge_pass =True
            for check_gate in self.pred_set[gate_to_merge] + self.com_set[gate_to_merge]:
                if (check_gate.et + check_gate.duration()
                    > c_gate.et
                    and 
                    check_gate != c_gate):

                    gate_to_merge_pass = False
                    break
            
            if not gate_to_merge_pass:
                continue

            for check_gate in list(set(self.suc_set[c_gate]) |set( self.com_set[c_gate])):
                if (check_gate.et > gate_to_merge.duration() + gate_to_merge.et
                    and 
                    check_gate != c_gate):

                    gate_to_merge_pass = False
                    break

            if not gate_to_merge_pass:
                continue
            
            if gate_to_merge not in c_merge_set:
                c_merge_set.append(gate_to_merge)

        return c_merge_set

    ########################################
    #    - This approach seems not promising
    #    - Save it here in case needed in the future

    #    - Bound the execution time of a gate by analyzing the commutation sets
    #      on wires. This way we can determine if two gates can be merged without
    #      scheduling and losing commutation information.
    #    - Helps speedup the merge
    #    - Commutation_bound is the bound for next commutation_set.
    #    - Execution_bound is the bound for the current gate.
    #    - The bounds are sort of approximation.
        
    #   def _gate_execution_bound(self):
    #     execution_bound = defaultdict(lambda: -1)
    #     commutation_bound = defaultdict(lambda: [])

    #     com_set_empty = False
    #     current_com_set = {}
    #     com_set = copy.deepcopy(self.commutation_set)

    #      for wire in self.wires:

    #         current_com_set[wire] = com_set[wire].pop(0)
    #         commutation_bound[wire].append(0.0)
        #     
        #     for c_gate in current_com_set[wire]:

        #          c_gate_determined = True
        #         for c_wire in c_gate.wires:
        #             if c_gate not in curren_com_set[c_wire]:
        #                 c_gate_determined = False

        #          execution_bound[c_gate] = 0.0
        #         commutation_bound[wire][-1] += (execution_bound[c_gate]
        #                                            + c_gate.duration())

        #          if not c_gate_determined:
        #             continue

        #          current_com_set[wire].remove(c_gate)

        #  while not com_set_empty:
        #     
        #     for wire in self.wires:
        #         if (len(current_com_set[wire]) == 0
        #             and len(com_set[wire]) != 0):

        #              current_com_set[wire] = com_set[wire].pop(0)
        #             commutation_bound[wire][-1]

        #          for c_gate in current_com_set[wire]:
        #             
        #             c_gate_determined = True
        #             for c_wire in c_gate.wires:
        #                 if c_gate not in current_com_set[c_wire]:
        #                     c_gate_determined = False
        #             
        #             if len(commutation_bound[wire]) >= 2:

        #                  execution_bound[c_gate] = commutation_bound[wire][-2]

        #              if commutation_bound[c_wire][-1] < (execution_bound[c_gate]
        #                                                 + c_gate.duration()):
        #                 commutation_bound[wire][-1] = (execution_bound[c_gate]
        #                                                + c_gate.duration())
        #             if not c_gate_determined:
        #                 continue
        #             

        #              current_com_set[wire].remove(c_gate)

    def total_length(self):

        return sorted(self.execution_time.items(), key = lambda t: t[0].duration() + t[1])[-1][1]

###########################
# Global helper functions
###########################

def _simple_product(cgate1, cgate2):
    """
    - Helper functions
    - The product of the two input unitaries(of 1 or 2 qubits)
    """

    wires = sorted(set(cgate1.wires + cgate2.wires))
    temp_qc = GCirc(wires)

    for gate in cgate1.gate_list:
        temp_qc.add_gate(QGate(gate.name, gate.wires, rot=gate.rot))

    for gate in cgate2.gate_list:
        temp_qc.add_gate(QGate(gate.name, gate.wires, rot=gate.rot))
    
    return temp_qc.unitary()

def _matrix_commute(cgate1, cgate2):
    """
    - Commutation detection based on matrix multiplication
    """

    if set(cgate1.wires) & set(cgate2.wires) == set():

        return True

    unitary_fs = _simple_product(cgate1, cgate2)
    unitary_sf = _simple_product(cgate2, cgate1)

    if unitary_fs is None:
        return False

    # np.array_equal fails here
    return np.allclose(unitary_fs, unitary_sf, atol=1E-6, rtol=1E-4)

def commute(cgate1, cgate2):
    """
    - Determine if two gates commute
    """

    return _matrix_commute(cgate1, cgate2)

def is_diagonal(array):
    """
    - Determine if a numpy is diagonal
    """

    # return not bool(np.count_nonzero(array - np.diag(np.diagonal(array))))
    return np.allclose(array, np.diag(np.diagonal(array)))

def _uniq(lst):
    """
    - Helper function for sort and deduplicate
    """

    last = object()
    for item in lst:
        if item == last:
            continue
        yield item
        last = item

def _sort_and_deduplicate(array):
    """
    - Needed when the iterms are not iterables and can't use python set
    """
    return list(_uniq(sorted(array, reverse=True)))

def max_matching(gate_set):
    """
    - Find the maximal matching of the current interaction graph set        
    """

    G = nx.MultiGraph()
    
    wire_set = []
    for gate in gate_set:
        wire_set += gate.wires
    
    wire_set = list(set(wire_set))

    for wire_id, wire in enumerate(wire_set):
        
        G.add_node(wire_id)
        G.node[wire_id]['name'] = wire
    
    for gate_id, gate in enumerate(gate_set):

        if len(gate.wires) == 1:
            
            G.add_edge(gate.wires[0], gate.wires[0], name = gate_id)
        
        if len(gate.wires) == 2:

            G.add_edge(gate.wires[0], gate.wires[1], name = gate_id)

    p_matching = _maximal_matching(G)
    rt_set = list(map(lambda x: gate_set[x[2]], p_matching))

    return rt_set

def _maximal_matching(G):
    """
    - Helper function to calculate the maximal matching
    """

    matching = set()
    nodes = set()
    for u, v, d in G.edges(data = True):

	# If the edge isn't covered, add it to the matching
	# then remove neighborhood of u and v from consideration.
	if u not in nodes and v not in nodes:

            l = d['name']
	    matching.add((u, v, l))
	    nodes.add(u)
	    nodes.add(v)

    return matching
