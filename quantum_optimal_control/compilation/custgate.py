# Aggregated gate class
from gcirc import *
from collections import defaultdict
class custQGate:

    def __init__(self, gate_list):

        self.wires = []
        self.gate_list = gate_list
        self.suc = defaultdict(lambda: [])
        self.pred = defaultdict(lambda: [])
        self.et = 0.0
        self.fence = {}
        
        for gate in self.gate_list:
            for g_wire in gate.wires:
                if not g_wire in self.wires:

                    self.wires.append(g_wire)

        self.gcirc = GCirc(self.wires)

        for qgate in self.gate_list:

            self.gcirc.add_gate(qgate)

        for wire in self.wires:
            self.fence[wire] = self.gcirc.fence(wire)

    def __eq__(self, other):

        if not isinstance(other, custQGate):
            return False

        if self.wires != other.wires:

            return False
        
        if self.gate_list != other.gate_list:

            return False

        return True

    def __hash__(self):

        return id(self)

    def duration(self):

        return self.gcirc.duration()

    def add_suc(self, wire, cgate):

        self.suc[wire].append(cgate)

        for qgate in self.gate_list:
            if wire in qgate.wires:
                if (len(qgate.suc[wire]) == 0 or
                    isinstance(qgate.suc[wire][0], custQGate)):

                    qgate.suc[wire].append(cgate)
    
    def add_pred(self, wire, cgate):
        self.pred[wire].append(cgate)
        for qgate in self.gate_list:
            if wire in qgate.wires:
                if len(qgate.pred[wire]) == 0 or isinstance(qgate.pred[wire][0], custQGate):
                    qgate.pred[wire].append(cgate)

    def remove_suc(self, wire, cgate):
        self.suc[wire].remove(cgate)
        for qgate in self.gate_list:
            if wire in qgate.wires and cgate in qgate.suc[wire]:
                    qgate.suc[wire].remove(cgate)

    def remove_pred(self, wire, cgate):
        self.pred[wire].remove(cgate)
        for qgate in self.gate_list:
            if wire in qgate.wires and cgate in qgate.pred[wire]:
                    qgate.pred[wire].remove(cgate)

    def unitary(self):
        return self.gcirc.unitary()
