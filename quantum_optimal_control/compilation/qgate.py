# Quantum gate class for intrinsic gates
from simulator import _GateMasterDef

class QGate():

    def __init__(self, op, args, rot = 0.):

        self.name = op # Gate name
        self.wires = args
        self.id = -2
        self.rot = rot
        self.layer = 0 
        self.length = 0
        self.time = 0.0
        self.matrix = _GateMasterDef(name = op, para = rot)
        self.suc = {}
        self.pred = {}

        for wire in self.wires:
            self.suc[wire] = []
            self.pred[wire] = []

    def unitary(self):
        return _self.matrix

