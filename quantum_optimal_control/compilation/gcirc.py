"""
 - Circuit class assigned to each customized gate, contains intrinsic quantum gate
 - Separated from QCirc to avoid loop import in custQGate
 - Assume to be small enough to be simulated
 - Layers are for simulation, not for scheduling
 - Omit commutation because customized gates are either
   for diagonal matrix of 2 qubits or for optimal control,
   in both case, the commutation doesn't matter too much.
 - Require config/gate.config file for gate initialization
"""

from quantum_optimal_control.compilation.simulator import Simulator
from quantum_optimal_control.compilation.errors import do_error
from collections import defaultdict

class GCirc(object):
    """
    - Circuit class for each customized gate
    """

    def __init__(self, wires, filename = "../config/gate.config"):

        self.wires = wires
        self.qb_tab = {}  # keep track of the gates on each wire
        self.qb2idx = {}  # helpful when calculating the unitary

        for wire_num, wire in enumerate(self.wires):
            self.qb_tab[wire] = []
            self.qb2idx[wire] = wire_num + 1

        self.op_tab = []
        self.layers = [] # helpful when calculating the unitary, make the simulator faster.
        self.time = defaultdict(lambda: [])
        
        self.gate_master_duration = {}

        with open(filename, "r") as f:
            for line in f:

                dat = line.split(' ')
                if dat[0] == '#':
                    continue

                self.gate_master_duration[dat[0]] = float(dat[1])
        
    def add_gate(self, qgate):
        """
        - For reading from the parser
        """

        self.op_tab.append(qgate) # Put the gate into the gate list
        qgate.id = len(self.op_tab) - 1 # Give the gate a unique ID number.
        qgate.length = self.gate_master_duration[qgate.name]

        for wire in qgate.wires:
            if wire not in self.qb_tab:
                error_msg = (wire, qgate.id, qgate.name + " " + wire)
                do_error('[QCircuit] No qubit {0} for gate # {1}: {2}'\
                         .format(error_msg[0], error_msg[1], error_msg[2]))

            # Determine the layer for the gate, mainly for generating unitary matrix faster
            # Also determine the timing of the gate
            if len(self.qb_tab[wire]) == 0:
                layer = 1

                timing = 0.0

            else:

                last_gate = self.op_tab[self.qb_tab[wire][-1]]

                layer = last_gate.layer + 1
                last_gate.suc[wire].append(qgate)
                qgate.pred[wire].append(last_gate)
                
                last_gate_length = self.gate_master_duration[last_gate.name]
                timing = last_gate.time + last_gate_length

            self.qb_tab[wire].append(qgate.id)

            if layer > qgate.layer:
                qgate.layer = layer

            if timing > qgate.time:
                qgate.time = timing

        if qgate.layer > len(self.layers):

            self.layers.append([])

        self.layers[qgate.layer -1].append(qgate.id)
        
        self.time[qgate.time].append(qgate)

    def fence(self, wire):
        """
        - Return the fences on each wire.
        """

        return (self.qb_tab[wire][0], self.qb_tab[wire][-1])

    def output_layer(self):	
        """
        - Print the layers of the circuit
        """

        for n_layer, layer in enumerate(self.layers):	# loop over timesteps
            print "%%  Layer %02d:" % (n_layer + 1)
            for gate in layer:		# loop over events in this timestep
                q_op = self.op_tab[gate]
                print "%%    Gate %02d %s(%s)" % (q_op.id,
                                                  q_op.name, ','.join(q_op.wires))

    def duration(self):
        """
        - Return the duration of the circuit
        """
         
        last_et, last_gates = sorted(self.time.items(), key = lambda t: t[0])[-1]
        
        duration = 0.0
        
        max_gate = max(last_gates, key = lambda x: x.length)
        
        return max_gate.length + last_et


    def unitary(self):
        """
        - Return the unitary of the circuit
        """

        return Simulator.calculate_unitary(self)
