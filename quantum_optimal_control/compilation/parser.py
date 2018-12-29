"""
 - Parse qasm files
 - No white spaces between qubit arguments but a comma
"""

import re
from quantum_optimal_control.compilation.qgate import QGate
from quantum_optimal_control.compilation.errors import do_error

class QasmParser(object):
    """
    - Parser class
    """

    # pylint: disable=anomalous-backslash-in-string
    # pylint: disable=too-few-public-methods

    def __init__(self, filename):
        self.qubits = []           #list of qubit names
        self.gates = []            #list of gates
        self.comments = ''         #string to store the comments

        error_line = 0             #for printing out error line

        with open(filename) as f_handler:

            for line in f_handler:

                error_line += 1

                if line[0] == '#':
                    self.comments += line
                    continue

                # qubit declaration, ignore cbits here
                match = re.compile('\s*qubit\s+(\S+)').search(line)
                if match:
                    self.qubits.append(match.group(1))	# add name
                    continue

                # gate
                match = re.compile('^([^Rq]+)\s+(\S+)').search(line)
                if match:
                    op_name = match.group(1)
                    args = match.group(2)
                    gate_check(op_name, args, error_line)
                    self.gates.append(QGate(op_name, args.split(',')))

                # rotations
                match = re.compile('^([Rxyz]+)\s+(\S+)\s+(\S+)').search(line)
                if match:
                    op_name = match.group(1)
                    rot = match.group(2)
                    args = match.group(3)
                    gate_check(op_name, args, error_line)
                    self.gates.append(QGate(op_name, args.split(','), rot=float(rot)))


def gate_check(op_name, args, error_line):
    """
    - Sanity check for gates
    """

    if op_name not in GATE_DEF:
        msg = (error_line, op_name, args)
        do_error("Line %d unknown gate op %s on %s" % msg)

    n_wires, _ = GATE_DEF[op_name]

    if len(args.split(',')) != n_wires:
        msg = (error_line, op_name + " " + args)
        do_error("Line %d wrong number of qubits in %s" % msg)

    # Non-cloning check
    wires = args.split(',')
    if [wires.count(qb) for qb in wires].count(1) < len(wires):
        msg = (error_line, op_name + " " + args)
        do_error("Line %d duplicate qubits in %s" % msg)

    return True

GATE_DEF = {'CNOT'     : (2, 1),
            'Rx'       : (1, 0),
            'Ry'       : (1, 0),
            'Rz'       : (1, 0),
            'CZ'       : (2, 1),
            'measure'  : (1, 0),
            'h'        : (1, 0),
            'H'        : (1, 0),
            'X'        : (1, 0),
            'Y'        : (1, 0),
            'Z'        : (1, 0),
            'S'        : (1, 0),
            'T'        : (1, 0),
            'Tdag'     : (1, 0),
            'swap'     : (2, 0),
           }
