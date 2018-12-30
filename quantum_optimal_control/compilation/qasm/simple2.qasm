qubit reg0
qubit reg1
qubit reg2
qubit reg3

H reg0
CNOT reg0,reg1
H reg0
CNOT reg0,reg1
CNOT reg1,reg0
CNOT reg2,reg1
Rz -1.2 reg1
CNOT reg3,reg2
CNOT reg0,reg3
CNOT reg1,reg0
