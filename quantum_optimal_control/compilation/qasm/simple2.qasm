qubit reg0
qubit reg1
qubit reg2
qubit reg3
qubit reg4
H reg0
H reg1
H reg2
H reg3
H reg4
Rz -0.300000 reg0
Rz 0.300000 reg1
Rz -1.200000 reg1
CNOT reg0,reg1
Rz 0.600000 reg1
CNOT reg0,reg1
Rz -0.360000 reg2
Rz 0.360000 reg3
Rz -1.440000 reg3
CNOT reg2,reg3
Rz 0.720000 reg3
CNOT reg2,reg3
Rz -0.120000 reg4
Rz 0.260000 reg1
Rz -0.260000 reg2
Rz 1.040000 reg2
CNOT reg1,reg2
Rz -0.520000 reg2
CNOT reg1,reg2
Rz -0.260000 reg3
Rz 0.260000 reg4
Rz -1.040000 reg4
CNOT reg3,reg4
Rz 0.520000 reg4
CNOT reg3,reg4
