qubit reg0
qubit reg1
qubit reg2
qubit reg3
H reg0
H reg1
H reg2
H reg3
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
Rz 0.260000 reg1
Rz -0.260000 reg2
Rz 1.040000 reg2
CNOT reg1,reg2
Rz -0.520000 reg2
CNOT reg1,reg2
H reg0
Rz -1.920000 reg0
H reg0
Rz -0.288000 reg0
H reg1
Rz -1.920000 reg1
H reg1
Rz 0.864000 reg1
H reg2
Rz -1.920000 reg2
H reg2
Rz 1.152000 reg2
H reg3
