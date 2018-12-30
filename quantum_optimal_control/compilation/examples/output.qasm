qubit reg0
qubit reg1
qubit reg2
qubit reg3
0.0 H reg3
0.0 H reg1
0.0 H reg2
0.0 H reg0
0.4 Rz 0.3 reg1
0.4 Rz 0.36 reg3
0.4 Rz -0.36 reg2
0.4 Rz -0.3 reg0
0.7 Rz -1.44 reg3
0.7 Rz -1.2 reg1
0.7 Rz -0.26 reg2
1.0 Rz 1.04 reg2
1.0 CNOT reg0,reg1
1.3 CNOT reg2,reg3
2.2 Rz 0.6 reg1
2.5 Rz 0.72 reg3
2.5 CNOT reg0,reg1
2.8 CNOT reg2,reg3
3.7 Rz 0.26 reg1
3.7 H reg0
4.0 CNOT reg1,reg2
4.0 H reg3
4.1 Rz -1.92 reg0
4.4 H reg0
4.8 Rz -0.288 reg0
5.2 Rz -0.52 reg2
5.5 CNOT reg1,reg2
6.7 H reg1
6.7 H reg2
7.1 Rz -1.92 reg1
7.1 Rz -1.92 reg2
7.4 H reg1
7.4 H reg2
7.8 Rz 0.864 reg1
7.8 Rz 1.152 reg2
