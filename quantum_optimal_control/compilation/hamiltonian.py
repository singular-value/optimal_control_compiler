import numpy as np
import scipy as sp
import scipy.linalg
import re
import sys
import numpy.random

def Q_x(qubit_state_num): return np.diag(np.sqrt(np.arange(1,qubit_state_num)),1)+np.diag(np.sqrt(np.arange(1,qubit_state_num)),-1) 
def Q_y(qubit_state_num): return (0+1j) *(np.diag(np.sqrt(np.arange(1,qubit_state_num)),1)-np.diag(np.sqrt(np.arange(1,qubit_state_num)),-1))
def Q_z(qubit_state_num): return np.array([[1.0,0.0],[0.0,-1.0]], dtype = np.complex)
def I_q(qubit_state_num): return np.identity(qubit_state_num)

def rz(theta):
    return np.array([[np.exp(-1j * theta / 2), 0],[0, np.exp(1j * theta / 2)]])
def rx (theta):
    return np.array([[np.cos(theta / 2), -1j * np.sin(theta / 2)], [-1j * np.sin(theta / 2), np.cos(theta / 2)]])
def ry(theta):
    return np.array([[np.cos(theta / 2), - np.sin(theta / 2)], [np.sin(theta / 2), np.cos(theta / 2)]], dtype = np.complex)

def kron_list(args):
    ret = args[0]
    for item in args[1:]:
        ret = np.kron(ret, item)
    return ret

def ion_trap_Hamiltonian(qubit_num, qubit_state_num, amp1, amp2):
    
    Q_x = np.array([[0.0, 1.0],[1.0,0.0]], dtype = np.complex)
    Q_y = np.array([[0.0, 1.0j],[-1.0j,0.0]], dtype = np.complex)
    I_q = np.identity(2)

    ret_list, str_list, amp_list = [], [],[]
    str_dict = {'i':I_q, 'x':Q_x, 'y':Q_y}
    for i in range(qubit_num):
        str_x, str_y  = ['i']*qubit_num, ['i']*qubit_num
        str_x[i],str_y[i] = 'x','y'
        ind_x, ind_y = [str_dict[st] for st in str_x], [str_dict[st] for st in str_y]

        ret_list.append(kron_list(ind_x))
        ret_list.append(kron_list(ind_y))
        str_list.append(''.join(str_x))
        str_list.append(''.join(str_y))
        
        amp_list.append(amp1)
        amp_list.append(amp1)

        for j in range(i, qubit_num):
            str_xx, str_yy  = ['i']*qubit_num, ['i']*qubit_num
            str_xx[i],str_xx[j],str_yy[i], str_yy[j] = 'x','x','y','y'
            ind_xx, ind_yy = [str_dict[st] for st in str_xx], [str_dict[st] for st in str_yy]

            ret_list.append(kron_list(ind_xx))
            ret_list.append(kron_list(ind_yy))
            str_list.append(''.join(str_xx))
            str_list.append(''.join(str_yy))

            amp_list.append(amp2)
            amp_list.append(amp2)

    return ret_list, str_list, amp_list

def xy_hamiltonian_2D(qubit_num, qubit_state_num, amp1, amp2, qubit_mapping):
    # mapping is in the form of a list of tuples
    Q_x = np.array([[0.0, 1.0],[1.0,0.0]], dtype = np.complex)
    Q_y = np.array([[0.0, 1.0j],[-1.0j,0.0]], dtype = np.complex)
    Q_y = np.array([[1.0, 0.0],[0.0,-1.0]], dtype = np.complex)
    I_q = np.identity(2)
    ret_list, str_list, amp_list = [], [],[]
    # note the use of 'x' for 'xy' coupling
    str_dict = {'i':I_q, 'x':Q_x, 'y':Q_y, 'z':Q_z}
    
    # keep track of connections already made
    # improvement comes through further using hash tables
    visited_connections = []
    visited_nodes       = {qubit_mapping[i]: 0 for i in range(len(qubit_mapping))}

    for i in range(qubit_num):
        str_x, str_y, str_z  = ['i']*qubit_num, ['i']*qubit_num, ['i']*qubit_num
        str_x[i],str_y[i], str_z[i] = 'x','y', 'z'
        ind_x, ind_y, ind_z = [str_dict[st] for st in str_x], [str_dict[st] for st in str_y],[str_dict[st] for st in str_z]


        ret_list.append(kron_list(ind_x))
        ret_list.append(kron_list(ind_y))
    #    ret_list.append(kron_list(ind_z))
        str_list.append(''.join(str_x))
        str_list.append(''.join(str_y))
    #    str_list.append(''.join(str_z))
        
        amp_list.append(amp1)
        amp_list.append(amp1)
    #    amp_list.append(amp1)

    for index_0 in range(len(qubit_mapping)):

        coord = qubit_mapping[index_0]

        if not visited_nodes[coord] == 1:

            visited_nodes[coord] = 1
            adj_0 = (coord[0] + 1, coord[1])
            adj_1 = (coord[0], coord[1] + 1)
            adj_2 = (coord[0] - 1, coord[1])
            adj_3 = (coord[0], coord[1] - 1)
            adj_coords = [adj_0, adj_1, adj_2, adj_3]

            # check adjacent nodes for connections
            for elem in adj_coords:

                if elem in qubit_mapping:
                    # determine position in kronecker
                    index_1 = qubit_mapping.index(elem)

                    # check whether this connection has been made
                    if not set([index_0, index_1]) in visited_connections:
                        # connect coord and elem
                        # should qubit_num be changed to len(qubit_mapping)
                        con_str_x = ['i']*len(qubit_mapping) # qubit_num
                        con_str_y = ['i']*len(qubit_mapping) # qubit_num
                        con_str_x[index_0] = 'x'
                        con_str_x[index_1] = 'x'
                        con_str_y[index_0] = 'y'
                        con_str_y[index_1] = 'y'

                        op_list_x = [str_dict[ch] for ch in con_str_x]
                        op_list_y = [str_dict[ch] for ch in con_str_y]
                        op_xy     = kron_list(op_list_x) + kron_list(op_list_y)

                        ret_list.append(op_xy)
                        str_list.append(''.join(con_str_x)) # x represents xy
                        amp_list.append(amp2)

                        connection = set([index_0, index_1])
                        visited_connections.append(connection)
                    else:
                        continue
                else:
                    continue
        else:
            continue

    return ret_list, str_list, amp_list

def jj_hamiltonian_2D(qubit_num, qubit_state_num, amp1, amp2, qubit_mapping):

    # mapping is in the form of a list of tuples
    Q_x = np.array([[0.0, 1.0],[1.0,0.0]], dtype = np.complex)
    Q_y = np.array([[0.0, 1.0j],[-1.0j,0.0]], dtype = np.complex)
    Q_z = np.array([[1.0, 0.0],[0.0,-1.0]], dtype = np.complex)
    I_q = np.identity(2)
    ret_list, str_list, amp_list = [], [],[]
    # note the use of 'x' for 'xy' coupling
    str_dict = {'i':I_q, 'x':Q_x, 'y':Q_y, 'z':Q_z}
    
    # keep track of connections already made
    # improvement comes through further using hash tables
    visited_connections = []
    visited_nodes       = {qubit_mapping[i]: 0 for i in range(len(qubit_mapping))}

    for index_0 in range(len(qubit_mapping)):

        coord = qubit_mapping[index_0]

        if not visited_nodes[coord] == 1:

            visited_nodes[coord] = 1
            adj_0 = (coord[0] + 1, coord[1])
            adj_1 = (coord[0], coord[1] + 1)
            adj_2 = (coord[0] - 1, coord[1])
            adj_3 = (coord[0], coord[1] - 1)
            adj_coords = [adj_0, adj_1, adj_2, adj_3]

            # check adjacent nodes for connections
            for elem in adj_coords:

                if elem in qubit_mapping:
                    # determine position in kronecker
                    index_1 = qubit_mapping.index(elem)

                    # check whether this connection has been made
                    if not set([index_0, index_1]) in visited_connections:
                        # connect coord and elem
                        # should qubit_num be changed to len(qubit_mapping)
                        con_str_x = ['i']*len(qubit_mapping) # qubit_num
                        con_str_y = ['i']*len(qubit_mapping) # qubit_num
                        con_str_z = ['z']*len(qubit_mapping) # qubit_num
                        con_str_x[index_0] = 'x'
                        con_str_x[index_1] = 'x'
                        con_str_y[index_0] = 'y'
                        con_str_y[index_1] = 'y'
                        con_str_z[index_0] = 'z'
                        con_str_z[index_1] = 'z'

                        op_list_x = [str_dict[ch] for ch in con_str_x]
                        op_list_y = [str_dict[ch] for ch in con_str_y]
                        op_list_z = [str_dict[ch] for ch in con_str_z]
                        op_xyz     = kron_list(op_list_x) + kron_list(op_list_y) + kron_list(op_list_z)

                        ret_list.append(op_xyz)
                        str_list.append(''.join(con_str_x)) # x represents xyz
                        amp_list.append(amp1)

                        connection = set([index_0, index_1])
                        visited_connections.append(connection)
                    else:
                        continue
                else:
                    continue
        else:
            continue

    return ret_list, str_list, amp_list

def zz_hamiltonian_2D(qubit_num, qubit_state_num, amp1, amp2, qubit_mapping):

    # mapping is in the form of a list of tuples
    Q_x = np.array([[0.0, 1.0],[1.0,0.0]], dtype = np.complex)
    Q_y = np.array([[0.0, 1.0j],[-1.0j,0.0]], dtype = np.complex)
    Q_z = np.array([[1.0, 0.0],[0.0,-1.0]], dtype = np.complex)
    I_q = np.identity(2)
    ret_list, str_list, amp_list = [], [],[]
    # note the use of 'x' for 'xy' coupling
    str_dict = {'i':I_q, 'x':Q_x, 'y':Q_y, 'z':Q_z}
    
    # keep track of connections already made
    # improvement comes through further using hash tables
    visited_connections = []
    visited_nodes       = {qubit_mapping[i]: 0 for i in range(len(qubit_mapping))}

    for index_0 in range(len(qubit_mapping)):

        coord = qubit_mapping[index_0]

        if not visited_nodes[coord] == 1:

            visited_nodes[coord] = 1
            adj_0 = (coord[0] + 1, coord[1])
            adj_1 = (coord[0], coord[1] + 1)
            adj_2 = (coord[0] - 1, coord[1])
            adj_3 = (coord[0], coord[1] - 1)
            adj_coords = [adj_0, adj_1, adj_2, adj_3]

            # check adjacent nodes for connections
            for elem in adj_coords:

                if elem in qubit_mapping:
                    # determine position in kronecker
                    index_1 = qubit_mapping.index(elem)

                    # check whether this connection has been made
                    if not set([index_0, index_1]) in visited_connections:
                        # connect coord and elem
                        # should qubit_num be changed to len(qubit_mapping)
                        con_str_z = ['i']*len(qubit_mapping) # qubit_num
                        con_str_z[index_0] = 'z'
                        con_str_z[index_1] = 'z'

                        op_list_z = [str_dict[ch] for ch in con_str_z]
                        op_z     = kron_list(op_list_z)

                        ret_list.append(op_z)
                        str_list.append(''.join(con_str_z))
                        amp_list.append(amp1)

                        connection = set([index_0, index_1])
                        visited_connections.append(connection)
                    else:
                        continue
                else:
                    continue
        else:
            continue

    return ret_list, str_list, amp_list

def sc_Hamiltonian(qubit_num, qubit_state_num, amp1, amp2):
    
    Q_x = np.array([[0.0, 1.0],[1.0,0.0]], dtype = np.complex)
    Q_y = np.array([[0.0, 1.0j],[-1.0j,0.0]], dtype = np.complex)
    I_q = np.identity(2)
    ret_list, str_list, amp_list = [], [],[]
    str_dict = {'i':I_q, 'x':Q_x+Q_y, 'y':Q_y}
    for i in range(qubit_num):
        str_x, str_y  = ['i']*qubit_num, ['i']*qubit_num
        str_x[i],str_y[i] = 'x','y'
        ind_x, ind_y = [str_dict[st] for st in str_x], [str_dict[st] for st in str_y]

        ret_list.append(kron_list(ind_x))
        ret_list.append(kron_list(ind_y))
        str_list.append(''.join(str_x))
        str_list.append(''.join(str_y))
        
        amp_list.append(amp1)
        amp_list.append(amp1)

        for j in range(i, qubit_num):
            if i!=j: 
                str_xx  = ['i']*qubit_num
                str_xx[i],str_xx[j]= 'x','x'
                ind_xx = [str_dict[st] for st in str_xx]

                ret_list.append(kron_list(ind_xx))
                str_list.append(''.join(str_xx))

                amp_list.append(amp2)

    return ret_list, str_list, amp_list
