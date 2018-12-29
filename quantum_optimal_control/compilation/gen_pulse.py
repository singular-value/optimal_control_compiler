import numpy as np
import scipy.linalg as la
import os,sys,inspect, re
import argparse
import itertools
from qutip import *

data_path = '/scratch/midway2/yunong/quantum_optimal_control/pulses/outputpulses/'

from quantum_optimal_control.helper_functions.grape_functions import *
from quantum_optimal_control.main_grape.grape import *
from quantum_optimal_control.core.convergence import *
from quantum_optimal_control.core.analysis import *
from quantum_optimal_control.core.run_session import *
from quantum_optimal_control.main_grape.grape import *
from quantum_optimal_control.compilation.hamiltonian import *
from quantum_optimal_control.compilation.custgate import *
from quantum_optimal_control.compilation.qcircuit import *
from quantum_optimal_control.compilation.qgate import *
from quantum_optimal_control.compilation.gcirc import *

def run_circuit(qasm_file_name, qubit_num, final_unitary, mapping, time, if_qutip = False):
    #Defining time scales
    total_time = time
    steps = qubit_num * 1000 

    # Choose optimizing State transfer or Unitary gate
    state_transfer = False
    
    # Choose whether include intermediate state evolution as part of the graph optimizatio
    use_inter_vecs = True
    
    # Defining U0 (Initial)
    q_identity = np.identity(2**qubit_num)
    U0= q_identity
    
    # Defining dressed info
    is_dressed = False
#     w_c, v_c, dressed_id = get_dressed_info(H0)
#     dressed_info = {'dressed_id':dressed_id, 'eigenvectors':v_c, 'eigenvalues':w_c,'is_dressed':is_dressed}
    
    dressed_info = None
    
    #Defining Concerned states (starting states)

    psi0 =  concerned(qubit_num,2)
    
    #Defining states to include in the drawing of occupation
    states_draw_list = range(2**qubit_num)
    states_draw_names =[]
    for ii in states_draw_list:
        states_draw_names.append(Basis(ii,qubit_num, 2))
    
    #
    H0 = np.identity(2**qubit_num)
    # Defining final U
    SWAP = np.array([[1.0,0.0,0.0,0.0],[0.0,0.0,1.0,0.0],[0.0,1.0,0.0,0.0],[0.0,0.0,0.0,1.0]], dtype=np.complex)
    U = final_unitary
#    print U
    Hops = []
    Hnames = []
    ops_max_amp = []
    Hops, Hnames, ops_max_amp = xy_hamiltonian_2D(qubit_num, 2, 5, 1, mapping)

    #Hops, Hnames, ops_max_amp = two_by_three(a, 0.2, 0.05)
    
    #Defining convergence parameters
    max_iterations = 2000
    decay = 1000 #max_iterations/2
    convergence = {'rate':0.01, 'update_step':10 ,'max_iterations':max_iterations,\
                   'conv_target':1e-3,'learning_rate_decay':decay,'min_grad':1e-60}

    reg_coeffs = {'speedup':100}
    
    #import h5py
    #with h5py.File('/home/nelson/Simulations/GRAPE-Data/spin_chain_hadamard/00002_gpu_spin_chain_hadamard_N9.h5','r') as hf:
    #    u0 = np.array(hf.get('uks'))[-1]
    
    qasm_file_name = qasm_file_name
    file_name = qasm_file_name + 'time'+'-'+time+'qutip'+str(if_qutip)
    
    u0 = None
    uks,U_final = Grape(H0,Hops,Hnames,U,total_time,steps,psi0,convergence=convergence, 
                        draw = [states_draw_list,states_draw_names], show_plots = True, 
                        use_gpu = True,sparse_H=False,sparse_U=False,state_transfer = state_transfer, 
                        use_inter_vecs = use_inter_vecs, unitary_error = 1e-8, 
                        maxA=ops_max_amp,Taylor_terms =[20,0],initial_guess=u0,
                        dressed_info = dressed_info, method = 'ADAM', reg_coeffs=reg_coeffs,
                        file_name=file_name,  
                        data_path = data_path)


if __name__ == "__main__":
    run_circuit(qubit_num, final_U, mapping, time, if_qutip=qut)
