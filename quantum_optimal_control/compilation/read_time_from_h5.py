import re
import os
import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys 

# Change this path
path = '/scratch/midway2/yunong/quantum_optimal_control/pulses/outputpulses/'

def read_block_duration(filename):

    with h5py.File(path+filename,'r') as hf:  
        if hf.get('uks') != None:
            speed_up_uks = np.array(hf.get('uks'))[-1]/(2*np.pi)
            speed_up_traj = np.array(hf.get('inter_vecs_mag_squared'))[-1]
            
            inter_vecs_real = np.array(hf.get('inter_vecs_real'))[-1]
            inter_vecs_imag = np.array(hf.get('inter_vecs_imag'))[-1]
            print "shapes"
            total_time = np.array(hf.get('total_time'))
            
            err = np.array(hf.get('error'))
            print 'error:', err[-1]
            target_U = np.array(hf.get('U'))
            print("total time:", total_time)
            if(total_time != None):
                
                time = np.linspace( float(total_time)/inter_vecs_real.shape[-1],float(total_time),inter_vecs_real.shape[-1])
                
                psi = inter_vecs_real + 1j * inter_vecs_imag
                
                all_U_arr = []
               
                for ii in range(psi.shape[0]): 
                    init_vec = np.zeros(psi.shape[0])
                    init_vec[ii] = 1
                    
                    U_list = []
                    for jj in range(psi.shape[-1]):
                        U_jj = np.outer(psi[ii,:,jj],init_vec)
                    
                        U_list.append(U_jj)
                        
                    all_U_arr.append(np.array(U_list))


                all_U = np.sum(all_U_arr,axis=0)

                all_U.shape

                fidelity_tlist = []

                for ii in range(all_U.shape[0]):
                    fidelity = (abs(np.trace(np.dot(np.transpose(np.conjugate(target_U)),all_U[ii]))))/all_U.shape[1]
                    
                    fidelity_tlist.append(fidelity)

                target_fidelity = 0.999
                fidelity_tarr = np.array(fidelity_tlist)
                target_reached_tarr = fidelity_tarr > 0.999

                target_reached_tarr_ind = [i for i, x in enumerate(target_reached_tarr) if x]
                #print(inp)
                if(len(target_reached_tarr_ind)>0):
                    target_reached_shortest_t = time[target_reached_tarr_ind[0]]
                    print "shortest t:", target_reached_shortest_t 
                else: print "not converged"

    return target_reached_shortest_t
