import numpy as np
import os
import time

from config import vehicle
from dynamic_st import Dynamic_ST
from kinematic_st import Kinematic


# Load data (states & inputs)
def load_data(file_path, idx, models):
    if idx == -1:
        data     = np.loadtxt(file_path, delimiter = ',',comments = '#')
    else:
        data     = np.loadtxt(file_path, delimiter = ',',comments = '#',max_rows=idx)
    wheel_v  = data[:,10:14] * 0.277778		   # wheel speed at each time step (should be a list or array)
    states   = data[:,1:9]
    inputs   = data[:,8:10]                    # ax(m/s^2),deltadelta(rad/s)
    # gpk_inputs  = np.zeros(data[:,8:11].shape)
    ekin_inputs = np.zeros(data[:,7:11].shape)
    roll   = data[:,14] #* -1	
    times  = data[:,0]
    states = np.stack((states[:,0],states[:,1],states[:,2],states[:,3],states[:,4],states[:,5],states[:,6],states[:,7],wheel_v[:,0],wheel_v[:,1],wheel_v[:,2],wheel_v[:,3],roll),axis=1)			   # Adding roll map data
    # x(m),y(m),vx(m/s),vy(m/s),phi(rad),delta(rad),omega(rad/s),ax(m/s^2),deltadelta(rad/s),wheel_fl/fr/ll/lr
    return times, states, inputs, ekin_inputs # gpk_inputs, 



    
    
    



