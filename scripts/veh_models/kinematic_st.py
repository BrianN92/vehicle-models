"""	
Kinematic single track model.
x is a 5x1 vector: [x, y, psi, v, 𝛿]^T
u is a 2x1 vector: [acc,Δ𝛿]^T
""" 

import numpy as np
from config import vehicle


class Kinematic:
	def __init__(self, params = None):
		if params == None:
			vehicle_params = vehicle.Vehicle_params()
		else:
			vehicle_params = params
		self.g      = vehicle_params['g']
		self.lf     = vehicle_params['lf']
		self.lr     = vehicle_params['lr']
		self.mass   = vehicle_params['mass']
		self.Iz     = vehicle_params['Iz']
		self.hcog   = vehicle_params['hcog']
		self.mu     = vehicle_params['mu']
		self.min_v  = vehicle_params['min_v']
		self.max_v  = vehicle_params['max_v']
		self.h_aero = vehicle_params['haero']
		self.tire_p = vehicle_params['l_pressure']
		# Physical constraints
		self.max_acc     = vehicle_params['max_acc']
		self.min_acc     = vehicle_params['min_acc']
		self.max_steer   = vehicle_params['max_steer']
		self.min_steer   = vehicle_params['min_steer']
		self.max_steer_v = vehicle_params['max_steer_vel']
		self.min_steer_v = vehicle_params['min_steer_vel']
		self.width       = vehicle_params['width']
		self.length      = vehicle_params['length']
		self.n_states    = 5
		self.n_inputs    = 2

	def sim_continuous(self, x0, u, t, wheel_v=None, roll=None):
		"""	simulates the nonlinear continuous model with given input vector
			by numerical integration using 6th order Runge Kutta method
			x0 is the initial state of size 5x1 [x, y, psi, v, 𝛿]^T
			u is the input vector of size 2x1 [acc (or ΔT),Δ𝛿]^T
			t is the time vector of size 2x1 [0, Ts]^T
		"""
		n_steps    = len(t)-1 # 1
		x          = np.zeros([n_steps+1, self.n_states]) # size: 2x5
		dxdt       = np.zeros([n_steps+1, self.n_states]) # size: 2x5
		dxdt[0, :] = self.derivative_eqs(None, x0, u[:,0])
		x[0, :]    = x0
		for ids in range(1, n_steps+1):
			x[ids, :]  = self.odeintRK4(x[ids-1, :], [t[ids-1],t[ids]], u[:,0])
			# if using odeintRK6 method: x[ids, :]  = self.odeintRK6(x[ids-1, :], t[ids], u)
			dxdt[ids, :] = self.derivative_eqs(None, x[ids, :], u[:,0])
		
		return x, dxdt

	def sim_continuous_multistep(self, x0, u, t, step, step_ls, wheel_v = None, roll=None):
		"""	simulates the nonlinear continuous model with given input vector
			by numerical integration using Runge Kutta method
			propagating sim for n steps
		"""
		x          = np.zeros([step+1, self.n_states]) # size: (step+1)x6
		dxdt       = np.zeros([step+1, self.n_states]) # size: (step+1)x6
		x[0, :]    = x0
		for ids in range(1, step+1):
			# x[ids, :]  = self.odeintScipy(x[ids-1, :], [t[0],t[1]], u[ids-1, :])
			x[ids, :]  = self.odeintRK4(x[ids-1, :], [t[0],t[1]], u[ids-1, :])
			dxdt[ids, :] = self.derivative_eqs(None, x[ids-1, :], u[ids-1, :])
			step_ls.append(step)
		return x, dxdt, step_ls

	def derivative_eqs(self, t, x, u):
		steer_v = u[1]
		acc     = u[0]
		psi     = x[2]
		v       = x[3] 
		𝛿       = x[4]
		β 		= np.arctan(self.lr*np.tan(𝛿)/(self.lf + self.lr))
		size    = len(x)
		dxdt    = np.zeros(size)
		dxdt[0] = v * np.cos(psi + β)
		dxdt[1] = v * np.sin(psi + β)
		dxdt[2] = np.sin(β) * v / self.lr
		dxdt[3] = acc
		dxdt[4] = steer_v 
		return dxdt
		
	# Numerical integration methods 
    # 3/8-rule, variation of Runge-Kutta 4th Order Method
	def odeintRK4(self, y0, t, u):
		A   = np.asarray([ [1/3],
						   [-1/3, 1],
						   [1, -1, 1]], dtype=object)
		B   = np.asanyarray([1/8, 3/8, 3/8, 1/8])
		C   = np.asarray([1/3, 2/3, 1])
		fun = self.derivative_eqs
		y_next = np.zeros([len(t)-1, len(y0)])
		K      = np.zeros((len(B), len(y0)))
		for i in range(len(t)-1):
			h      = t[i+1] - t[i]
			K[0]   = h * fun(t[i], y0, u)
			K[1]   = h * fun(C[0] * h, y0+A[0][0]*K[0], u)
			K[2]   = h * fun(C[1] * h, y0+A[1][0]*K[0]+A[1][1]*K[1], u)
			K[3]   = h * fun(C[2] * h, y0+A[2][0]*K[0]+A[2][1]*K[1]+A[2][2]*K[2], u)
			y_next[i, :] = y0 + B@K
			y0 = y_next[i, :]
		return y_next

