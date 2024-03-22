"""	
Dynamic single track model.
x is a 6x1 vector: [x, y, vx, vy, phi, 𝛿, ω]^T
u is a 4x1 vector: [ax, Δ𝛿]^T
""" 

import numpy as np
import math
from config import vehicle

class Dynamic_ST:
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
		self.Csf    = vehicle_params['Csf']
		self.Csr    = vehicle_params['Csr']
		self.hcog   = vehicle_params['hcog']
		self.mu     = vehicle_params['mu']
		self.min_v  = vehicle_params['min_v']
		self.max_v  = vehicle_params['max_v']
		self.h_aero = vehicle_params['haero']
		self.tire_p = vehicle_params['l_pressure']
		# Physical constraints
		self.switch_v    = 4.0
		self.max_acc     = vehicle_params['max_acc']
		self.min_acc     = vehicle_params['min_acc']
		self.max_steer   = vehicle_params['max_steer']
		self.min_steer   = vehicle_params['min_steer']
		self.max_steer_v = vehicle_params['max_steer_vel']
		self.min_steer_v = vehicle_params['min_steer_vel']
		self.width       = vehicle_params['width']
		self.length      = vehicle_params['length']
		self.MF_long     = vehicle_params['MF_long']
		self.MF_lat	 	 = vehicle_params['MF_lat']
		self.n_states    = 7
		self.n_inputs    = 2
		# Force coefficients
		self.rho   = 1.225							     # Mass density of air [kg/m^3]
		self.Af    = 1.6 + 0.00056 * (self.mass - 765)   # Frontal area  (ref: Book VDC Eq 4.3)
		self.Cd    = 0.725 								 # Aerodynamic drag coefficient: 0.7 to 1.1 typical values for Formula Once car
		self.theta = 0.21								 # Banked curve	[rad]
		self.tire_prs   = 200568.49

	def sim_continuous(self, x0, u, wheel_v, t, roll=None):
		"""	simulates the nonlinear continuous model with given input vector
			by numerical integration using 4th order Runge Kutta method
		"""
		n_steps    = len(t)-1 
		x          = np.zeros([n_steps+1, self.n_states]) 
		dxdt       = np.zeros([n_steps+1, self.n_states]) 
		dxdt[0, :] = self.derivative_eqs(None, x0, wheel_v, u[:,0], roll)
		x[0, :]    = x0
		for ids in range(1, n_steps+1):
			x[ids, :]  = self.odeintRK4(x[ids-1, :], [t[ids-1],t[ids]], u[:,ids-1], wheel_v, roll)
			dxdt[ids, :] = self.derivative_eqs(None, x[ids, :], u[:,ids-1], wheel_v, roll)
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
			x[ids, :]  = self.odeintRK4(x[ids-1, :], [t[0],t[1]], u[ids-1, :], wheel_v[ids-1], roll[ids-1])
			dxdt[ids, :] = self.derivative_eqs(None, x[ids-1, :], u[ids-1, :], wheel_v[ids-1], roll[ids-1])
			step_ls.append(step)
		return x, dxdt, step_ls

	def calc_forces(self, x, u, wheel_v, roll):
		acc   = u[0] 
		vx    = x[2] 
		vy    = x[3]
		psi   = x[4]
		𝛿     = x[5] 
		ω     = x[6] 
		froll = 0.029					# rolling resistance coefficient f varies in the range 0.01 to 0.04. (Wong, 2001)
		'''
		Rolling resistance force = Cr * Normal tire force
		Cr = coefficient of rolling resistance,
		Normal tire force = Weight of the vehicle * (g + a * sin(θ)) / number of tires
		g = gravitational acceleration (9.8 m/s^2 on the surface of the Earth)
		a = lateral acceleration due to the centripetal force
		θ = banking angle
		'''
		Fznorm = self.mass * (self.g + acc * math.sin(roll)) / 4
		Rx     = froll * Fznorm
		#######################################
		# Calculate forces on tires
		#######################################
		MF_params_lat = self.MF_lat
		Bf = MF_params_lat['Bf']
		Cf = MF_params_lat['Cf']
		Df = MF_params_lat['Df']
		Ef = MF_params_lat['Ef']
		Br = MF_params_lat['Br']
		Cr = MF_params_lat['Cr']
		Dr = MF_params_lat['Dr']
		Er = MF_params_lat['Er']
		Shfy = MF_params_lat['Shfy']
		Shry = MF_params_lat['Shry']
		Svfy = MF_params_lat['Svfy']
		Svry = MF_params_lat['Svry']
		MF_params_long = self.MF_long
		long_B = MF_params_long['B']
		long_C = MF_params_long['C']
		long_D = MF_params_long['D']
		long_E = MF_params_long['E']
		# Slip angles & slip ratio
		kappa, alpha_f, alpha_r = self.calc_slips(x, u, wheel_v)
		alpha_f += Shfy
		alpha_r += Shry
		# Lateral tire forces
		Fyf = Svfy + Df * math.sin(Cf * math.atan(Bf * alpha_f) - Ef * (Bf * alpha_f - math.atan(Bf * alpha_f)))
		Fyr = Svry + Dr * math.sin(Cr * math.atan(Br * alpha_r) - Er * (Br * alpha_r - math.atan(Br * alpha_r)))
		# Longitudinal tire forces (Assume RWD)
		Fxr = long_D * np.sin(long_C * np.arctan(long_B * kappa - long_E * (long_B * kappa - np.arctan(long_B * kappa))))
		Fcx = self.mass * vy * ω							 		
		Fbx = self.mass * self.g * math.sin(math.radians(roll)) * math.sin(psi)
		# Aerodynamic drag force
		Faero = 1/2 * self.rho * self.Cd * self.Af * vx ** 2
		# Longitudinal total forces
		Fxl = Fxr - Fcx - Fbx - Fyf * np.sin(𝛿)
		# Lateral force total
		Fby = self.mass * self.g * math.sin(math.radians(roll)) * math.cos(psi)     
		Fcy = self.mass * vx * ω							  
		Fyl = Fyr + Fyf * math.cos(𝛿) - Fby - Fcy
		return Fyr, Fyf, Fby, Fcy, Fxr

	def calc_slips(self, x, u, wheel_v):
		acc   = u[0] 
		vx    = x[2] 
		vy    = x[3]
		𝛿     = x[5] 
		ω     = x[6] 
		wheel_v += 1e-12 
		ACCELERATION = True if acc > 0 else False
		# Longitudinal slip ratio is defined differently during acceleration and braking
		if ACCELERATION:
			kappa = (wheel_v - vx) / wheel_v    			# Eq. 4.11
		else:
			kappa = (wheel_v - vx) / vx						# Eq. 4.12
		# Slip angles of front wheel and rear wheel
		β = math.atan(self.lr*np.tan(𝛿)/(self.lf + self.lr))
		theta_vf = math.atan((vy + self.lf * ω) / abs(vx))          # Eq. 2.29  
		theta_vr = math.atan((vy - self.lr * ω) /  abs(vx))			# Eq. 2.30  
		alpha_f  = 𝛿 - theta_vf					 		         # Eq. 2.23
		alpha_r  = -1 * theta_vr								    # Eq. 2.24
		if vx <= 3.0:
			alpha_f = 0.0
			alpha_r = 0.0
			kappa   = 0.0
		# print("front slip: {:.4f}, rear slip: {:.4f}".format(alpha_f, alpha_r))
		return kappa, alpha_f, alpha_r

	def derivative_eqs(self, t, x, u, wheel_v, roll):
		ax   = u[0]
		𝛿_vel = u[1]
		vx    = x[2] 
		vy    = x[3]
		phi   = x[4]
		𝛿     = x[5] 
		ω     = x[6] 

		Fyr, Fyf, _, _, Fxr = self.calc_forces(x, u, wheel_v, roll)
		size    = len(x)
		dxdt    = np.zeros(size)
		dxdt[0] = vx * math.cos(phi) - vy * math.sin(phi)
		dxdt[1] = vx * math.sin(phi) + vy * math.cos(phi)
		dxdt[2] = 1./self.mass * (Fxr - Fyf*math.sin(𝛿)) + vy * ω
		dxdt[3] = 1./self.mass * (Fyr + Fyf * math.cos(𝛿)) - vx * ω 
		dxdt[4] = ω 
		dxdt[5] = 𝛿_vel 
		dxdt[6] = 1/self.Iz * (self.lf * Fyf * math.cos(𝛿) - self.lr * Fyr)
		if vx <= self.switch_v: # switch to simple model 
			if vx <= -0.01:
				ax = 0
			dxdt[0] = vx * math.cos(phi) - vy * math.sin(phi)
			dxdt[1] = vx * math.sin(phi) + vy * math.cos(phi)
			dxdt[2] = ax
			dxdt[3] = ax * np.tan(𝛿)
			dxdt[4] = ω 
			dxdt[5] = 𝛿_vel 
			dxdt[6] = (1 / (self.lf + self.lr)) * ax * ω 
		
		return dxdt
		
	# Numerical integration methods 
	# 3/8-rule, variation of Runge-Kutta 4th Order Method

	def odeintRK4(self, y0, t, u, wheel_v, roll):
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
			K[0]   = h * fun(t[i], y0, u, wheel_v, roll)
			K[1]   = h * fun(C[0] * h, y0+A[0][0]*K[0], u, wheel_v, roll)
			K[2]   = h * fun(C[1] * h, y0+A[1][0]*K[0]+A[1][1]*K[1], u, wheel_v, roll)
			K[3]   = h * fun(C[2] * h, y0+A[2][0]*K[0]+A[2][1]*K[1]+A[2][2]*K[2], u, wheel_v, roll)
			y_next[i, :] = y0 + B@K
			y0 = y_next[i, :]
		return y_next

