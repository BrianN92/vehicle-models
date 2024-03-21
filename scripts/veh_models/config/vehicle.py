def Vehicle_params():
	mu = 1.0489 		    # Surface friction coefficient
	Csf = 4.718				# Cornering stiffness coefficient, front [1/rad]
	Csr = 5.4562			# Cornering stiffness coefficient, rear [1/rad]
	lf = 1.7    			# Distance from center of gravity to front axle [m]
	lr = 1.2 			    # Distance from center of gravity to rear axle [m]
	hcog = 0.4		    	# Height of center of gravity [m]
	haero = 0.25			# height of the location at which the equivalent aerodynamic force acts 
	mass = 800      		# Total mass of the vehicle [kg]
	Iz = 1000   			# Moment of inertial of the entire vehicle about the z axis [kgm^2]
	width = 2.2		    	# width of the vehicle [m]
	length = 5.1            # length of the vehicle [m]
	g = 9.81  				# Gravitational acceleration [m/s^2]
	p = 151685              # Tire pressure of left front and rear tires (Pa)
	min_steer = -0.3 	    # Minimum steering angle constraint [rad]
	max_steer = 0.3		    # Maximum steering angle constraint [rad]
	max_steer_vel = 3   	# Maximum steering velocity constraint [rad/s]
	min_steer_vel = -3   	# Minimum steering velocity constraint [rad/s]
	min_v = -5.0			# Minimum longitudinal velocity [m/s]
	max_v = 89.408			# Maximum longitudinal velocity [m/s]
	max_acc = 9.51 			# max acceleration [m/s^2]
	min_acc = -13.26 		# max deceleration [m/s^2]
	#######################
	# Tire model parameters
	#######################
	MF_Long = {
		'B': 18,
		'C': 2,
		'D': 1.7,
		'E': 1
	}
	MF_Lat = {
		'Bf': 17,
		'Cf': 1.67,
		'Df': 2651,
		'Ef': 0.2,
		'Br': 17.4,
		'Cr': 1.7,
		'Dr': 3452,
		'Er': -0.08,
		'Shfy': 0.007,
		'Shry': 0.01,
		'Svfy': 134.7,
		'Svry': 194.07
	}
	tire_re = 0.3118		# radius [m]
	max_inputs = [max_acc, max_steer]
	min_inputs = [min_acc, min_steer]

	max_rates = [None, max_steer_vel]
	min_rates = [None, min_steer_vel]	

	params_dict = {
		'lf': lf,
		'lr': lr,
		'mass': mass,
		'Iz': Iz,
		'Csf': Csf,
		'Csr': Csr,
		'hcog': hcog,
		'haero':haero,
		'mu': mu,
		'min_v': min_v,
		'max_v': max_v,
		# 'switch_v': switch_v,
		'max_acc': max_acc,
		'min_acc': min_acc,
		'max_steer': max_steer,
		'min_steer': min_steer,
		'max_steer_vel': max_steer_vel,
		'min_steer_vel': min_steer_vel,
		'width': width,
		'length': length,
		'max_inputs': max_inputs,
		'min_inputs': min_inputs,
		'max_rates': max_rates,
		'min_rates': min_rates,
		'g': g,
		'l_pressure': p, 
		'MF_long': MF_Long,
		'MF_lat': MF_Lat
		}
	return params_dict
