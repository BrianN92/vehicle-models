"""	GPKIN
Generate data by simulating (kinematic/dynamic) bicycle model.
"""
import math
import numpy as np
from config import vehicle
from dynamic_st import Dynamic_ST
from kinematic_st import Kinematic
from utils import load_data
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle, Circle
from matplotlib.lines import Line2D
from matplotlib.animation import FuncAnimation
from collections import deque
import time
import os

class Vehicle:
    def __init__(self, params):
        self.params = params
        self.true_states = None
        self.pred_states = None
        self.inputs      = None
        self.ekin_inputs = None
        self.roll_data   = None
        self.wheel_data  = None
        self.states      = None

    def gen_data(self, path, dataset, models):
        data_limit = -1
        times, states, self.inputs, self.ekin_inputs = load_data(path+'{}.csv'.format(dataset), idx = data_limit, models=models)
        # ground truth vehicle states
        self.true_states = states[:, :7]  # x(m),y(m),vx(m/s),vy(m/s),psi(rad),delta(rad),omega(rad/s
        # model prediction vehicle states
        self.pred_states = np.zeros(self.true_states.shape)
        self.pred_states[0, :] = self.true_states[0, :]
        self.wheel_data = states[:, -5:-1]
        self.roll_data  = states[:, -1]
        self.states = states

class Simulation:
    def __init__(self, vehicle, model):
        self.vehicle = vehicle
        self.N_SAMPLES = len(vehicle.inputs)
        self.SAMPLING_TIME = 0.04
        self.states = vehicle.states
        self.inputs = vehicle.inputs
        self.true_states = vehicle.true_states
        self.states_pred = vehicle.pred_states
        self.wheel_data = vehicle.wheel_data
        if model != 'single_track':
            self.wheel_data = np.mean(vehicle.wheel_data, axis=1, keepdims=True)
        self.roll_data  = vehicle.roll_data
        self.dxdt_ek = np.zeros(vehicle.true_states.shape)
        self.base_horizon = 10
        self.steer_inputs = vehicle.ekin_inputs[:, 2]
        self.veh_model = model

    def calc_inrange(self, arr, a, b):
        lower_diff = arr - a 
        lower_diff[lower_diff<0] = 0

        upper_diff = b - arr
        upper_diff[upper_diff<0] = 0
        in_range_diff = np.minimum(lower_diff, upper_diff)
        return in_range_diff
    
    # Function to calculate out-range difference
    def calc_steer_diff(self, steering, std_dev):
        # Check if the current steering angle is out of range based on standard deviation
        if abs(steering) > std_dev:
            # Calculate the out range difference between the current steering angle and standard deviation
            out_range_difference = abs(steering) - std_dev
            in_range_difference  = 0
        else:
            # Calculate the in range difference between the current steering angle and standard deviation
            in_range_difference  = std_dev - abs(steering) 
            out_range_difference = 0
        return in_range_difference, out_range_difference
    
    def open_sim(self, multi_step, step_base, c1, c2, horizon_fixed=True):
        # Simulation logic here
        self.base_horizon = step_base
        step = step_base
        if self.veh_model == 'dynamic':
            model = Dynamic_ST()
        elif self.veh_model == 'kinematic':
            model = Kinematic()
            # x: [x, y, psi, v, 𝛿]^T
			# u: [acc, Δ𝛿]^T
            kin_states = np.zeros((len(self.true_states), model.n_states))
            kin_inputs = np.zeros((len(self.steer_inputs), model.n_inputs))
            kin_states[:, :2] = self.states[:,:2]
            kin_states[:, 2] = self.states[:,4]
            kin_states[:, 3] = self.states[:,2]
            kin_states[:, 4] = self.states[:,5]
            kin_inputs[:, 0] = self.states[:,7]
            kin_inputs[1:, 1] = np.diff(self.steer_inputs)
            self.true_states = kin_states
            self.inputs = kin_inputs
            self.states_pred = np.zeros(kin_states.shape)
            self.states_pred[0,:] = kin_states[0, :]
            self.dxdt_ek = np.zeros(kin_states.shape)
        std_steer = np.std(self.steer_inputs)
        c1_scaler = c1
        c2_scaler = c2
        self.step_ls = [step]
        if multi_step:
            # print("Multi-step mode")
            idn = 0
            if horizon_fixed:
                while idn < (self.N_SAMPLES-1):
                    if idn > self.N_SAMPLES - 1 - step:
                        step = self.N_SAMPLES - 1 - idn
                    inputs = np.vstack((self.inputs[idn:idn+step, :]))
                    x_next, dxdt_next, self.step_ls = model.sim_continuous_multistep(self.true_states[idn, :], inputs, [0, self.SAMPLING_TIME], step, self.step_ls, wheel_v=self.wheel_data[idn:idn+step, :] , roll=self.roll_data[idn:idn+step])
                    self.states_pred[idn+1:idn+step+1, :] = x_next[1:, :]
                    self.dxdt_ek[idn+1:idn+step+1, :] = dxdt_next[1:, :]
                    idn += step # step idx update 
            else:
                while idn < (self.N_SAMPLES-1):
                    step  = step_base
                    curr_steer = self.inputs[idn, 2]
                    in_range_diff, out_range_diff = self.calc_steer_diff(curr_steer, std_steer)
                    step += math.floor(in_range_diff*c1_scaler*step) - math.ceil(out_range_diff*c2_scaler*step)
                    if idn > self.N_SAMPLES - 1 - step:
                        step = self.N_SAMPLES - 1 - idn
                    step = 5 if step < 5 else step # Lower bound of step
                    # step_ls.append(step)
                    inputs = np.vstack((self.inputs[idn:idn+step, :]))
                    x_next, dxdt_next, self.step_ls = model.sim_continuous_multistep(self.true_states[idn, :], inputs, [0, self.SAMPLING_TIME], step, self.step_ls, wheel_v=self.wheel_data[idn:idn+step,:] , roll=self.roll_data[idn:idn+step])
                    self.states_pred[idn+1:idn+step+1, :] = x_next[1:, :]
                    self.dxdt_ek[idn+1:idn+step+1, :] = dxdt_next[1:, :]
                    idn += step
        else:
            print("Single-step mode")
            for idn in range(self.N_SAMPLES-1):
                inputs = np.vstack((self.inputs[idn, :], self.inputs[idn+1, :])).T
                x_next, dxdt_next   = model.sim_continuous(self.vehicle.true_states[idn, :], inputs, [0, self.SAMPLING_TIME])
                self.states_pred[idn+1, :] = x_next[-1, :]
                self.dxdt_ek[idn+1, :]   = x_next[-1, :]
    
class Plotting:
    def __init__(self, simulation, vx_start, left_np, right_np, pit_np, Whole_track = False):
        self.simulation = simulation
        model = simulation.veh_model
        horizon = simulation.base_horizon + 2
        self.vx = simulation.states[:, 2]
        self.steering_angle = simulation.steer_inputs
        self.step_ls = simulation.step_ls
        self.id = np.where(self.vx>vx_start)[0][0]

        self.x = simulation.states[self.id:, 0]
        self.y = simulation.states[self.id:, 1]
        self.vx = self.vx[self.id:]
        self.steering_angle = self.steering_angle[self.id:]
        self.x_pred = simulation.states_pred[self.id:, 0]
        self.y_pred = simulation.states_pred[self.id:, 1]
        self.angles = simulation.states[self.id:, 4]
        if model == 'kinematic':
            self.angles_pred = simulation.states_pred[self.id:, 2]
        else:
            self.angles_pred = simulation.states_pred[self.id:, 4]

        # Initialize the figure and axis
        self.fig, self.ax = plt.subplots(figsize=(12,8))
        self.ax.plot(left_np[:,0], left_np[:,1], label='left bound')
        self.ax.plot(right_np[:, 0], right_np[:,1], label='right bound')
        if pit_np.any():
            self.ax.plot(pit_np[:, 0], pit_np[:,1], label='pit')
        self.ax.set_xlabel("$X [m]$", fontsize=14)
        self.ax.set_ylabel("$Y [m]$", fontsize=14)
        self.ax.set_title("Animated Car Movement: Ground Truth vs. Model [{}]".format(model), fontsize=18)
        # Initialize a text label
        # text_label = ax.text(0.1, 0.1, 'Car Speed: ${:.2f} [m/s]$'.format(vx[0]), transform=ax.transAxes)

        self.focus_on = False if Whole_track else True
        if Whole_track:
            self.ax.set_xlim(min(self.x)-10, max(self.x)+10)
            self.ax.set_ylim(min(self.y)-10, max(self.y)+10)
            self.ax.set_aspect('equal', 'box')
        # Initialize arrow storage and maximum number of arrows to keep
        self.arrow1_storage = deque(maxlen=horizon*2)  # keep only arrows during horizons * 2
        self.arrow2_storage = deque(maxlen=horizon*2)


        # Initialize arrow (vehicle)
        self.length = 0.5
        arrow = plt.Arrow(self.x[0], self.y[0], self.length * np.cos(self.angles[0]), self.length * np.sin(self.angles[0]), width=0.1, color='green', label='True position')
        arrow_pred = plt.Arrow(self.x_pred[0], self.y_pred[0], self.length * np.cos(self.angles_pred[0]), self.length * np.sin(self.angles_pred[0]), width=0.1, color='red', label='Ekin position')
        self.ax.add_patch(arrow)
        self.ax.add_patch(arrow_pred)
        self.ax.legend()

        # Initialize car patches
        self.true_patches = self.draw_car(self.x[0], self.y[0], self.angles[0], carID='true')
        self.pred_patches = self.draw_car(self.x[0], self.y[0], self.angles[0], carID='pred')
        for true_patch, pred_patch in zip(self.true_patches,self.pred_patches):
            self.ax.add_patch(pred_patch)
            self.ax.add_patch(true_patch)
            
        self.last_time = None

        # Initialize steering wheel
        self.steering_wheel_radius = 15.0  # Adjust as needed
        self.steering_wheel_offset_x = 0  # Horizontal offset
        self.steering_wheel_offset_y = 35  # 10 units above the car
        # Define the offset from the steering wheel's position
        self.text_offset_x = 0  # No horizontal offset
        self.text_offset_y = -self.steering_wheel_radius - 10  # 5 units below the steering wheel
        self.draw_steering_wheel(self.ax, self.x[0] + self.steering_wheel_offset_x, self.y[0] + self.steering_wheel_offset_y, self.steering_wheel_radius, self.angles[0])
        self.steering_angle_text = self.ax.text(self.x[0] + self.steering_wheel_offset_x, self.y[0] + self.steering_wheel_offset_y - self.steering_wheel_radius - 10, '', ha='center')


    def draw_car(self, x, y, angle, carID):
        patches = []
    
        # Define the basic body of the car as a triangle
        body_x = [-2, -2, 2.5]
        body_y = [-1., 1., 0]
        
        # Define the front wing
        front_wing_x = [1.5, 2, 2, 1.5]
        front_wing_y = [-1.2, -1.2, 1.2, 1.2]
        
        # Define the rear wing
        rear_wing_x = [-2.5, -2, -2, -2.5]
        rear_wing_y = [-0.7, -0.7, 0.7, 0.7]

        # Tires (drawn as circles; you may replace with rectangles or other shapes)
        tire_positions = [(-1.5, 0.9), (-1.5, -0.9), (1.2, 0.4), (1.2, -0.4)]

        # Rotate and translate
        angle_rad = angle
        rot_matrix = np.array([
            [np.cos(angle_rad), -np.sin(angle_rad)],
            [np.sin(angle_rad), np.cos(angle_rad)]
        ])
        
        body = np.dot(rot_matrix, np.array([body_x, body_y]))
        front_wing = np.dot(rot_matrix, np.array([front_wing_x, front_wing_y]))
        rear_wing = np.dot(rot_matrix, np.array([rear_wing_x, rear_wing_y]))
        
        # Translate
        body[0, :] += x
        body[1, :] += y
        front_wing[0, :] += x
        front_wing[1, :] += y
        rear_wing[0, :] += x
        rear_wing[1, :] += y

        tire_positions = np.dot(rot_matrix, np.array(tire_positions).T)
        tire_positions[0, :] += x
        tire_positions[1, :] += y
        
        if carID == 'true':
            patches.append(Polygon(body.T, closed=True, facecolor='green', edgecolor='black'))
            patches.append(Polygon(front_wing.T, closed=True, facecolor='orange', edgecolor='black'))
            patches.append(Polygon(rear_wing.T, closed=True, facecolor='orange', edgecolor='black'))
        elif carID == 'pred':
            patches.append(Polygon(body.T, closed=True, facecolor='red', edgecolor='black'))
            patches.append(Polygon(front_wing.T, closed=True, facecolor='blue', edgecolor='black'))
            patches.append(Polygon(rear_wing.T, closed=True, facecolor='blue', edgecolor='black'))
        
        for pos in tire_positions.T:
            patches.append(Circle((pos[0], pos[1]), radius=0.2, facecolor='black'))

        return patches

    def draw_steering_wheel(self, ax, center_x, center_y, radius, angle):
        # Draw outer circle
        outer_circle = Circle((center_x, center_y), radius, fill=False, color='blue')
        outer_circle.is_steering_wheel_component = True  # Custom attribute
        ax.add_patch(outer_circle)
        
        # Draw inner circle
        inner_circle = Circle((center_x, center_y), radius * 0.2, fill=True, color='blue')
        inner_circle.is_steering_wheel_component = True  # Custom attribute
        ax.add_patch(inner_circle)
        
        # Draw horizontal line (spoke), rotated by the steering angle
        angle_rad = np.radians(angle)
        dx = radius * np.cos(angle_rad)
        dy = radius * np.sin(angle_rad)
        spoke = Line2D([center_x - dx, center_x + dx], [center_y - dy, center_y + dy], color='blue')
        spoke.is_steering_wheel_component = True  # Custom attribute
        ax.add_line(spoke)

    def update(self, frame):
        # Remove the oldest path for ground truth
        if len(self.arrow1_storage) == self.arrow1_storage.maxlen:
            oldest_arrow1 = self.arrow1_storage.popleft()  
            oldest_arrow1.remove()
        # Remove the oldest path for ekin
        if len(self.arrow2_storage) == self.arrow2_storage.maxlen:
            oldest_arrow2 = self.arrow2_storage.popleft()  
            oldest_arrow2.remove()
        # Path of ground truth
        dx = self.length * np.cos(self.angles[frame])
        dy = self.length * np.sin(self.angles[frame])
        arrow = plt.Arrow(self.x[frame], self.y[frame], dx, dy, width=0.1, color='green')
        self.ax.add_patch(arrow)
        # Path of model prediction
        dx_pred = self.length * np.cos(self.angles_pred[frame])
        dy_pred = self.length * np.sin(self.angles_pred[frame])
        arrow_pred = plt.Arrow(self.x_pred[frame], self.y_pred[frame], dx_pred, dy_pred, width=0.1, color='red')
        self.ax.add_patch(arrow_pred)
        self.arrow1_storage.append(arrow)
        self.arrow2_storage.append(arrow_pred)
        # Create and add new car position
        new_patches_pred = self.draw_car(self.x_pred[frame], self.y_pred[frame], self.angles_pred[frame], carID='pred')
        new_patches_true = self.draw_car(self.x[frame], self.y[frame], self.angles[frame], carID='true')

        for old_patch, new_patch in zip(self.pred_patches, new_patches_pred):
            if isinstance(old_patch, Circle) and isinstance(new_patch, Circle):
                old_patch.set_center(new_patch.center)
            elif isinstance(old_patch, Polygon) and isinstance(new_patch, Polygon):
                old_patch.set_xy(new_patch.get_xy())

        for old_patch, new_patch in zip(self.true_patches, new_patches_true):
            if isinstance(old_patch, Circle) and isinstance(new_patch, Circle):
                old_patch.set_center(new_patch.center)
            elif isinstance(old_patch, Polygon) and isinstance(new_patch, Polygon):
                old_patch.set_xy(new_patch.get_xy())
        speed_text = f'Car Speed: {self.vx[frame]:.2f} [m/s]'
        horizon_text = f'Horizon: {self.step_ls[frame+self.id]} [step]'

        # Remove existing steering wheel components
        patches_to_remove = [patch for patch in self.ax.patches if hasattr(patch, 'is_steering_wheel_component') and patch.is_steering_wheel_component]
        lines_to_remove = [line for line in self.ax.lines if hasattr(line, 'is_steering_wheel_component') and line.is_steering_wheel_component]

        for patch in patches_to_remove:
            patch.remove()
        for line in lines_to_remove:
            line.remove()

        # Zoomed in mode
        if self.focus_on:
            # Update the plot limits to focus on the arrow
            window_size = 50
            x_center = self.x[frame]
            y_center = self.y[frame]
            x_min = x_center - window_size * 3 / 4
            x_max = x_center + window_size * 3 / 4
            y_min = y_center - window_size * 1 / 3 
            y_max = y_center + window_size * 1 / 3
            self.ax.set_xlim(x_min, x_max)
            self.ax.set_ylim(y_min, y_max)
            self.ax.set_aspect('equal', 'box')
            # Set steering wheel icon
            self.steering_wheel_radius = window_size * 0.05
            steering_wheel_center_x = self.x[frame] + window_size/10
            steering_wheel_center_y = self.y[frame] + window_size/5
            text_x = steering_wheel_center_x 
            text_y = steering_wheel_center_y - window_size/10
            self.steering_angle_text.set_position((text_x, text_y))
        else:
            # Draw new steering wheel based on the current steering angle
            steering_wheel_center_x = self.x[frame] + self.steering_wheel_offset_x
            steering_wheel_center_y = self.y[frame] + self.steering_wheel_offset_y
            text_x = steering_wheel_center_x + self.text_offset_x
            text_y = steering_wheel_center_y + self.text_offset_y
            self.steering_angle_text.set_position((text_x, text_y))
        if pit_np.any():
            self.ax.legend(['left bound', 'right bound', 'pitlane', 'True position', 'Pred position', speed_text, horizon_text])
        else:
            self.ax.legend(['left bound', 'right bound', 'True position', 'Pred position', speed_text, horizon_text])
        current_steering_angle = self.steering_angle[frame] / 0.0666 / (math.pi/180.0) # Steering  wheel angle data
        self.draw_steering_wheel(self.ax, steering_wheel_center_x, steering_wheel_center_y, self.steering_wheel_radius, current_steering_angle)
        self.steering_angle_text.set_text(f'Steering wheel Angle: {current_steering_angle:.2f}°')
        # Calculate and display FPS
        current_time = time.time()
        if self.last_time is not None:
            fps = 1.0 / (current_time - self.last_time)
            # print(f"FPS: {fps}")
        self.last_time = current_time

    def animate(self):
        # self.fig, self.ax = plt.subplots(figsize=(12, 8))
        ani = FuncAnimation(self.fig, self.update, frames=range(len(self.x)), interval=50)
        plt.show()

if __name__ == "__main__":
    # Load racetrack boundaries
    racetracks   = 'tms'
    script_dir   = os.path.dirname(os.path.abspath(__file__))
    bounds_path    = os.path.join(script_dir, '..', 'data', 'racetracks', '{}/'.format(racetracks))
    left_bounds  = 'inner_bound.csv' #'left_line.csv'
    right_bounds = 'outer_bound.csv' #'right_line.csv'
    pitlane      = 'pit_through_local.csv'
    pit_np   = np.loadtxt(bounds_path+pitlane, delimiter=',', dtype=np.float64)
    left_np  = np.loadtxt(bounds_path+left_bounds, delimiter=',', dtype=np.float64)
    right_np = np.loadtxt(bounds_path+right_bounds, delimiter=',', dtype=np.float64)
    
    # Vehicle model & params
    params = vehicle.Vehicle_params()    
    vehicle = Vehicle(params)
    modelID = 0
    models = ['dynamic','kinematic']
    path   = os.path.join(script_dir, '..', 'data', 'example_runs/')
    dataset = 'autoverse_tms_purepursuit2'
    vehicle.gen_data(path, dataset, models)
    simulation = Simulation(vehicle, models[modelID])
    multi_step = True
    base_step  = 25 
    c1_scaler  = 9
    c2_scaler  = 2
    simulation.open_sim(multi_step, base_step, c1_scaler, c2_scaler, horizon_fixed=True)
    vx_start = 55 # start from the point vx >= [] m/s
    # Dynamic plot
    plotting = Plotting(simulation, vx_start, left_np, right_np, pit_np, Whole_track = False)
    plotting.animate()

