# Vehicle Model

Welcome to the "vehicle-model" repository. This repository includes implementations of single-track vehicle models, both dynamic and kinematic, tailored for simulation and analysis of vehicle dynamics.

## Repository Structure

The repository is organized into a primary `scripts` folder, which contains three key subdirectories:

- `data/`: Contains data collected from practice runs by a full-sized Indy racecar, as well as racetrack reference lines for simulation purposes.
- `tire_model/`: Includes the implementation of the Pacejka tire model (focused on lateral dynamics only). For a complete model, please refer to [Pacejka Tire Model](https://github.com/BrianN92/Pacejka-tire-model).
- `veh_models/`: Consists of a dynamic single-track model and a kinematic single-track model. Additionally, a validation script `dynamic_validation.py` is provided to visualize the correlation between model predictions and ground truth measurements.

## Usage

To use the models and scripts provided in this repository:

1. Navigate to the `scripts` directory.
2. Explore the `data` directory for real-world racecar data and track references.
3. For tire dynamics, look into the `tire_model` directory.
4. The `veh_models` directory contains the vehicle models and the validation script.

Run the validation script to compare the dynamic model predictions with actual data:

```bash
python dynamic_validation.py
