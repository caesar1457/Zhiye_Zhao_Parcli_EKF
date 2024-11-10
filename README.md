# Zhiye_Zhao_Parcli_EKF

This repository contains MATLAB code for implementing an Extended Kalman Filter (EKF) to estimate the state of a robot based on IMU and encoder data.

## Folder Structure and File Descriptions

### IMU_EKF Folder
Contains files for EKF implementation using IMU data:

- **EKF_parcli.m**: Main EKF script that estimates the robot's state using IMU data.
- **calculate_acceleration.m**: Helper function for calculating acceleration.
- **compute_Fk.m**: Function for computing the state transition matrix \( F_k \).
- **compute_Hq.m**: Function for computing the observation matrix \( H_q \).
- **expmap.m**: Exponential map function for quaternion transformations.
- **simulation_data.mat**: Data file providing necessary inputs for the EKF process, including state data and sensor readings.

### Kine_EKF_demo Folder
Contains files for a kinematic EKF implementation using encoder data only:

- **EKF_Kine_demo.m**: Example of a kinematic EKF implementation based solely on encoder data.
- **parody_mark1.urdf**: URDF file for the robot model, used for MATLAB simulation.
