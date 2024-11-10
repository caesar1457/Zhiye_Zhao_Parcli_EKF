# Zhiye_Zhao_Parcli_EKF

IMU_EKF Folder
•	EKF_parcli.m: Main EKF implementation file that estimates the robot's state using IMU data.
•	calculate_acceleration.m: Helper function for calculating acceleration.
•	compute_Fk.m: Function for computing the state transition matrix FkF_kFk.
•	compute_Hq.m: Function for computing the observation matrix HqH_qHq.
•	expmap.m: Exponential map function for quaternion transformations.
•	simulation_data.mat: Data file providing necessary inputs for the EKF process, including state data and sensor readings.

Kine_EKF_demo Folder
•	EKF_Kine_demo.m: Example of a kinematic EKF implementation based solely on encoder data.
•	parody_mark1.urdf: URDF file for the robot model, used for MATLAB simulation.
