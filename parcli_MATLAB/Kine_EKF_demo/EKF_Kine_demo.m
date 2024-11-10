clc
clear
close
% Load the robot model
robot = importrobot('parody_mark1.urdf');
% show(robot, 'Visuals', 'on', 'Collisions', 'off');
% hold on;

% Get the number of joints in the robot
numJoints = numel(homeConfiguration(robot));

% Initialize the state vector [position; velocity; angle]
state = zeros(3 * numJoints, 1);  % Assume each joint has position, velocity, and angle states
P = eye(3 * numJoints);  % Covariance matrix
Q = 0.01 * eye(3 * numJoints);  % Process noise covariance
R = 0.1 * eye(3 * numJoints);   % Measurement noise covariance

% Set simulation parameters
dt = 0.1;  % Time step
T = 10;    % Simulation duration
N = T/dt;  % Number of time steps

% Create variables to store the states
position_true = zeros(numJoints, N);
velocity_true = zeros(numJoints, N);
angle_true = zeros(numJoints, N);
position_measure = zeros(numJoints, N);
velocity_measure = zeros(numJoints, N);
angle_measure = zeros(numJoints, N);
position_est = zeros(numJoints, N);
velocity_est = zeros(numJoints, N);
angle_est = zeros(numJoints, N);

% Get the initial configuration
config = homeConfiguration(robot);

% Main loop: Simulate true states, perform EKF estimation, and visualize motion
figure;
for k = 2:N
    % Simulate true joint states (random motion)
    for j = 1:numJoints
        velocity_true(j, k) = velocity_true(j, k-1) + 0.1 * randn;
        position_true(j, k) = position_true(j, k-1) + velocity_true(j, k) * dt;
        angle_true(j, k) = angle_true(j, k-1) + 0.05 * randn;
    end
    
    % Measurements (with added noise)
    velocity_measure(:, k) = velocity_true(:, k) + 0.2 * randn(numJoints, 1);
    position_measure(:, k) = position_true(:, k) + 0.2 * randn(numJoints, 1);
    angle_measure(:, k) = angle_true(:, k) + 0.1 * randn(numJoints, 1);
    
    % EKF prediction step
    state_prior = state;  % Assume a linear motion model
    F = eye(3 * numJoints);  % State transition matrix
    P = F * P * F' + Q;  % Predict covariance
    
    % EKF update step
    H = eye(3 * numJoints);  % Observation matrix
    z = [velocity_measure(:, k); position_measure(:, k); angle_measure(:, k)];
    y = z - H * state_prior;  % Calculate measurement residual
    S = H * P * H' + R;  % Calculate residual covariance
    K = P * H' / S;  % Calculate Kalman gain
    
    state_posterior = state_prior + K * y;  % Update state estimate
    P = (eye(3 * numJoints) - K * H) * P;  % Update covariance matrix
    
    % Store posterior estimates
    velocity_est(:, k) = state_posterior(1:numJoints);
    position_est(:, k) = state_posterior(numJoints+1:2*numJoints);
    angle_est(:, k) = state_posterior(2*numJoints+1:end);
    
    % Update the state
    state = state_posterior;
    
    % Visualize the current state of the robot
    for j = 1:numJoints
        config(j).JointPosition = angle_est(j, k);  % Update configuration with estimated joint angles
    end
    show(robot, config, 'Visuals', 'on', 'Collisions', 'off');
    drawnow;  % Update the figure in real-time
end

for j = 1:numJoints
    figure;

    % Plot position estimation
    subplot(3, 1, 1);  % Create 3 subplots in one figure
    hold on;
    plot(position_true(j, :), '-', 'LineWidth', 1);
    plot(position_measure(j, :), '--', 'LineWidth', 1);
    plot(position_est(j, :), ':', 'LineWidth', 1);
    legend('True', 'Measured', 'Estimated');
    xlabel('Time Step');
    ylabel(['Position of Joint ', num2str(j)]);
    title(['Position Estimation for Joint ', num2str(j)]);

    % Plot velocity estimation
    subplot(3, 1, 2);
    hold on;
    plot(velocity_true(j, :), '-', 'LineWidth', 1);
    plot(velocity_measure(j, :), '--', 'LineWidth', 1);
    plot(velocity_est(j, :), ':', 'LineWidth', 1);
    legend('True', 'Measured', 'Estimated');
    xlabel('Time Step');
    ylabel(['Velocity of Joint ', num2str(j)]);
    title(['Velocity Estimation for Joint ', num2str(j)]);

    % Plot angle estimation
    subplot(3, 1, 3);
    hold on;
    plot(angle_true(j, :), '-', 'LineWidth', 1);
    plot(angle_measure(j, :), '--', 'LineWidth', 1);
    plot(angle_est(j, :), ':', 'LineWidth', 1);
    legend('True', 'Measured', 'Estimated');
    xlabel('Time Step');
    ylabel(['Angle of Joint ', num2str(j)]);
    title(['Angle Estimation for Joint ', num2str(j)]);
end
