clc;
clear all;
close all;

% Load simulation data
data = load('simulation_data.mat');
% data = readtable('20241107163958.csv');
state_time = data.state_time;
state_data = data.state_data;

% Define parameters
dt = 0.001;
total_time = 4;
num_steps = floor(total_time / dt);

% Generate time series
time_series = 0:dt:(total_time - dt);

% Extract true states from data
positions_data = state_data(1:num_steps, 5:7)';
% velocities_data = state_data(1:num_steps, 25:27)';
velocities_data = state_data(1:num_steps, 25:27)';
quaternions_data = state_data(1:num_steps, 1:4)';
% angular_velocity_true = state_data(1:num_steps, [23, 24, 25])';
angular_velocity_true = state_data(1:num_steps, 22:24)';

% % Plot angular velocity over time
% figure;
% subplot(3,1,1);
% plot(time_series, angular_velocity_true(1,:), 'r', 'LineWidth', 1.5);
% xlabel('Time (s)');
% ylabel('Angular Velocity X (rad/s)');
% title('Angular Velocity X over Time');
% 
% subplot(3,1,2);
% plot(time_series, angular_velocity_true(2,:), 'g', 'LineWidth', 1.5);
% xlabel('Time (s)');
% ylabel('Angular Velocity Y (rad/s)');
% title('Angular Velocity Y over Time');
% 
% subplot(3,1,3);
% plot(time_series, angular_velocity_true(3,:), 'b', 'LineWidth', 1.5);
% xlabel('Time (s)');
% ylabel('Angular Velocity Z (rad/s)');
% title('Angular Velocity Z over Time');
% 
% % Plot linear velocity over time
% figure;
% subplot(3,1,1);
% plot(time_series, velocities_data(1,:), 'r', 'LineWidth', 1.5);
% xlabel('Time (s)');
% ylabel('Velocity X (m/s)');
% title('Linear Velocity X over Time');
% 
% subplot(3,1,2);
% plot(time_series, velocities_data(2,:), 'g', 'LineWidth', 1.5);
% xlabel('Time (s)');
% ylabel('Velocity Y (m/s)');
% title('Linear Velocity Y over Time');
% 
% subplot(3,1,3);
% plot(time_series, velocities_data(3,:), 'b', 'LineWidth', 1.5);
% xlabel('Time (s)');
% ylabel('Velocity Z (m/s)');
% title('Linear Velocity Z over Time');

% Initialize matrices for storing true and estimated states
positions_true = zeros(3, num_steps);
velocities_true = zeros(3, num_steps);
quaternions_true = zeros(4, num_steps);
acceleration_true = zeros(3, num_steps);
g = [0; 0; -9.81];

% Calculate true acceleration
for i = 1:num_steps
    v = velocities_data(:, i);
    omega = angular_velocity_true(:, i);
    alpha = [0; 0; 0];
    r = positions_data(:, i);
    acceleration_true(:, i) = calculate_acceleration(v, omega, alpha, r, g);
end

% test1 = angular_velocity_true(:, 3500)

% Initial true and estimated states
r_k_true = positions_data(:,1);
v_k_true = velocities_data(:,1);
q_k_true = quaternions_data(:,1);
b_a_k_true = [0; 0; 0];
b_w_k_true = [0; 0; 0];
x_true = [r_k_true; v_k_true; q_k_true; b_a_k_true; b_w_k_true];

r_k_est = r_k_true;
v_k_est = v_k_true;
q_k_est = q_k_true;
b_a_k_est = [0; 0; 0];
b_w_k_est = [0; 0; 0];

% Initialize covariance and noise matrices
P_k = eye(16);
Q = diag([0.01*ones(1,3), 0.01*ones(1,3), 0.01*ones(1,4), 0.001*ones(1,3), 0.001*ones(1,3)]);
R = diag([0.001*ones(1,3), 0.001*ones(1,3)]);

positions_est = zeros(3, num_steps);
velocities_est = zeros(3, num_steps);
quaternions_est = zeros(4, num_steps);
noise_t = 0%0.0000001;
% Predefine angular velocity matrix
angular_velocity_est = zeros(3, num_steps - 1);
% dt = time_series(2) - time_series(1); % Assume time step is the same

for k = 1:num_steps -1
    % Get current and next quaternions
    q_current = quaternions_data(:, k);
    q_next = quaternions_data(:, k + 1);

    % Compute quaternion rate of change
    dq = quatmultiply(q_next', quatinv(q_current'))'; % Compute relative quaternion between q_next and q_current

    % Convert relative quaternion to angular velocity
    theta = 2 * acos(dq(1)); % Rotation angle
    if abs(theta) > 1e-9
        % Compute angular velocity vector
        omega = (dq(2:4) / sin(theta / 2)) * (theta / dt);
    else
        % If theta is very small, assume angular velocity is zero
        omega = [0; 0; 0];
    end

    % Store angular velocity
    angular_velocity_est(:, k) = omega;
end

% % Plot angular velocity over time
% figure;
% subplot(3,1,1);
% plot(time_series(1:end-1), angular_velocity_est(1,:), 'r', 'LineWidth', 1.5);
% xlabel('Time (s)');
% ylabel('Angular Velocity X (rad/s)');
% title('Computed Angular Velocity X over Time');
% 
% subplot(3,1,2);
% plot(time_series(1:end-1), angular_velocity_est(2,:), 'g', 'LineWidth', 1.5);
% xlabel('Time (s)');
% ylabel('Angular Velocity Y (rad/s)');
% title('Computed Angular Velocity Y over Time');
% 
% subplot(3,1,3);
% plot(time_series(1:end-1), angular_velocity_est(3,:), 'b', 'LineWidth', 1.5);
% xlabel('Time (s)');
% ylabel('Angular Velocity Z (rad/s)');
% title('Computed Angular Velocity Z over Time');

% Compute the error between estimated angular velocity and true angular velocity
angular_velocity_true_interp = angular_velocity_true(:, 1:num_steps - 1); % Align time steps
angular_velocity_error = sqrt(sum((angular_velocity_est - angular_velocity_true_interp).^2, 1));

% Plot the angular velocity error over time
figure;
plot(time_series(1:end-1), angular_velocity_error, '-m', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Angular Velocity Error (rad/s)');
title('Angular Velocity Error Between Computed and True Values Over Time');
legend('Angular Velocity Error');
grid on;

% Main EKF loop
for k = 1:num_steps
    a_k_true = acceleration_true(:, k); %1.81109769577700e-07/-2.03492491944509e-07/1.75001909996699e-06
    w_k_true = angular_velocity_true(:, k); %0.00155552864772623/-0.000173182852259141/-0.000181119649000092
    a_k_est = a_k_true + noise_t * randn(3,1);
    w_k_est = w_k_true + noise_t * randn(3,1);
    u_est = [a_k_est; w_k_est];

    % Update true state
    R_qk_true = R_from_quat(q_k_true);
    r_k1_true = r_k_true + dt * v_k_true + 0.5 * dt^2 * (R_qk_true.' * (a_k_true - b_a_k_true));
    v_k1_true = v_k_true + dt * (R_qk_true.' * (a_k_true - b_a_k_true));
    q_k1_true = quatnormalize(quatmultiply(q_k_true', expmap(dt * (w_k_true - b_w_k_true))'))';

    % R_qk_true = R_from_quat(q_k_true);
    % r_k1_true = positions_data(:,k);
    % v_k1_true = velocities_data(:,k);
    % q_k1_true = quaternions_data(:,k);
    x_true = [r_k1_true; v_k1_true; q_k1_true; b_a_k_true; b_w_k_true];

    % EKF Prediction Step
    R_qk_est = R_from_quat(q_k_est);
    r_k1_est = r_k_est + dt * v_k_est + 0.5 * dt^2 * (R_qk_est.' * (a_k_est - b_a_k_est));
    v_k1_est = v_k_est + dt * (R_qk_est.' * (a_k_est - b_a_k_est));
    q_k1_est = quatnormalize(quatmultiply(q_k_est', expmap(dt * (w_k_est - b_w_k_est))'))';
    x_est = [r_k1_est; v_k1_est; q_k1_est; b_a_k_est; b_w_k_est];

    % Compute Jacobian and measurement matrices
    F_k = compute_Fk(x_est, u_est, dt);
    H_k_acc = [zeros(3,6), compute_Hq(q_k_est, (a_k_est - b_a_k_est - g)), -R_qk_est', zeros(3,3)];
    H_k_gyro = [zeros(3,6), zeros(3,4), zeros(3,3), -eye(3)];
    H_k = [H_k_acc; H_k_gyro];

    % Update covariance
    P_k = F_k * P_k * F_k' + Q;

    % Measurements with noise
    noise_a = noise_t * randn(3,1);
    z_acc = R_qk_true' * (a_k_true - b_a_k_true) + noise_a;
    noise_w = noise_t * randn(3,1);
    z_gyro = w_k_true - b_w_k_true + noise_w;
    z_k = [z_acc; z_gyro];

    % Predicted measurements
    z_acc_est = R_qk_est' * (a_k_est - b_a_k_est);
    z_gyro_est = w_k_est - b_w_k_est;
    z_k_est = [z_acc_est; z_gyro_est];

    % EKF Update Step
    y_k = z_k - z_k_est;
    S_k = H_k * P_k * H_k' + R;
    K_k = P_k * H_k' / S_k;
    delta_x = K_k * y_k;

    % Update estimated state
    r_k1_est = r_k1_est + delta_x(1:3);
    v_k1_est = v_k1_est + delta_x(4:6);
    delta_q = expmap(delta_x(7:9));
    q_k1_est = quatnormalize(quatmultiply(q_k1_est', delta_q'))';
    b_a_k_est = b_a_k_est + delta_x(10:12);
    b_w_k_est = b_w_k_est + delta_x(13:15);

    % Update covariance with regularization
    epsilon = 1e-8;
    P_k = (eye(16) - K_k * H_k) * P_k * (eye(16) - K_k * H_k)' + K_k * R * K_k' + epsilon * eye(16);

    % Store results
    positions_est(:,k) = r_k1_est;
    velocities_est(:,k) = v_k1_est;
    quaternions_est(:,k) = q_k1_est;
    positions_true(:, k) = r_k1_true;
    velocities_true(:, k) = v_k1_true;
    quaternions_true(:, k) = q_k1_true;

    % Update for next step
    r_k_est = r_k1_est;
    v_k_est = v_k1_est;
    q_k_est = q_k1_est;
    r_k_true = r_k1_true;
    v_k_true = v_k1_true;
    q_k_true = q_k1_true;
end

% Plot Position, Velocity, and Quaternion Errors
position_error_total = sqrt(sum((positions_true - positions_est).^2, 1));
velocity_error_total = sqrt(sum((velocities_true - velocities_est).^2, 1));
figure; plot(time_series, position_error_total, 'DisplayName', 'Total Position Error'); xlabel('Time (s)'); ylabel('Position Error (m)'); legend;
figure; plot(time_series, velocity_error_total, 'DisplayName', 'Total Velocity Error'); xlabel('Time (s)'); ylabel('Velocity Error (m/s)'); legend;

% Quaternion Error
quat_error = zeros(1, num_steps);
for k = 1:num_steps
    dq = quatmultiply(quaternions_true(:, k)', quatinv(quaternions_est(:, k)'))';
    dq(1) = max(min(dq(1), 1), -1); % Clamp dq(1) within [-1, 1]
    quat_error(k) = 2 * acos(abs(dq(1))); % Rotation angle error
end
figure; plot(time_series, quat_error); xlabel('Time (s)'); ylabel('Quaternion Error (rad)');

% True vs. Estimated Positions
figure;
subplot(3,1,1); plot(time_series, positions_true(1,:), 'b', time_series, positions_est(1,:), 'r--'); xlabel('Time (s)'); ylabel('Position X (m)'); legend('True', 'Estimated');
subplot(3,1,2); plot(time_series, positions_true(2,:), 'b', time_series, positions_est(2,:), 'r--'); xlabel('Time (s)'); ylabel('Position Y (m)'); legend('True', 'Estimated');
subplot(3,1,3); plot(time_series, positions_true(3,:), 'b', time_series, positions_est(3,:), 'r--'); xlabel('Time (s)'); ylabel('Position Z (m)'); legend('True', 'Estimated');


% % Align positions_data to the timestamps of positions_true using linear interpolation
% positions_data_interp = interp1(time_series', positions_data', time_series', 'linear', 'extrap')';
% 
% % Calculate the error (Euclidean distance) between positions_true and positions_data for each timestamp
% position_error = sqrt(sum((positions_true - positions_data_interp).^2, 1));
% 
% % Plot the position error over time
% figure;
% plot(time_series, position_error, '-r', 'LineWidth', 1.5);
% xlabel('Time (s)');
% ylabel('Position Error (m)');
% title('Position Error Between positions true and positions data Over Time');
% legend('Position Error');
% grid on;
% 
% % Align velocities_data to the timestamps of velocities_true using linear interpolation
% velocities_data_interp = interp1(time_series', velocities_data', time_series', 'linear', 'extrap')';
% 
% % Calculate the error (Euclidean distance) between velocities_true and velocities_data for each timestamp
% velocities_error = sqrt(sum((velocities_true - velocities_data_interp).^2, 1));
% 
% % Plot the velocity error over time
% figure;
% plot(time_series, velocities_error, '-r', 'LineWidth', 1.5);
% xlabel('Time (s)');
% ylabel('Velocity Error (m/s)');
% title('Velocities Error Between velocities true and velocities data Over Time');
% legend('Velocity Error');
% grid on;
% 
% % Align quaternions_data to the timestamps of quaternions_true using linear interpolation
% quaternions_data_interp = interp1(time_series', quaternions_data', time_series', 'linear', 'extrap')';
% 
% % Ensure interpolated quaternions are normalized
% for k = 1:num_steps
%     quaternions_data_interp(:, k) = quatnormalize(quaternions_data_interp(:, k)')';
% end
% 
% % Calculate quaternion error (rotation angle difference) between quaternions_true and quaternions_data
% quat_error = zeros(1, num_steps);
% for k = 1:num_steps
%     % Compute the relative quaternion (difference) between true and interpolated data
%     dq = quatmultiply(quaternions_true(:, k)', quatinv(quaternions_data_interp(:, k)'))';
% 
%     % Ensure the first element of dq is positive to avoid 180-degree ambiguity
%     dq(1) = max(min(dq(1), 1), -1); % Clamp dq(1) within [-1, 1] to avoid complex values in acos
%     quat_error(k) = 2 * acos(abs(dq(1))); % Calculate rotation angle error in radians
% end
% 
% % Plot quaternion error over time
% figure;
% plot(time_series, quat_error, '-b', 'LineWidth', 1.5);
% xlabel('Time (s)');
% ylabel('Quaternion Error (rad)');
% title('Quaternion Error Between quaternions true and quaternions data Over Time');
% legend('Quaternion Error');
% grid on;
