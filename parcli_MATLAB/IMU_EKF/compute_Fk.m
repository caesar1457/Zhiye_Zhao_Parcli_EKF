function F_k = compute_Fk(x, u, dt)
    % Compute the F_k matrix
    
    % Compute the F_c matrix
    F_c = compute_F_c(x, u, dt);

    % Define the identity matrix
    I = eye(16); % Changed to a 16x16 identity matrix
    % Calculate F_k = I + F_c * dt
    F_k = I + F_c * dt;
end

function Fc = compute_F_c(x, u, dt)
    % Compute the Jacobian matrix Fc = ∂f(x,u)/∂x
    % Inputs:
    %   x  - State vector [r_k; v_k; q_k; b_a_k; b_w_k]
    %   u  - Control input [a_k; omega_k]
    %   dt - Time step
    %   g  - Gravitational vector
    % Outputs:
    %   Fc - Jacobian matrix

    % Extract state variables
    r_k = x(1:3);
    v_k = x(4:6);
    q_k = x(7:10);   % Quaternion, assumed to be a unit quaternion [qw; qx; qy; qz]
    b_a_k = x(11:13);
    b_w_k = x(14:16);

    % Extract control inputs
    a_k = u(1:3);
    omega_k = u(4:6);

    % Compute rotation matrix
    R_q = quatToRotMat(q_k);

    % Compute the Jacobian matrix blocks
    % Initialize Jacobian matrix
    Fc = zeros(16,16);

    % Position partial derivatives with respect to state
    Fc(1:3,1:3) = eye(3);            % ∂r_{k+1}/∂r_k
    Fc(1:3,4:6) = dt * eye(3);       % ∂r_{k+1}/∂v_k

    % Velocity partial derivatives with respect to state
    Fc(4:6,4:6) = eye(3);            % ∂v_{k+1}/∂v_k
    % ∂v_{k+1}/∂q_k
    Fc(4:6,7:10) = compute_dvdq(q_k, a_k - b_a_k, dt);
    Fc(4:6,11:13) = -R_q' * dt;      % ∂v_{k+1}/∂b_{a,k}

    % Quaternion partial derivatives with respect to state
    Fc(7:10,7:10) = compute_dqdq(q_k, omega_k - b_w_k, dt);
    Fc(7:10,14:16) = compute_dqdbw(q_k, omega_k - b_w_k, dt);

    % Accelerometer bias partial derivatives
    Fc(11:13,11:13) = eye(3);        % ∂b_{a,k+1}/∂b_{a,k}

    % Gyroscope bias partial derivatives
    Fc(14:16,14:16) = eye(3);        % ∂b_{w,k+1}/∂b_{w,k}
end

function R = quatToRotMat(q)
    % Convert quaternion to rotation matrix
    % Input: q = [qw; qx; qy; qz]
    qw = q(1);
    qx = q(2);
    qy = q(3);
    qz = q(4);

    R = [1 - 2*(qy^2 + qz^2),     2*(qx*qy - qw*qz),     2*(qx*qz + qw*qy);
         2*(qx*qy + qw*qz),     1 - 2*(qx^2 + qz^2),     2*(qy*qz - qw*qx);
         2*(qx*qz - qw*qy),       2*(qy*qz + qw*qx),   1 - 2*(qx^2 + qy^2)];
end

function dvdq = compute_dvdq(q, acc, dt)
    % Compute ∂v_{k+1}/∂q_k
    % Inputs:
    %   q   - Current quaternion
    %   acc - (a_k - b_{a,k})
    %   dt  - Time step

    % Compute the derivative of the rotation matrix with respect to the quaternion
    % ∂R(q_k)^T/∂q_k
    dR_dq = compute_dRdq(q);

    % Compute ∂v_{k+1}/∂q_k = [∂(R(q_k)^T (a_k - b_{a,k}))/∂q_k] * dt
    dvdq = zeros(3,4);
    for i = 1:4
        dvdq(:,i) = dR_dq(:,:,i)' * acc * dt;
    end
end

function dR_dq = compute_dRdq(q)
    % Compute the derivative of the rotation matrix with respect to the quaternion
    % Output: dR_dq - A tensor of size 3x3x4 representing the derivative with respect to each quaternion component
    qw = q(1);
    qx = q(2);
    qy = q(3);
    qz = q(4);

    dR_dq = zeros(3,3,4);

    % Partial derivatives with respect to qw
    dR_dq(:,:,1) = 2 * [0,     -qz,      qy;
                        qz,       0,     -qx;
                       -qy,      qx,       0];

    % Partial derivatives with respect to qx
    dR_dq(:,:,2) = 2 * [0,       qy,      qz;
                        qy,   -2*qx,     -qw;
                        qz,      qw,   -2*qx];

    % Partial derivatives with respect to qy
    dR_dq(:,:,3) = 2 * [-2*qy,    qx,      qw;
                         qx,       0,      qz;
                        -qw,      qz,   -2*qy];

    % Partial derivatives with respect to qz
    dR_dq(:,:,4) = 2 * [-2*qz,   -qw,      qx;
                         qw,   -2*qz,      qy;
                         qx,      qy,       0];
end

function dqdq = compute_dqdq(q, omega, dt)
    % Use the quaternion update formula and the Jacobian matrix for quaternion multiplication
    omega_dt = omega * dt / 2;
    omega_norm = norm(omega_dt);
    epsilon = 1e-8; % Add a small offset to prevent division by zero

    if omega_norm < epsilon
        omega_dt_norm = omega_dt;
    else
        omega_dt_norm = omega_dt * (sin(omega_norm) / (omega_norm + epsilon));
    end

    delta_q = [cos(omega_norm); omega_dt_norm];

    dqdq = quatMultiplicationJacobian(q, delta_q);
    % display(dqdq)
end


function dqdbw = compute_dqdbw(q, omega, dt)
    % Compute ∂q_{k+1}/∂b_{w,k}
    % Inputs:
    %   q     - Current quaternion
    %   omega - (ω_k - b_{w,k})
    %   dt    - Time step

    omega_dt = (omega) * dt / 2;
    omega_norm = norm(omega_dt);
    epsilon = 1e-8; % Add a small offset to prevent division by zero

    if omega_norm < epsilon
        omega_dt_norm = omega_dt;
        cos_omega_norm = 1;
        sin_omega_norm = omega_norm;
    else
        omega_dt_norm = omega_dt / omega_norm;
        cos_omega_norm = cos(omega_norm);
        sin_omega_norm = sin(omega_norm);
    end

    % ∂delta_q/∂omega = Jacobian matrix
    ddeltaq_domega = zeros(4,3);
    for i = 1:3
        e_i = zeros(3,1);
        e_i(i) = 1;
        ddeltaq_domega(1,i) = - (omega_dt_norm' * e_i) * sin_omega_norm * dt / 2 / (omega_norm + epsilon);
        ddeltaq_domega(2:4,i) = (e_i * sin_omega_norm / (omega_norm + epsilon) + omega_dt_norm * (omega_dt_norm' * e_i) * (cos_omega_norm - sin_omega_norm / (omega_norm + epsilon))) * dt / 2;
    end

    % Compute ∂q_{k+1}/∂b_{w,k} = q_k ⊗ ∂delta_q/∂omega * (-1)
    dqdbw = -quatMultiplyJacobian_wrt_second(q, ddeltaq_domega);
    % disp(dqdbw)
end



function J = quatMultiplicationJacobian(q1, q2)
    % Compute the Jacobian matrix of quaternion multiplication q = q1 ⊗ q2 with respect to q1
    % Inputs:
    %   q1, q2 - Quaternions [qw; qx; qy; qz]
    % Output:
    %   J - 4x4 Jacobian matrix ∂(q1⊗q2)/∂q1

    qw1 = q1(1); qx1 = q1(2); qy1 = q1(3); qz1 = q1(4);
    qw2 = q2(1); qx2 = q2(2); qy2 = q2(3); qz2 = q2(4);

    J = [qw2, -qx2, -qy2, -qz2;
         qx2,  qw2, -qz2,  qy2;
         qy2,  qz2,  qw2, -qx2;
         qz2, -qy2,  qx2,  qw2];
end

function J = quatMultiplyJacobian_wrt_second(q1, dq2_domega)
    % Compute the Jacobian matrix of quaternion multiplication q = q1 ⊗ q2 with respect to q2,
    % then multiply by ∂q2/∂omega
    % Inputs:
    %   q1         - Quaternion [qw; qx; qy; qz]
    %   dq2_domega - 4x3 matrix, ∂q2/∂omega
    % Output:
    %   J - 4x3 matrix, ∂(q1⊗q2)/∂omega

    qw1 = q1(1); qx1 = q1(2); qy1 = q1(3); qz1 = q1(4);

    T = [ qw1, -qx1, -qy1, -qz1;
          qx1,  qw1,  qz1, -qy1;
          qy1, -qz1,  qw1,  qx1;
          qz1,  qy1, -qx1,  qw1];

    J = T * dq2_domega;
end
