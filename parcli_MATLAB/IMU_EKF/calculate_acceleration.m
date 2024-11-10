function a = calculate_acceleration(v, omega, alpha, r, g)
    % calculate_acceleration calculates the acceleration given linear velocity and angular velocity, considering gravity.
    % Inputs:
    %   v      - Linear velocity vector (3x1)
    %   omega  - Angular velocity vector (3x1)
    %   alpha  - Angular acceleration vector (3x1)
    %   r      - Vector from the point to the rigid body's center of mass (3x1)
    %   g      - Gravitational acceleration vector (3x1), for example [0; 0; -9.81]
    % Output:
    %   a      - Total acceleration vector (3x1), including gravity

    % Calculate the centrifugal acceleration part
    centrifugal_acceleration = cross(omega, cross(omega, r));
    
    % Calculate the acceleration due to angular acceleration
    angular_acceleration_part = cross(alpha, r);
    
    % Add gravity and consider translational acceleration
    a = angular_acceleration_part + centrifugal_acceleration ;
end
