% q = [1; 0; 0; 0];
% a = [0; 0; -9.81];
% H_q = compute_Hqd(q, a)

function H_q = compute_Hq(q, a)
    % compute_Hq calculates the derivative of the acceleration observation model with respect to the quaternion q
    % Inputs:
    %   q - 4x1 quaternion, [qw; qx; qy; qz]
    %   a - 3x1 acceleration vector (typically gravity or linear acceleration)
    % Output:
    %   H_q - 3x4 derivative matrix, the partial derivatives of acceleration with respect to quaternion

    qw = q(1);
    qx = q(2);
    qy = q(3);
    qz = q(4);
    
    % Derivatives of the components of the rotation matrix with respect to the quaternion
    H_q = zeros(3, 4);
    
    % Compute partial derivatives
    H_q(1, :) = 2 * [ qw*a(1) + qy*a(3) - qz*a(2), ...
                      qx*a(1) + qy*a(2) + qz*a(3), ...
                     -qy*a(1) + qx*a(2) - qw*a(3), ...
                      qz*a(1) + qw*a(2) - qx*a(3) ];
                  

    H_q(2, :) = 2 * [ qz*a(1) - qx*a(3) + qw*a(2), ...
                      qy*a(1) + qw*a(3) - qx*a(2), ...
                      qx*a(1) + qy*a(2) + qz*a(3), ...
                     -qw*a(1) + qz*a(2) - qy*a(3) ];
                  

    H_q(3, :) = 2 * [ -qy*a(1) + qx*a(2) + qw*a(3), ...
                      qz*a(1) - qw*a(2) + qx*a(3), ...
                      qy*a(1) + qx*a(2) + qz*a(3), ...
                      qx*a(1) + qw*a(2) - qz*a(3) ];
end
