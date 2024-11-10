% Exponential map function for quaternion update
function q_exp = expmap(theta)
    theta_norm = norm(theta);
    if theta_norm > 1e-6
        half_theta = theta_norm / 2;  % Use half-angle
        q_exp = [cos(half_theta); sin(half_theta) * theta / theta_norm];
    else
        q_exp = [1; 0; 0; 0];  % Return the unit quaternion
    end
end
