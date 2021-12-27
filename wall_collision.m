% Kelvin Li
% Calculate the velocity of the ball after colliding with a wall

function [v_new, s_new] = wall_collision(v, s, wall_n, restitute_coef, friction_coef)
    % transform to the coordinate system with x axis aligned to the normal
    % direction of the wall
    A = get_transform(wall_n);
    vn = A * (s * v)';

    % new value for the velocity component normal to the wall
    vn(1) = restitute_coef * vn(1);

    % new value for the velocity component tangent to the wall
    delta_v = abs((1 - restitute_coef) * vn(1));
    sign = vn(2) / abs(vn(2));
    if friction_coef * delta_v < abs(vn(2))
        vn(2) = vn(2) - sign * friction_coef * delta_v;
    else
        vn(2) = 0;
    end

    % transform to original x-y coordinate
    Tv = [-1, 0; 0, 1];
    v_new = (A \ (Tv * vn))';
    s_new = norm(v_new);
    v_new = v_new / s_new;
    






