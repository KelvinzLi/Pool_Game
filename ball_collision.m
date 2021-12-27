% Kelvin Li
% Calculation velocities of free-to-move balls after collision

function [vn1, sn1, vn2, sn2] = ball_collision(v1, s1, v2, s2, collide_n, restitute_coef, friction_coef)
    % transform into the coordinate system with x-axis aligned with the
    % normal direction of the 2 balls
    A = get_transform(collide_n);
    vs1 = A * (s1 * v1)';
    vs2 = A * (s2 * v2)';
    vn1 = vs1;
    vn2 = vs2;

    % new value for the velocity component normal to the two balls
    vn1(1) = (restitute_coef * (vs2(1) - vs1(1)) + vs1(1) + vs2(1)) / 2;
    vn2(1) = (restitute_coef * (vs1(1) - vs2(1)) + vs1(1) + vs2(1)) / 2;

    % new value for the velocity component tangent to both balls
    delta_v = abs(vn1(1) - vs1(1));
    sign = (vs1(2) - vs2(2)) / abs(vs1(2) - vs2(2));
    if friction_coef * delta_v < abs(vs1(2) - vs2(2)) / 2
        vn1(2) = vn1(2) - sign * friction_coef * delta_v;
        vn2(2) = vn2(2) + sign * friction_coef * delta_v;
    else
        vn1(2) = (vs1(2) + vs2(2)) / 2;
        vn2(2) = (vs1(2) + vs2(2)) / 2;
    end

    % transform back to original x-y coordinate
    vn1 = (A \ vn1)';
    vn2 = (A \ vn2)';
    sn1 = norm(vn1);
    sn2 = norm(vn2);
    vn1 = vn1 / sn1;
    vn2 = vn2 / sn2;
end