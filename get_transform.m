% Kelvin Li
% generate transformation matrix from the line between centers of balls

function A = get_transform(u)
    T = [0, -1; 1, 0];
    A = eye(2) / [u', T * u'];