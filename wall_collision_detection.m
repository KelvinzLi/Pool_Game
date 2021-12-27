% Kelvin Li
% detect collision between a moving ball and the walls

function [wall_id, dist] = wall_collision_detection(p, v, r, walls)
    % walls is the x/y coordinates of the 4 walls in anti-clockwise order
    % from the rightmost wall.

    % determine which walls the ball will possibly run into
    test_id = [1, 2, 3, 4];
    if v(1) >= 0 && v(2) >= 0
        test_id = [1, 2];
    elseif v(1) < 0 && v(2) >= 0
        test_id = [3, 2];
    elseif v(1) < 0 && v(2) < 0
        test_id = [3, 4];
    elseif v(1) >= 0 && v(2) < 0
        test_id = [1, 4];
    end

    % calculate required distance to travel before colliding with
    % respective walls
    x_dist = [0 0];
    count = 0;
    for id = test_id
        count = count + 1;
        if mod(id, 2) == 1
            if v ~= [0, 1]
                x_dist(count) = abs(walls(id) - p(1)) - r;
            else
                x_dist(count) = nan;
            end
        else
            x_dist(count) = abs((abs(walls(id) - p(2)) - r) * v(1) / v(2));
        end
    end

    % determine the wall that the ball will collide with and the distance
    % to travel before that happens
    [M, I] = min(x_dist);
    wall_id = test_id(I);
    if M > 0
        if v ~= [0, 1]
            dist = M / abs(v * [1, 0]');
        else
            dist = abs(walls(2) - p(2));
        end
    else
        dist = 0;
    end