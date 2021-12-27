% Kelvin Li
% Detect possible collisions between balls

function [dist_array, time_array] = ball_collision_detection(ball_id, test_ids, p_array, v_array, s_array, a, r, walls, check_pos_range)
    ball_count = size(p_array);
    ball_count = ball_count(1);

    % initializing data recording arrays
    dist_array = nan(ball_count, 2);
    time_array = nan(ball_count, 1);

    for now_id = test_ids
        if now_id ~= ball_id
            % checker checks if the velocity and position in x-y directions
            % allow the two balls to possibly collide
            p1 = p_array(ball_id, :);
            p2 = p_array(now_id, :);
            v1 = v_array(ball_id, :);
            v2 = v_array(now_id, :);
            s1 = s_array(ball_id);
            s2 = s_array(now_id);

            % case 1: 2 moving balls collide with each other
            if s1 > 0 && s2 > 0
                % exclude certain situations where collision is 100%
                % impossible to avoid unnecessary calculation
                if ~(v1 * (p2 - p1)' < 0 && v2 * (p2 - p1)' > 0)
                    % coefficient for the quartic equation to calculate
                    % collision time
                    c4 = a ^ 2 * (v2 - v1) * (v2 - v1)' / 4;
                    c3 = -a * (v2 - v1) * (s2 * v2 - s1 * v1)';
                    c2 = (s2 * v2 - s1 * v1) * (s2 * v2 - s1 * v1)' - a * (v2 - v1) * (p2 - p1)';
                    c1 = 2 * (s2 * v2 - s1 * v1) * (p2 - p1)';
                    c0 = (p2 - p1) * (p2 - p1)' - 4 * r ^ 2;
    
                    % some indicators for whether real roots exist
                    [delta, P, D] = quartic_flags(c4, c3, c2, c1, c0);
    
                    if ~(delta > 0 && (P > 0 || D > 0))
                        solutions = roots([c4, c3, c2, c1, c0]);
    
                        t = nan;
                        % find the real root
                        for solution = solutions'
                            if solution == real(solution) && (isnan(t) || solution < t)
                                t = solution;
                            end
                        end
    
                        if t >= 0
                            d1 = s1 * t - a * t ^ 2 / 2;
                            d2 = s2 * t - a * t ^ 2 / 2;
    
                            collide_pos1 = p1 + d1 * v1;
                            collide_pos2 = p2 + d2 * v2;
   
                            % check if the collision positions are inside
                            % the walls if required by "check_in_range"
                            if (collide_pos1(1) < walls(1) - r && collide_pos1(1) > walls(3) + r && collide_pos1(2) < walls(2) - r && collide_pos1(2) > walls(4) + r ...
                                && collide_pos2(1) < walls(1) - r && collide_pos2(1) > walls(3) + r && collide_pos2(2) < walls(2) - r && collide_pos2(2) > walls(4) + r) || ~check_pos_range
    
                                dist_array(now_id, 1) = d1;
                                dist_array(now_id, 2) = d2;
    
                                time_array(now_id) = t;
                            end
                        end
                    end
                end

            % case 2: a moving ball collide with a still ball
            elseif s1 > 0 || s2 > 0
                % determine the moving ball
                if s1 > 0
                    moveID = 1;
                    stillID = 2;
                    dp = p2 - p1;
                    p = p1;
                    v = v1;
                    s = s1;
                elseif s2 > 0
                    moveID = 2;
                    stillID = 1;
                    dp = p1 - p2;
                    p = p2;
                    v = v2;
                    s = s2;
                end

                if v * dp' > 0
                    proj = dp * v';
                    normal_v = dp - proj * v;
                    normal = norm(normal_v);
    
                    % distance to move before collision
                    d = proj - sqrt((2 * r) ^ 2 - normal ^ 2);

                    % precaution for some accumulated computational errors
                    overlap_flag = false;
                    if (d < 0 && d >= -1e-10) % || (norm(dp) < 2 * r && proj > 0)
                        d = 0;
                        if p(1) < walls(1) - r || p(1) > walls(3) + r || p(2) < walls(2) - r || p(2) > walls(4) + r
                            overlap_flag = true;
                        end
                    end
    
                    collide_pos = p + v * d;
    
                    if (collide_pos(1) < walls(1) - r && collide_pos(1) > walls(3) + r && collide_pos(2) < walls(2) - r && collide_pos(2) > walls(4) + r) || overlap_flag || ~check_pos_range
                        if d >= 0
                            c = [-a/2, s, -d];
                            ts = roots(c);
            
                            t = nan;
    
                            for test_t = ts'
                                if (test_t >= 0 && test_t == real(test_t)) && (isnan(t) || test_t < t)
                                    t = test_t;
                                end
                            end
            
                            if ~isnan(t)
                                dist_array(now_id, moveID) = d;
                                dist_array(now_id, stillID) = 0;
                                time_array(now_id) = t;
                            end
                        end
                    end
                end
            end
        end
    end
end 

function [delta, P, D] = quartic_flags(c4, c3, c2, c1, c0)
    % indicators for the root distribution of quartic equations

    delta = 256 * c4^3 * c0^3 - 192 * c4^2 * c3 * c1 * c0^2 - 128 * c4^2 * c2^2 * c0^2 ...
            + 144 * c4^2 * c2 * c1^2 * c0 - 27 * c4^2 * c1^4 ...
            + 144 * c4 * c3^2 * c2 * c0^2 - 6 * c4 * c3^2 * c1^2 * c0 - 80 * c4 * c3 * c2^2 * c1 * c0 ...
            + 18 * c4 * c3 * c2 * c1^3 + 16 * c4 * c2^4 * c0 ...
            - 4 * c4 * c2^3 * c1^2 - 27 * c3^4 * c0^2 + 18 * c3^3 * c2 * c1 * c0 ...
            - 4 * c3^3 * c1^3 - 4 * c3^2 * c2^3 *c0 + c3^2 * c2^2 * c1^2;
    P = 8 * c4 * c2 - 3 * c3^2;
    D = 64 * c4^3 * c0 - 16 * c4^2 * c2^2 + 16 * c4 * c3^2 *c2 - 16 * c4^2 * c3 * c1 - 3 * c3^4;
end