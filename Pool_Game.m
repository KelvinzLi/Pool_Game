% Kelvin Li
% a game of pool (as the file name suggested)

function Pool_Game()
    % get the current position of the mouse on the figure
    function pos = get_mouse_pos()
        pos = get(gca, 'CurrentPoint');
        pos = pos(1, 1:2);
    end

    function OnMouseDown(~, ~)
        % if no ball is in motion
        if sum(MotionTracker) == 0
            % if player isn't allowed to freely place ball
            if ~toPlaceBall    
                pos = get_mouse_pos();
                if (norm(pos - pos_array(1, :)) <= ball_r)
                    isDragging = true;
                end

            % if allowed
            else
                pos = get_mouse_pos();
                flag = place_cue_ball(pos);
    
                if flag
                    toPlaceBall = false;
                end
            end
        end
    end

    function OnMouseUp(~, ~)
        if isDragging
            % delete auxiliary lines and initiate collision detection
            delete(past_line)
            delete(auxiliary_plot_array)
            isDragging = false;
            MotionTracker(1) = 1;

            justStopped = true;
            stroke_ball_id = nan;
            potted_ball_ids = [];

            CollisionManager([1]);
        end
    end

    function OnMouseMove (~, ~)
        if isDragging
            mouse_pos = get_mouse_pos();

            % translate current draging position into initial velocity of
            % the cue ball, and updating the arrays
            v = (mouse_pos - pos_array(1, :));
            v_unit = -1 * v / norm(v);
            v_array(1, :) = v_unit;

            coef = min(max_drag, norm(v)) / max_drag;

            speed_array(1) = max_speed * coef;
            color = (1 - coef) * (([255, 204, 207] / 255) - [1 0 0]) + [1 0 0];

            % calculate details of auxiliary line and drawing it
            delete(past_line)
            drag_pos = pos_array(1, :) + (v / norm(v)) * min(max_drag, norm(v));
            past_line = plot([pos_array(1, 1) drag_pos(1)], [pos_array(1, 2) drag_pos(2)], "Color", color, 'LineWidth',2);
            hold on

            if auxiliary_flag
                delete(auxiliary_plot_array)
                rotA = [0, -1; 1, 0];
                auxiliary_plot_array = gobjects(2, 1);
                count = 0;
                for ii = [1, -1]
                    count = count + 1;
                    p_side = pos_array(1, :) + ii * ball_r * (rotA * v_unit')';
                    [~, d] = wall_collision_detection(p_side, v_unit, 0, walls);
                    p_end = p_side + d * v_unit;
                    auxiliary_plot_array(count) = plot([p_side(1), p_end(1)], [p_side(2), p_end(2)], "Color", [0, 0, 0]);
                    hold on
                end
            end
        end
    end

    function OnKeyDown(~, ko)
        % click to restart
        if ko.Key == "r"
            player_colors = zeros(2, 1);

            now_player = 1;
            winner = nan;
            score = zeros(2, 1);
            title([])
            xlabel([])

            init_setup;
            xlabel(player_names(now_player) + ": Drag and release to shoot the cue ball!")
        end

        % click to tur on or off the auxiliary line
        if ko.Key == "h"
            if auxiliary_flag
                delete(auxiliary_plot_array)
            end
            auxiliary_flag = ~auxiliary_flag;
        end
    end

    % resume intial game setup
    function init_setup()
        % Setting up balls
        pos_array = zeros(ball_count, 2);
        v_array = zeros(ball_count, 2);
        speed_array = zeros(ball_count, 1);
        v_next_array = zeros(ball_count, 2);
        speed_next_array = zeros(ball_count, 1);
        
        dist_array = zeros(ball_count, 1);
        time_array = nan(ball_count, 1);

        collide_ball_array = zeros(ball_count, 1);

        pos_array(1, :) = ball_cue_pos;

        on_table_id = 1: ball_count;
        enter_hole_tracker = zeros(ball_count, 1);

        MotionTracker = zeros(ball_count, 1);

        % calculating intial ball positions
        id = 1;
        for ii = 0: ball_rows
            x = - sqrt(3) * (2 * ii * ball_r + ii * ball_margin) / 2;
            y_base = - (2 * ii * ball_r + ii * ball_margin) / 2;
            for jj = 0: ii
                id = id + 1;
                y = y_base + jj * (2 * ball_r + ball_margin);
                pos_array(id, :) = [x, y] + ball_tri_pos;
            end
        end

        for id = on_table_id
            add_ball_plot(id);
        end

        % Plotting takes much more time than calculating the path. Plotting
        % using "rectangle" would result in bad animation quality
        % Through experiment, I found that this is the quickest way to draw "circles"
        delete(ball_patch)
        ball_patch = patch(ball_plot_x, ball_plot_y, color_array, 'EdgeColor', "none");
    end

    % delete all information of the given ball (used when the ball is potted)
    function delete_ball(id)
        on_table_id = on_table_id(on_table_id ~= id);

        pos_array(id, :) = [nan, nan];
        v_array(id, :) = [nan, nan];
        speed_array(id) = nan;
        v_next_array(id) = nan;
        speed_next_array(id) = nan;
        
        dist_array(id) = nan;
        time_array(id) = nan;

        MotionTracker(id) = 0;
    end

    % detect if the player intended position to place the ball overlap with
    % other balls or the wall; otherwise, place the ball there
    function flag = place_cue_ball(pos)
        flag = in_box_checker(pos, kitchen + [-1, -1, 1, 1] * (ball_r + 1e-10));
        if flag
            for id = on_table_id
                if norm(pos_array(id, :) - pos) < 2 * ball_r
                    flag = false;
                    break
                end
            end
        end

        if flag
            pos_array(1, :) = pos;
            on_table_id = [1, on_table_id];

            add_ball_plot(1)

            delete(ball_patch)
            ball_patch = patch(ball_plot_x(:, on_table_id), ball_plot_y(:, on_table_id), color_array(on_table_id, :, :), 'EdgeColor', "none");
        end
    end

    function OnFigureClose(~, ~)
        isRunning = false;
        delete(f)
    end

    % plot the pool table
    function plot_table_w_hole()
        angles = [pi/2, pi/4, -3*pi/4, pi/4, 5*pi/4, pi; ...
                    pi, 3*pi/4, -pi/4, 3*pi/4, 7*pi/4, 3*pi/2; ...
                    3*pi/2, pi, 0, pi, 0, -pi/2; ...
                    3*pi/2, 5*pi/4, pi/4, 5*pi/4, pi/4, 0; ...
                    0, -pi/4, 3*pi/4, 7*pi/4, 3*pi/4, pi/2; ...
                    pi/2, 0, -pi, 0, pi, pi/2];
        table_x = [];
        table_y = [];

        for ii = 1: 6
            for jj = 1: 3
                if jj == 2
                    now_r = hole_r;
                    angle_step = pi / 16;
                else
                    angle_step = - pi / 16;
                    if mod(ii, 3) ~= 0
                        now_r = bound_r;
                    else
                        now_r = hole_r;
                    end
                end
                th = angles(ii, 2 * jj - 1): angle_step : angles(ii, 2 * jj);
                table_x = [table_x, (now_r * cos(th) + bound_center_pos(ii, 2*jj-1))];
                table_y = [table_y, (now_r * sin(th) + bound_center_pos(ii, 2*jj))];
            end
        end

        rectangle('Position',[-margin, -margin, (y_length + margin) * 2, y_length + margin * 2], 'FaceColor',[0 0 0], 'Curvature', 0.2);
        hold on
        pgon = polyshape(table_x, table_y, "Simplify", false);
        pgon_plot = plot(pgon);
        pgon_plot.FaceAlpha = 1;
        pgon_plot.FaceColor = [0 0.75 0.75];

        drawnow
        hold on

        plot([1.5 * y_length, 1.5 * y_length], [0 + 1e-5, y_length - 1e-5], "w", 'LineWidth', 0.001)
        hold on
    end

    % update the position of a ball on the plot
    function add_ball_plot(id)
        th = linspace(0, 2*pi * (circle_sides - 1) / circle_sides, circle_sides);
        ball_plot_x(:, id) = (ball_r * cos(th) + pos_array(id, 1))';
        ball_plot_y(:, id) = (ball_r * sin(th) + pos_array(id, 2))';
    end

    % check if the ball is still inside the boundary of the pool table
    % including the slots; if the ball is in the slots, return which slot
    % it is in
    function [flag, slot_id] = in_boundary_checker(pos)
        flag = false;
        slot_id = 0;
        % if the ball is in side the table but not in the slot
        if in_box_checker(pos, walls + [-1, -1, 1, 1] * (ball_r - 1e-8))
            flag = true;
        else
            % iterate over slots to check which one it's inside
            for ii = 1: 6
                if in_box_checker(pos, hole_region(ii, :) + [-1, -1, 1, 1] * (-1e-10))
                    if mod(ii, 3) ~= 0
                        now_bound_r = bound_r;
                    else
                        now_bound_r = hole_r;
                    end

                    if (min(norm(pos - bound_center_pos(ii, 1: 2)), norm(pos - bound_center_pos(ii, 5: 6))) > now_bound_r + ball_r - 1e-10)
                        flag = true;
                        slot_id = ii;
                    else
                        flag = false;
                        break
                    end
                end
            end
        end
    end

    % detect which part of the curves of the slot will the ball
    % collide with
    function [collide_d, normal, bound_ball_id] = slot_collision_detection(p, v, s, slot_id)
        bound_d_array = nan(3, 1);
        % iterate over curves to check how far the ball needs to travel
        % before colliding with them
        for ii = 1: 3
            dp = bound_center_pos(slot_id, 2 * ii - 1: 2 * ii) - p;
            if ii == 2
                r_sum = hole_r; % = (hole_r - ball_r) + ball_r
            else
                if mod(slot_id, 3) ~= 0
                    r_sum = bound_r + ball_r;
                else
                    r_sum = hole_r + ball_r;
                end
            end
            bound_d_array(ii) = still_ball_collision_distance(dp, v, r_sum);
        end

        % found out which curve it will collide with and return the normal
        % direction of the collision
        [collide_d, bound_ball_id] = min(bound_d_array);
        if ~isnan(collide_d) && in_box_checker(p + collide_d * v, hole_region(slot_id, :) + [-1, -1, 1, 1] * (-1e-10))
            if collide_d > s ^ 2 / (2 * acc)
                collide_d = s ^ 2 / (2 * acc);
                normal = [0, 0];
            else
                normal = bound_center_pos(slot_id, 2 * bound_ball_id - 1: 2 * bound_ball_id) - (p + collide_d * v);
                normal = normal / norm(normal);
            end
        else
            collide_d = nan;
            normal = nan(2, 1);
        end
    end

    % calculate how far the ball will travel before colliding with the
    % boundary of the pool table or stop
    function [collide_d, normal, enter_hole_flag] = boundary_collision_detection(p, v, s)
        flag = false;
        enter_hole_flag = false;

        collide_d = nan;
        normal = nan(1, 2);

        % if the ball is inside a slot
        if ~in_box_checker(p, walls + [-1, -1, 1, 1] * (ball_r - 1e-5))
            [~, slot_id] = in_boundary_checker(p);
            if slot_id ~= 0
                [collide_d, normal, bound_ball_id] = slot_collision_detection(p, v, s, slot_id);
                if ~isnan(collide_d)
                    % flag is turned on if the ball will directly collide
                    % with the slot it's in
                    flag = true;
    
                    if bound_ball_id == 2
                        enter_hole_flag = true;
                    end
                end
            end
        end

        if ~flag
            [wall_id, d] = wall_collision_detection(p, v, ball_r, walls);

            % if the ball won't stop before colliding
            if d < s ^ 2 / (2 * acc)
                collide_p = p + d * v;
                slot_id = nan;
                wall_type = mod(wall_id, 2) + 1;

                % check if the ball will enter a slot
                for ii = 1: 3
                    if collide_p(wall_type) >= slot_pos_on_wall(wall_type, 2 * ii - 1) && collide_p(wall_type) <= slot_pos_on_wall(wall_type, 2 * ii)
                        slot_id = slot_pos_to_index(wall_id, ii);
                        break
                    end
                end

                if ~isnan(slot_id)
                    [collide_d, normal, bound_ball_id] = slot_collision_detection(p, v, s, slot_id);
                    if bound_ball_id == 2
                        enter_hole_flag = true;
                    end
                else
                    collide_d = d;
                    if mod(wall_id, 2) == 1
                        normal = [1 0];
                    else
                        normal = [0 1];
                    end
                end
            else
                collide_d = s ^ 2 / (2 * acc);
                normal = [0, 0];
            end
        end
    end

    % find out which wall/ball the given ball will collide with first
    function [flag, collision_id] = BallCollisionManager(update_id, test_ids, bound_d_array)
        % calculate all time and distances needed to travel to collide with other balls
        [now_dists, now_times] = ball_collision_detection(update_id, test_ids, pos_array, v_array, speed_array, acc, ball_r, walls, false);
        [sort_time, sort_id] = sort(now_times);

        flag = false;
        collision_id = 0;
        if ~isnan(sort_time(1))
            for count = 1: ball_count
                ball_id = sort_id(count);
                % if the given ball will collide with the current examining
                % ball quicker than all examined balls before and before the balls collide with the boundary
                if ~isnan(now_times(ball_id)) && (now_times(ball_id) < time_array(ball_id) || isnan(time_array(ball_id)))
                    if (now_dists(ball_id, 1) <= bound_d_array(update_id) || isnan(bound_d_array(update_id))) ...
                            && (now_dists(ball_id, 2) <= bound_d_array(ball_id) || isnan(bound_d_array(ball_id)))

                        % update information
                        dist_array(update_id) = now_dists(ball_id, 1);
                        dist_array(ball_id) = now_dists(ball_id, 2);
    
                        time_array(update_id) = now_times(ball_id);
                        time_array(ball_id) = now_times(ball_id);
    
                        normal = (pos_array(update_id, :) + dist_array(update_id) * v_array(update_id, :)) ...
                                    - (pos_array(ball_id, :) + dist_array(ball_id) * v_array(ball_id, :));
                        normal = normal / norm(normal);
    
                        [v1, s1, v2, s2] = ball_collision(v_array(update_id, :), speed_array(update_id, :) - acc * time_array(update_id), ...
                                                            v_array(ball_id, :), speed_array(ball_id, :) - acc * time_array(ball_id), normal, ...
                                                            restitute_coef, friction_coef);
    
                        v_next_array(update_id, :) = v1;
                        v_next_array(ball_id, :) = v2;
    
                        speed_next_array(update_id) = s1;
                        speed_next_array(ball_id) = s2;
    
                        if s1 > 0
                            MotionTracker(update_id) = 1;
                        end
                        if s2 > 0
                            MotionTracker(ball_id) = 1;
                        end

                        if enter_hole_tracker(update_id) == 1
                            enter_hole_tracker(update_id) = 0;
                        end
                        if enter_hole_tracker(ball_id) == 1
                            enter_hole_tracker(ball_id) = 0;
                        end
    
                        flag = true;
                        collision_id = ball_id;
                        break
                    end
                end
            end
        end
    end

    function CollisionManager(update_array)
        next_update_array = [];
        test_ids = on_table_id;
        bound_d_array = nan(ball_count, 1);
        normal_array = nan(ball_count, 2);
        hole_flag_array = boolean(zeros(ball_count, 1));
        % calculate the details of all balls colliding with the boundary
        for id = on_table_id
            if speed_array(id) > 0
                [bound_d_array(id), normal_array(id, :), hole_flag_array(id)] = boundary_collision_detection(pos_array(id, :), v_array(id, :), speed_array(id));
            end
        end

        for update_id = update_array
            test_ids = test_ids(test_ids ~= update_id);

            [flag, collision_id] = BallCollisionManager(update_id, test_ids, bound_d_array);
            if flag
                if collision_ball_array(collision_id) ~= 0
                    next_update_id = collision_ball_array(collision_id);
                    if sum(update_id == next_update_id) == 0
                        time_array(next_update_id) = nan;
                        v_next_array(next_update_id, :) = [0, 0];
                        speed_next_array(next_update_id) = 0;
    
                        next_update_array = [next_update_array, next_update_id];
                    end
                end

                collision_ball_array(update_id) = collision_id;
                collision_ball_array(collision_id) = update_id;
            else
                collision_ball_array(update_id) = 0;
                if speed_array(update_id) > 0
                    bound_d = bound_d_array(update_id);
                    normal = normal_array(update_id, :);
                    enter_hole_flag = hole_flag_array(update_id);

                    % calculate velocity information after collision
                    if normal == [0, 0]
                        t = speed_array(update_id) / acc;
                        dist_array(update_id) = (speed_array(update_id) ^ 2) / (2 * acc);
                        time_array(update_id) = t;
    
                        v_next_array(update_id, :) = [0, 0];
                        speed_next_array(update_id) = 0;
                    else
                        c = [-acc/2, speed_array(update_id), -bound_d];
                        
                        ts = roots(c);
                    
                        t = nan;
                        for test_t = ts'
                            if (test_t >= 0 && test_t == real(test_t)) && (isnan(t) || test_t < t)
                                t = test_t;
                            end
                        end

                        dist_array(update_id) = bound_d;
                        time_array(update_id) = t;

                        collision_speed = speed_array(update_id) - acc * t;
                        [v_next_array(update_id, :), speed_next_array(update_id)] = wall_collision(v_array(update_id, :), collision_speed, normal, restitute_coef, friction_coef);

                        if enter_hole_flag
                            enter_hole_tracker(update_id) = 1;
                        end
                    end
                end
            end
        end

        % detect the first ball that the cue ball collide with
        if isnan(stroke_ball_id)
            stroke_ball_ids = find(time_array >= 0);
            if length(stroke_ball_ids) > 1
                stroke_ball_id = stroke_ball_ids(2);
            end
        end

        if ~isempty(next_update_array)
            CollisionManager(next_update_array)
        end
    end

    % update ball information in the game loop
    function update_balls(t)
        update_record = zeros(ball_count, 1);
        for id = on_table_id
            if time_array(id) >= 0
                if time_array(id) > t
                    if speed_array(id) ~= 0 || ((speed_next_array(id) == 0) && (speed_array(id) <= acc * t + 1e-20))
                        % if the ball is still moving normally
                        if speed_array(id) > acc * t + 1e-20
                            now_d = calc_d(speed_array(id), acc, t);
                            speed_array(id) = speed_array(id) - acc * t;
                            time_array(id) = time_array(id) - t;
                        % if the ball is coming to a rest
                        else
                            now_d = speed_array(id) ^ 2 / (2 * acc);
                            time_array(id) = nan;
                            speed_array(id) = 0;
                            MotionTracker(id) = 0;
    
                            update_record(id) = 1;
                        end
    
                        pos_array(id, :) = pos_array(id, :) + v_array(id, :) * now_d;
                    % if the ball is at rest now but will collide with some
                    % ball later
                    else
                        time_array(id) = time_array(id) - t;
                    end
                % if the ball collide with something inside the time step
                else
                    % if it won't enter the hole
                    if enter_hole_tracker(id) == 0
                        now_d_before = calc_d(speed_array(id), acc, time_array(id));
                        now_d_after = calc_d(speed_next_array(id), acc, t - time_array(id));
        
                        pos_array(id, :) = pos_array(id, :) + v_array(id, :) * now_d_before + v_next_array(id, :) * now_d_after;
    
                        v_array(id, :) = v_next_array(id, :);
                        speed_array(id) = speed_next_array(id) - acc * (t - time_array(id));
    
                        update_record(id) = 1;
                    % if it enters the hoel
                    else
                        delete_ball(id)

                        potted_ball_ids = [potted_ball_ids, id];
                    end
                end

                % In some weird situations (probably due to the precision
                % of computer calculation), the program will not realize
                % that the ball will stop. We need the following "if" to
                % double check
                if speed_array(id) == 0 && speed_next_array(id) == 0
                    time_array(id) = nan;
                    MotionTracker(id) = 0;
                    
                    update_record(id) = 1;
                end
            % precaution
            elseif time_array(id) < 0
                time_array(id) = nan;
                speed_array(id) = 0;
                MotionTracker(id) = 0;
        
                update_record(id) = 1;
            end
        end

        % if any ball experience collision or come to rest,
        % we will need to update this information to all other balls
        if sum(update_record) > 0
            update_array = 1: ball_count;
            update_array = update_array(update_record == 1);
            CollisionManager(update_array);
        end
    end

    % Settings
    y_length = 10;
    margin = 1.8;
    circle_sides = 16;

    ball_r = 0.32;
    hole_r = 1.6 * ball_r;
    ball_margin = 0.1 * ball_r; % distance between boundary of balls in the intial triangle

    ball_rows = 5; % rows of balls in the intial triangle
    ball_count = 1 + sum(1: ball_rows);

    ball_tri_pos = [0.5 * y_length + sqrt(3) * (2 * ball_r + ball_margin), 0.5 * y_length]; % initial position of the tip of the triangle
    ball_cue_pos = [1.5 * y_length, 0.5 * y_length];

    acc = 12; % acceleration
    time_step = 0.025; % update time-step
    max_speed = 48;
    max_drag = y_length / 2;

    % physics coefficient used when colliding
    restitute_coef = 0.8;
    friction_coef = 0.1;

    auxiliary_flag = true;

    % information used in calculation
    bound_r = (3 + 2 * sqrt(2)) * hole_r;
    bound_side = bound_r - 2 * hole_r;
    walls = [2 * y_length, y_length, 0, 0];
    kitchen = [2 * y_length, y_length, 1.5 * y_length, 0];
    hole_anchor = [2 * y_length, 0; ...
                    2 * y_length, y_length; ...
                    y_length, y_length; ...
                    0, y_length; ...
                    0, 0; ...
                    y_length, 0];
    hole_region = [2 * hole_r, bound_side, -bound_side, -2 * hole_r; ...
                    2 * hole_r, 2 * hole_r, -bound_side, -bound_side; ...
                    2 * hole_r, 2 * hole_r,  -2 * hole_r, -ball_r
                    bound_side, 2 * hole_r, -2 * hole_r, -bound_side; ...
                    bound_side, bound_side, -2 * hole_r, -2 * hole_r; ...
                    2 * hole_r, ball_r, -2 * hole_r, -2 * hole_r] ...
                    + [hole_anchor, hole_anchor];
    bound_center_pos = [-bound_side, -bound_r, hole_r, -hole_r, bound_r, bound_side; ...
                    bound_r, -bound_side, hole_r, hole_r, -bound_side, bound_r; ...
                    2 * hole_r, hole_r, 0, hole_r, -2 * hole_r, hole_r; ...
                    bound_side, bound_r, -hole_r, hole_r, -bound_r, -bound_side; ...
                    -bound_r, bound_side, -hole_r, -hole_r, bound_side, -bound_r; ...
                    -2 * hole_r, -hole_r, 0, -hole_r, 2 * hole_r, -hole_r] ...
                    + [hole_anchor, hole_anchor, hole_anchor];

    slot_pos_on_wall = [0, bound_side, y_length - 2 * hole_r, y_length + 2 * hole_r, 2 * y_length - bound_side, 2 * y_length; ... % for wall 2, 4
                        0, bound_side, nan, nan, y_length - bound_side, y_length]; % for wall 1, 3 (use mod(wall_id, 2) + 1)
    slot_pos_to_index = [1, nan, 2;
                        4, 3, 2;
                        5, nan, 4;
                        5, 6, 1];

    % Setting up balls
    pos_array = zeros(ball_count, 2);
    v_array = zeros(ball_count, 2);
    speed_array = zeros(ball_count, 1);
    v_next_array = zeros(ball_count, 2);
    speed_next_array = zeros(ball_count, 1);
    
    dist_array = zeros(ball_count, 1);
    time_array = nan(ball_count, 1);

    collision_ball_array = zeros(ball_count, 1);


    ball_patch = gobjects(ball_count, 1);

    on_table_id = 1: ball_count;
    enter_hole_tracker = zeros(ball_count, 1);

    % Setting Callbacks
    isDragging = false;
    isRunning = true;
    toPlaceBall = false;
    justStopped = false;

    stroke_ball_id = nan;
    potted_ball_ids = [];

    MotionTracker = zeros(ball_count, 1);

    %-----------------------------------------------

    ball_plot_x = zeros(16, ball_count);
    ball_plot_y = zeros(16, ball_count);

    ball_type_array = zeros(ball_count, 1);
    type2color = ["1 0 0", "1 0.99 0"];

    % iterate over the balls to setup their positions and colors
    color_array = zeros(ball_count, 3);
    color_array(1, :) = [1 1 1];
    for ii = 2: ball_count
        if ii > 6
            delta = 1;
        else
            delta = 0;
        end
        if mod(ii - delta, 2) == 0
            color_array(ii, :) = [1 0 0];
            ball_type_array(ii) = 1;
        else
            color_array(ii, :) = [1 0.99 0];
            ball_type_array(ii) = 2;
        end
    end

    if ball_count >= 7
        color_array(6, :) = [0 0 0];
        ball_type_array(6) = -1;
    end

    color_array = reshape(color_array, [ball_count, 1, 3]);

    player_names = strings(2, 1);
    player_colors = zeros(2, 1);

    disp("Welcome to the Pool game!" + newline + newline + ...
        "Your goal here is to pot the balls in the pockets using the cue (white) ball" + newline + ...
        "You'll be assigned a color to hit once one player pots a ball," + newline + ...
        "and you'll need to pot the black ball in the end to win the game" + newline + ...
        "Don't pot the black ball unless you potted all balls of the color you're assigned to," + newline + ...
        "otherwise you'll be considered losing the game" + newline + newline + ...
        "The following situations are considered as a foul" + newline + ...
        "1. you potted the cue (white) ball" + newline + ...
        "2. your cue ball didn't hit any other ball" + newline + ...
        "3. your cue ball collide with the other player's ball first after shooting (this is after you're assigned a color" + newline + newline + ...
        "In the case of a foul, the other player is allowed to freely place a ball behind the white line, " + newline + ...
        "So try not to make a foul" + newline + newline + ...
        "Now, please enter your names!")

    player_names(1) = input("Please enter the name of player 1: ", 's');
    player_names(2) = input("Please enter the name of player 2: ", 's');

    now_player = 1;
    winner = nan;
    score = zeros(2, 1);

    disp(newline + "Let the game begin!")

    % Setting up pool table
    f = figure(1);

    plot_table_w_hole();
    axis equal

    set(gca,'xtick',[],'ytick', []);

    hold on

    set(f, 'WindowButtonDownFcn', @OnMouseDown)
    set(f, 'WindowButtonUpFcn', @OnMouseUp)
    set(f, 'WindowButtonMotionFcn', @OnMouseMove);
    set(f, 'KeyPressFcn', @OnKeyDown);
    set(f, 'CloseRequestFcn', @OnFigureClose);

    init_setup();

    xlabel(player_names(now_player) + ": Drag and release to shoot the cue ball!")

    past_line = plot(0, 0);
    auxiliary_plot_array = gobjects(2, 1);

    max_time = 0;

    % the game iteration
    while isRunning
        tic;
        if sum(MotionTracker) ~= 0
            time_left = time_step;

            while min(time_array) < time_left
                [min_time, ~] = min(time_array);
                update_balls(min_time);
                time_left = time_left - min_time;
            end

            update_balls(time_left);

            for id = on_table_id
                add_ball_plot(id);
            end

            if sum(v_array(on_table_id)) == 0
                MotionTracker = zeros(ball_count, 1);
            end

            delete(ball_patch)
            ball_patch = patch(ball_plot_x(:, on_table_id), ball_plot_y(:, on_table_id), color_array(on_table_id, :, :), 'EdgeColor', "none");
        else
            if justStopped
                justStopped = false;
                foul = false;
                foul_reason = "";
                potted = false;
                now_player_reason = "";

                % detect fouls
                if sum(potted_ball_ids == 1) == 1 
                    foul = true;
                    foul_reason = player_names(now_player) + " potted the cue ball!";
                elseif isnan(stroke_ball_id)
                    foul = true;
                    foul_reason = player_names(now_player) + " didn't hit any ball!";
                elseif ball_type_array(stroke_ball_id) ~= player_colors(now_player) && sum(player_colors) ~= 0
                    foul = true;
                    foul_reason = player_names(now_player) + " hit the other player's ball first!";
                end

                if length(potted_ball_ids) >= 1
                    score_ball_ids = potted_ball_ids(potted_ball_ids ~= 1 & potted_ball_ids ~= 6);

                    % detect game-overs
                    if sum(potted_ball_ids == 6)
                        if score(now_player) == (ball_count - 2) / 2
                            winner = now_player;
                        else
                            winner = 2 - mod(now_player + 1, 2);
                        end
                        break
                    end

                    if ~isempty(score_ball_ids)
                        % assign colors
                        if sum(player_colors) == 0
                            player_colors(now_player) = ball_type_array(score_ball_ids(1));
                            player_colors(2 - mod(now_player + 1, 2)) = 2 - mod(ball_type_array(score_ball_ids(1)) + 1, 2);
                        end

                        if sum(ball_type_array(score_ball_ids) == player_colors(now_player)) > 0
                            potted = true;
                        end

                        score(1) = score(1) + sum(ball_type_array(score_ball_ids) == player_colors(1));
                        score(2) = score(2) + sum(ball_type_array(score_ball_ids) == player_colors(2));
                    end

                    if sum(player_colors) ~= 0
                        title("\color[rgb]{" + type2color(player_colors(1)) + "}" + player_names(1) + ": " + string(score(1)) + " \color{black}vs " + ...
                            "\color[rgb]{" + type2color(player_colors(2)) + "}" + player_names(2) + ": " + string(score(2)))
                    end
                end

                if ~potted || foul
                    now_player = 2 - mod(now_player + 1, 2);
                else
                    now_player_reason = "Wow! Nice shot! Keep it on!" + newline;
                end

                if foul
                    xlabel("Foul! " + foul_reason + newline + player_names(now_player) + ": click to place the cue ball" + newline + "on the right of the line")

                    if ismember(1, on_table_id)
                        delete_ball(1)
                        delete(ball_patch)
                        ball_patch = patch(ball_plot_x(:, on_table_id), ball_plot_y(:, on_table_id), color_array(on_table_id, :, :), 'EdgeColor', "none");
                    end
                    toPlaceBall = true;
                else
                    xlabel(now_player_reason + player_names(now_player) + ": Drag and release to shoot the cue ball!")
                end
            end
        end

        elapsed = toc;
        if elapsed > max_time
            max_time = elapsed;
        end

        pause(time_step - elapsed)
    end

    if ~isnan(winner)
        xlabel("Congratulations! The winner is " + player_names(winner) + "!")
    end
end

function d = still_ball_collision_distance(dp, v, r_sum)
    proj = dp * v';
    normal_v = dp - proj * v;
    normal = norm(normal_v);

    d = proj - sqrt(r_sum ^ 2 - normal ^ 2);

    if d < 0 || d ~= real(d)
        d = nan;
    end
end

function d = calc_d(s, a, t)
    d = s * t - a * t^2 / 2;
end

function flag = in_box_checker(pos, box)
    flag = (pos(1) < box(1) && pos(1) > box(3) && pos(2) < box(2) && pos(2) > box(4));
end