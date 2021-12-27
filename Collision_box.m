% Kelvin Li
% Simple 2D ball collision simulation inside a box without any pool setup

function Collision_box()
    function plot_table()
        rectangle('Position',[-margin, -margin, (y_length + margin) * 2, y_length + margin * 2], 'FaceColor',[0 0 0], 'Curvature', 0.2);
        rectangle('Position',[0, 0, y_length * 2, y_length], 'FaceColor',[0 0.75 0.75]);
        drawnow
        hold on
    end

    function add_ball_plot(id)
        th = linspace(0, 2*pi * (circle_sides - 1) / circle_sides, circle_sides);
        ball_plot_x(:, id) = (ball_r * cos(th) + pos_array(id, 1))';
        ball_plot_y(:, id) = (ball_r * sin(th) + pos_array(id, 2))';
    end

    function pos = get_mouse_pos()
        pos = get(gca, 'CurrentPoint');
        pos = pos(1, 1:2);
    end

    function OnMouseDown(~, ~)
        if sum(MotionTracker) == 0
            pos = get_mouse_pos();
            if (norm(pos - pos_array(1, :)) <= ball_r)
                isDragging = true;
            end
        end
    end

    function OnMouseUp(~, ~)
        if isDragging
            delete(past_line)
            isDragging = false;
            MotionTracker(1) = 1;
            CollisionManager([1]);
        end
    end

    function OnMouseMove (~, ~)
        if isDragging
            mouse_pos = get_mouse_pos();

            v = (mouse_pos - pos_array(1, :));
            v_array(1, :) = -1 * v / norm(v);

            coef = min(max_drag, norm(v)) / max_drag;

            speed_array(1) = max_speed * coef;
            color = (1 - coef) * (([255, 204, 207] / 255) - [1 0 0]) + [1 0 0];

            delete(past_line)
            drag_pos = pos_array(1, :) + (v / norm(v)) * min(max_drag, norm(v));
            past_line = plot([pos_array(1, 1) drag_pos(1)], [pos_array(1, 2) drag_pos(2)], "Color", color, 'LineWidth',2);
            hold on
        end
    end

    function OnFigureClose(~, ~)
        isRunning = false;
        delete(f)
    end

    function CollisionManager(update_array)
        test_ids = 1: ball_count;

        for update_id = update_array
            test_ids = test_ids(test_ids ~= update_id);

            [now_dists, now_times] = ball_collision_detection(update_id, test_ids, pos_array, v_array, speed_array, acc, ball_r, walls, true);
            [sort_time, sort_id] = sort(now_times);

            flag = 0;
            if ~isnan(sort_time(1))
                for count = 1: ball_count
                    ball_id = sort_id(count);
                    if ~isnan(now_times(ball_id)) && (sort_time(count) < time_array(ball_id) || isnan(time_array(ball_id)))
                        dist_array(update_id) = now_dists(ball_id, 1);
                        dist_array(ball_id) = now_dists(ball_id, 2);

                        time_array(update_id) = now_times(ball_id);
                        time_array(ball_id) = now_times(ball_id);

                        normal = (pos_array(update_id, :) + dist_array(update_id) * v_array(update_id, :)) ...
                                    - (pos_array(ball_id, :) + dist_array(ball_id) * v_array(ball_id, :));
                        normal = normal / norm(normal);

                        [v1, s1, v2, s2] = ball_collision(v_array(update_id, :), speed_array(update_id, :) - acc * time_array(update_id), ...
                                                            v_array(ball_id, :), speed_array(ball_id, :) - acc * time_array(ball_id), normal, ...
                                                            restitute_coef, tangent_coef);
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

                        flag = 1;
                        break
                    end
                end
            end

            if flag == 0
                if speed_array(update_id) > 0
                    [wall_id, d] = wall_collision_detection(pos_array(update_id, :), v_array(update_id, :), ball_r, walls);
    
                    c = [-acc/2, speed_array(update_id), -d];
                    
                    ts = roots(c);
                
                    t = nan;
                    for test_t = ts'
                        if (test_t >= 0 && test_t == real(test_t)) && (isnan(t) || test_t < t)
                            t = test_t;
                        end
                    end
    
                    if isnan(t)
                        t = speed_array(update_id) / acc;
                        dist_array(update_id) = d;
                        time_array(update_id) = t;
    
                        v_next_array(update_id, :) = [0, 0];
                        speed_next_array(update_id) = 0;
                    else
                        dist_array(update_id) = d;
                        time_array(update_id) = t;
        
                        if mod(wall_id, 2) == 1
                            normal = [1 0];
                        else
                            normal = [0 1];
                        end

%                         disp(string(update_id) + " will collide with wall")
%         
                        collision_speed = speed_array(update_id) - acc * t;
                        [v_next_array(update_id, :), speed_next_array(update_id)] = wall_collision(v_array(update_id, :), collision_speed, normal, restitute_coef, tangent_coef);
                    end
                end
            end
        end

        close_collision_record = max(0, close_collision_record - 0.04);
    end

    function update_balls(t)
        update_record = zeros(ball_count, 1);
        for id = 1: ball_count
            if time_array(id) >= 0
                if time_array(id) > t || ((speed_next_array(id) == 0) && (speed_array(id) <= acc * t + 1e-20))
                    if speed_array(id) ~= 0
                        if speed_array(id) > acc * t + 1e-20
                            now_d = calc_d(speed_array(id), acc, t);
                            speed_array(id) = speed_array(id) - acc * t;
                            time_array(id) = time_array(id) - t;
                        else
                            now_d = speed_array(id) ^ 2 / (2 * acc);
                            time_array(id) = nan;
                            speed_array(id) = 0;
                            MotionTracker(id) = 0;
    
                            update_record(id) = 1;
                        end
    
                        pos_array(id, :) = pos_array(id, :) + v_array(id, :) * now_d;
                    else
                        time_array(id) = time_array(id) - t;
                    end
                else
%                     disp(string(id) + " has undergone collision")

                    now_d_before = calc_d(speed_array(id), acc, time_array(id));
                    now_d_after = calc_d(speed_next_array(id), acc, t - time_array(id));
    
                    pos_array(id, :) = pos_array(id, :) + v_array(id, :) * now_d_before + v_next_array(id, :) * now_d_after;

                    v_array(id, :) = v_next_array(id, :);
                    speed_array(id) = speed_next_array(id) - acc * (t - time_array(id));

                    update_record(id) = 1;
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
            end
        end
        if sum(update_record) > 0
            update_array = [1: ball_count];
            update_array = update_array(update_record == 1);
            CollisionManager(update_array);
        end
    end


    % Settings
    y_length = 10;
    margin = 1.5;
    circle_sides = 16;

    ball_rows = 5;
    ball_count = 1 + sum(1: ball_rows);

    ball_r = 0.4;
    ball_margin = 0.1 * ball_r;

    ball_tri_pos = [0.75 * y_length, 0.5 * y_length];

    acc = 12;
    time_step = 0.025;
    max_speed = 48;
    max_drag = y_length / 2;

    restitute_coef = 0.8;
    tangent_coef = 0.8;

    walls = [2 * y_length, y_length, 0, 0];

    % Setting up pool table
    f = figure(1);

    plot_table();
    axis equal

    hold on

    % Setting up balls
    pos_array = zeros(ball_count, 2);
    v_array = zeros(ball_count, 2);
    speed_array = zeros(ball_count, 1);
    v_next_array = zeros(ball_count, 2);
    speed_next_array = zeros(ball_count, 1);
    
    dist_array = zeros(ball_count, 1);
    time_array = nan(ball_count, 1);

    pos_array(1, :) = [1.5 * y_length, 0.5 * y_length];

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

    % Setting Callbacks
    isDragging = false;
    isRunning = true;

    MotionTracker = zeros(ball_count, 1);

    past_line = plot(0, 0);

    set(f, 'WindowButtonDownFcn', @OnMouseDown)
    set(f, 'WindowButtonUpFcn', @OnMouseUp)
    set(f, 'WindowButtonMotionFcn', @OnMouseMove);
    set(f, 'CloseRequestFcn', @OnFigureClose);

    ball_plot_x = zeros(16, ball_count);
    ball_plot_y = zeros(16, ball_count);

    color_array = zeros(ball_count, 3);
    color_array(1, :) = [1 1 1];

    for ii = 2: ball_count
        color_array(ii, :) = [1 0 0];
    end
%     color_array(6, :) = [0 0 0];

    color_array = reshape(color_array, [ball_count, 1, 3]);

    for id = 1: ball_count
        add_ball_plot(id);
    end

    % Plotting takes much more time than calculating the path. Plotting
    % using "rectangle" would result in bad animation quality
    % Through experiment, I found that this is the quickest way to draw "circles"
    ball_patch = patch(ball_plot_x, ball_plot_y, color_array, 'EdgeColor', "none");

    disp("Ready to go!")

    while isRunning
        tic;
        if sum(MotionTracker) ~= 0
            time_left = time_step;
            while min(time_array) < time_left
                [min_time, ~] = min(time_array);
                update_balls(min_time);

                if break_flag == 1
                    break
                end

                time_left = time_left - min_time;
            end

            update_balls(time_left);

            for id = 1: ball_count
                add_ball_plot(id);
            end

            if break_flag == 1
                break
            end

            delete(ball_patch)
            ball_patch = patch(ball_plot_x, ball_plot_y, color_array, 'EdgeColor', "none");
        end

        elapsed = toc;
        title(elapsed)

        pause(time_step - elapsed)
    end
end


function d = calc_d(s, a, t)
    d = s * t - a * t^2 / 2;
end