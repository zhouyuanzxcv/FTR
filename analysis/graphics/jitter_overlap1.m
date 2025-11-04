function handles = jitter_overlap1(handles, b_spacing_multiplier, thresh)
if nargin < 3
    thresh = [0.3, 0.01];  % Balanced thresholds
end

if nargin < 2
    b_spacing_multiplier = 0.1;
end

y_positions = arrayfun(@(h) h.YData, handles);
thresh(2) = max(y_positions) * thresh(2);

% First, detect the natural spacing in the data to set appropriate thresholds
%     x_positions = arrayfun(@(h) h.XData, handles);
%     if length(unique(x_positions)) > 1
%         min_x_spacing = min(diff(sort(unique(x_positions))));
% %         thresh(1) = min(thresh(1), min_x_spacing * 0.8);  % Adaptive x-threshold
%     end

% Create groups based on both x-proximity and y-overlap potential
groups = create_xy_groups(handles, thresh);

% Process each group
for group_idx = 1:length(groups)
    group_indices = groups{group_idx};
    if length(group_indices) > 1
        group_handles = handles(group_indices);

        % Check if any bars in this group actually overlap with reasonable threshold
        %             if check_group_overlap(group_handles, thresh(2) * 0.5)  % Reduced threshold for actual check
        group_handles = jitter_group_balanced(group_handles, b_spacing_multiplier, thresh);
        handles(group_indices) = group_handles;
        %             end
    end
end
end





function groups = create_xy_groups(handles, thresh)
% Create groups with reasonable criteria
n_handles = length(handles);
x_positions = arrayfun(@(h) h.XData, handles);
y_intervals = arrayfun(@(h) get_y_interval(h, thresh(2)), handles, 'UniformOutput', false);

% Create adjacency matrix with reasonable criteria
adjacency = false(n_handles);
for i = 1:n_handles
    for j = i+1:n_handles
        % X-proximity check
        x_close = abs(x_positions(i) - x_positions(j)) <= thresh(1);

        % Y-overlap potential with reasonable threshold
        y_potential = intervals_overlap_single(y_intervals{i}, y_intervals{j});

        if x_close && y_potential
            adjacency(i, j) = true;
            adjacency(j, i) = true;
        end
    end
end

% Find connected components
groups = {};
visited = false(1, n_handles);

for i = 1:n_handles
    if ~visited(i)
        current_group = i;
        visited(i) = true;
        queue = find(adjacency(i, :));

        while ~isempty(queue)
            current = queue(1);
            queue(1) = [];

            if ~visited(current)
                current_group = [current_group, current];
                visited(current) = true;
                queue = [queue, find(adjacency(current, :))];
            end
        end

        groups{end+1} = current_group;
    end
end

% sort each group
for i = 1:length(groups)
    groups{i} = sort(groups{i});
end
end



function y_interval = get_y_interval(h, thresh)
y_interval = [h.YData - h.YNegativeDelta - thresh, ...
    h.YData + h.YPositiveDelta + thresh];
end

function overlap = intervals_overlap_single(A, B)
overlap = A(1) <= B(2) && B(1) <= A(2);
end

function group_handles = jitter_group_balanced(group_handles, spacing_multiplier, thresh)
n_bars = length(group_handles);

if n_bars == 1
    return; % No jittering needed for single bar
end

% Get original positions and y-intervals
original_x_positions = arrayfun(@(h) h.XData, group_handles);
y_centers = arrayfun(@(h) h.YData, group_handles);
y_intervals = arrayfun(@(h) get_y_interval(h, thresh(2)), group_handles, 'UniformOutput', false);
lengths = cellfun(@(ys) abs(ys(1) - ys(2)), y_intervals);

group_center_x = mean(original_x_positions);

% Sort by y-center for processing
%     [~, sort_idx] = sort(lengths, 'descend');
%     y_centers_sorted = y_centers(sort_idx);
%     y_intervals_sorted = y_intervals(sort_idx);

% Calculate minimum required x-separation for each pair
min_separations = calculate_min_separations(y_intervals, spacing_multiplier);

% Use constraint satisfaction to find optimal x-positions
x_positions = optimize_x_positions(y_centers, min_separations, ...
    group_center_x, spacing_multiplier);

% Apply the optimized x-positions
for i = 1:n_bars
    %         group_handles(sort_idx(i)).XData = x_positions(i);
    group_handles(i).XData = x_positions(i);
end
end

function min_separations = calculate_min_separations(y_intervals, base_spacing)
n_bars = length(y_intervals);
min_separations = zeros(n_bars, n_bars);

for i = 1:n_bars
    for j = i+1:n_bars
        if intervals_overlap_single(y_intervals{i}, y_intervals{j})
            % Bars overlap, need minimum separation
            overlap = min(y_intervals{i}(2), y_intervals{j}(2)) - max(y_intervals{i}(1), y_intervals{j}(1));
            min_separations(i, j) = 1;
            min_separations(j, i) = min_separations(i, j);
        else
            % Bars don't overlap, can be closer together
            gap = max(y_intervals{i}(1), y_intervals{j}(1)) - min(y_intervals{i}(2), y_intervals{j}(2));
            min_separations(i, j) = 0;
            min_separations(j, i) = min_separations(i, j);
        end
    end
end
end

function x_positions = optimize_x_positions(y_centers, min_separations, ...
    center_x, spacing_multiplier)

n_bars = length(y_centers);

% Initialize with greedy placement
x_positions = zeros(1, n_bars);
x_positions(1) = center_x;

% Place bars sequentially while respecting constraints
for i = 2:n_bars
    % Find minimum required position based on all previous bars
    x_candidates = unique(x_positions(1:i-1));

    find_empty_x = false;
    for j = 1:length(x_candidates)
        ind = (x_candidates(j) == x_positions(1:i-1));
        if ~any(min_separations(i,ind))
            x_positions(i) = x_candidates(j);
            find_empty_x = true;
            break
        end
    end

    if ~find_empty_x
        x_positions(i) = max(x_candidates) + spacing_multiplier/2;
    end
end

% Center the entire configuration around the original center
current_center = mean(x_positions);
x_positions = x_positions - (current_center - center_x);

% Optional: Refine with iterative relaxation
%x_positions = refine_positions(x_positions, min_separations, center_x);
end

function x_positions = refine_positions(x_positions, min_separations, center_x)
% Simple iterative relaxation to minimize spread
n_bars = length(x_positions);
max_iterations = 100;
tolerance = 1e-3;

for iter = 1:max_iterations
    moved = false;

    for i = 1:n_bars
        % Calculate forces from all other bars
        force = 0;
        for j = 1:n_bars
            if i ~= j
                current_distance = abs(x_positions(i) - x_positions(j));
                min_distance = min_separations(i, j);

                if current_distance < min_distance
                    direction = sign(x_positions(i) - x_positions(j));
                    force = force + direction * (min_distance - current_distance) * 0.5;
                end
            end
        end

        % Also add a force towards the center
        center_force = (center_x - x_positions(i)) * 0.1;

        total_force = force + center_force;

        if abs(total_force) > tolerance
            x_positions(i) = x_positions(i) + total_force;
            moved = true;
        end
    end

    if ~moved
        break;
    end
end
end

function required_spacing = calculate_required_spacing(group_handles, base_spacing)
% Calculate spacing needed to avoid overlaps with some tolerance
n_bars = length(group_handles);
y_intervals = zeros(n_bars, 2);

for i = 1:n_bars
    y_intervals(i, :) = get_y_interval(group_handles(i), 0);  % No threshold for actual calculation
end

required_spacing = base_spacing;

for i = 1:n_bars-1
    for j = i+1:n_bars
        if intervals_overlap_single(y_intervals(i, :), y_intervals(j, :))
            % Add spacing for each overlapping pair
            overlap = min(y_intervals(i, 2), y_intervals(j, 2)) - max(y_intervals(i, 1), y_intervals(j, 1));
            required_spacing = required_spacing + overlap * 0.3 + base_spacing * 0.2;
        end
    end
end

% Ensure minimum spacing for visual clarity
required_spacing = max(required_spacing, base_spacing * n_bars * 0.5);
end

function needs_jittering = check_group_overlap(group_handles, y_thresh)
% Check for overlap with reasonable threshold
n_bars = length(group_handles);
y_intervals = zeros(n_bars, 2);

for i = 1:n_bars
    y_intervals(i, :) = get_y_interval(group_handles(i), y_thresh);
end

for i = 1:n_bars
    for j = i+1:n_bars
        if intervals_overlap_single(y_intervals(i, :), y_intervals(j, :))
            needs_jittering = true;
            return;
        end
    end
end

needs_jittering = false;
end