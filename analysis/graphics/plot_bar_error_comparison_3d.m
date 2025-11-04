function plot_bar_error_comparison_3d(tensor_data, method_names, method_colors, xticklabels_cell, dataset_names, metric_names)
% PLOT_BAR_ERROR_COMPARISON_3D Creates bar plots comparing methods across datasets, dimensions, and metrics
%
% Inputs:
%   tensor_data: 4D tensor [datasets*methods x dimensions x runs x metrics]
%   method_names: cell array of method names
%   method_colors: cell array or matrix of RGB colors for each method
%   xticklabels_cell: 3x4 cell array with xticklabels for each dataset
%   dataset_names: cell array of dataset names
%   metric_names: cell array of metric names (length 2)

% Get dimensions
num_datasets = size(tensor_data, 1) / length(method_names);
num_methods = length(method_names);
num_dims = size(tensor_data, 2);
num_runs = size(tensor_data, 3);
num_metrics = size(tensor_data, 4);

% Validate inputs
validate_input(tensor_data, num_datasets, num_methods, method_colors,...
    num_metrics, xticklabels_cell, num_dims, dataset_names, metric_names);

% Calculate means and standard deviations (handle NaN values)
[means, stds] = calc_stats(num_datasets, num_methods, num_dims, num_metrics, tensor_data);




% Create main tiled layout with adjusted spacing
main_tile = tiledlayout(num_datasets, num_metrics, 'TileSpacing', ...
    'compact', 'Padding', 'compact');



% Create plots for each dataset and metric combination
for dataset = 1:num_datasets
    for metric = 1:num_metrics
        nexttile((dataset-1)*num_metrics + metric);
        hold on;

        % Prepare data for this dataset and metric
        bar_data = squeeze(means(dataset, :, :, metric))'; % [dims × methods]
        error_data = squeeze(stds(dataset, :, :, metric))'; % [dims × methods]

        % Create grouped bar plot
        h = bar(bar_data, 'grouped');

        % Set colors for each method and add error bars
        for method = 1:num_methods
            h(method).FaceColor = method_colors{method};
            h(method).EdgeColor = 'none';
            h(method).FaceAlpha = 0.5;

            [ngroups, nbars] = size(bar_data);
            groupwidth = min(0.8, nbars/(nbars + 1.5));

            x = (1:ngroups) - groupwidth/2 + (2*method-1) * groupwidth / (2*nbars);

            

            for j = 1:num_dims
                y = squeeze(tensor_data((dataset - 1) * num_methods + method, j, :, metric));
                swarmchart(repmat(x(j),num_runs,1), y, 20, method_colors{method}, 'filled', ...
                    'MarkerFaceAlpha', 0.6, 'XJitterWidth', 0.2);
            end

            % Add error bars for this method
            if num_runs > 1 % Only add error bars if there are multiple runs
                
                errorbar(x, bar_data(:, method), error_data(:, method), ...
                    'k.', 'LineWidth', 1, 'CapSize', 4, 'HandleVisibility', 'off');
            end
        end

        % Customize plot appearance
        set(gca, 'XTick', 1:num_dims);
        set(gca, 'XTickLabel', xticklabels_cell(dataset, :));
        xlim([0.5, num_dims + 0.5]);
        grid on;
        box on;

        % Add x-label only for bottom row
        if dataset == num_datasets
            xlabel('Number of regions');
        end

        % Add title only for top row (not bold)
        if dataset == 1
            title({'','',metric_names{metric}}, 'fontweight', 'normal');
        end

        % Add dataset name as y-label for first column only
        if metric == 1
            ylabel(dataset_names{dataset});
            ylim([0.5,1]);
        elseif metric == 2
            ylim([0,1]);
        end
    end
end

add_legend(main_tile, method_names);

end

function validate_input(tensor_data, num_datasets, num_methods, ...
    method_colors, num_metrics, xticklabels_cell, num_dims, dataset_names, metric_names)
if size(tensor_data, 1) ~= num_datasets * num_methods
    error('Number of rows in tensor must be num_datasets * num_methods');
end

if length(method_colors) ~= num_methods
    error('Number of colors must match number of methods');
end

if num_metrics ~= 2
    error('Number of metrics must be 2');
end

if ~isequal(size(xticklabels_cell), [num_datasets, num_dims])
    error('xticklabels_cell must be a %dx%d cell array', num_datasets, num_dims);
end

if length(dataset_names) ~= num_datasets
    error('Number of dataset names must match number of datasets');
end

if length(metric_names) ~= 2
    error('Number of metric names must be 2');
end
end

function [means, stds] = calc_stats(num_datasets, num_methods, num_dims, num_metrics, tensor_data)
means = zeros(num_datasets, num_methods, num_dims, num_metrics);
stds = zeros(num_datasets, num_methods, num_dims, num_metrics);

for dataset = 1:num_datasets
    for method = 1:num_methods
        row_idx = (dataset - 1) * num_methods + method;
        for metric = 1:num_metrics
            for dim = 1:num_dims
                data_slice = squeeze(tensor_data(row_idx, dim, :, metric));
                valid_data = data_slice(~isnan(data_slice));

                if ~isempty(valid_data)
                    means(dataset, method, dim, metric) = mean(valid_data);
                    stds(dataset, method, dim, metric) = std(valid_data);
                else
                    means(dataset, method, dim, metric) = NaN;
                    stds(dataset, method, dim, metric) = NaN;
                end
            end
        end
    end
end
end

function add_legend(main_tile, method_names)
% Get position of top-left subplot
top_left_ax = main_tile.Children(end); % First created axes is at the end
top_left_pos = top_left_ax.Position;
% 
% Create legend in the north tile but manually position it
lg = legend(method_names, 'Orientation', 'horizontal');

% Calculate proper position for left alignment
fig = gcf;
fig_pos = fig.Position;
ax_pos = top_left_ax.Position;

% Convert normalized positions to pixels
ax_left_pixels = ax_pos(1) * fig_pos(3);
ax_width_pixels = ax_pos(3) * fig_pos(3);

% Position legend left-aligned with top-left subplot
lg.Position = [ax_left_pixels/fig_pos(3), ... % Left position (normalized)
    0.91, ... % Top position (normalized)
    0.6, ... % Width (normalized)
    0.05];   % Height (normalized)
end