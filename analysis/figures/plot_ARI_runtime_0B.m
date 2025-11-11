function plot_ARI_runtime_0B()
close all
noises = [0.01, 0.1, 0.2];
x = {5:5:30, 5:5:30, 5:5:20};
generation_method = {'Sigmoid', 'SuStaIn'};
generation_label = {'Sigmoid function','Piecewise linear function'};
clustering_method = {'FTR_MCEM','FTR_kmeans', 'sustain'};
color_area = {[70 240 240]./255, [220 190 255]./255, [250 190 212]./255};
color_line = {[0 0.4470 0.7410], [145 30 180]./255, [0.8500 0.3250 0.0980]};
shape_line = {'-',':x',':|'};


%% runtime
close all
figure;
hold on
nsubp = 4;
options = [];
hs_all = [];
for idxG = 1:length(generation_method)
    g_method = generation_method{idxG};
    idxN = 0;
    for n = noises
        idxN = idxN + 1;
        subplot(3,3,nsubp);
        nsubp = nsubp + 1;
        hs = [];
        for idxC = 1:length(clustering_method)
            options.color_area = color_area{idxC};
            options.color_line = color_line{idxC};
            options.plot_type = 'semilogx';
            options.line_shape = shape_line{idxC};
            options.line_width = 1;
            c_method = clustering_method{idxC};
            y = [];
            for b = x{idxC}
                y_b = [];
                result_dir = ['output/Synthetic_data/synthetic_data_',g_method,'/',g_method,'_b',...
                    int2str(b), '_n', num2str(n)];
                for f = 1:5
                    result_path = [result_dir, '_f', int2str(f), '/', c_method, '_result.csv'];
                    t = readmatrix(result_path);
                    y_b =  [y_b, log10(t(2))];
                end
                y = [y; y_b];
            end

            h = plot_uncertainty(x{idxC}, y', options);
            hold on;
            % h = plot(x{idxC},mean(y,2),shape_line{idxC},'Color',color_line{idxC},'LineWidth',1,'MarkerSize',3);
            hs = [hs, h];

            xlabel('Number of biomarkers');
            if idxN == 1
                ylabel({generation_label{idxG},'Running time (min)'});
            else
                ylabel('Running time (min)');
            end
            if idxG == 1 && idxN == 1 % 仅在第一个 noise 和 generation_method 时保存 legend handles
                hs_all = hs;
            end
            if idxG == 1
                title({['Noise = ',num2str(n * 5)],''})
            end
        end
        ylim([0 3.2])
        yticks([0,1,2,3])
        yticklabels({'10^0','10^1','10^2','10^3'})
    end
end
% 获取 subplot(3,3,4) 和 subplot(3,3,6) 的位置
pos_left = get(subplot(3,3,4), 'Position');
pos_right = get(subplot(3,3,6), 'Position');

% 创建 legend 并设置位置
leg = legend(hs_all, {'FTR (MCEM)','FTR (kmeans)','SuStaIn'}, 'Orientation', 'horizontal');
legend_position = get(leg, 'Position');
legend_position(1) = pos_left(1); % 左对齐到 subplot(3,3,4)
legend_position(3) = pos_right(1) + pos_right(3) - pos_left(1); % 宽度跨越到 subplot(3,3,6)
legend_position(2) = legend_position(2) + 0.13; % 向上移动 0.05（可以根据需要调整）
set(leg, 'Position', legend_position);

% subplot(3,3,2)
% legend(hs, {'FTR (MCEM)','FTR (kmeans)','SuStaIn'},'Location','bestoutside','Orientation','horizontal');
export_fig './figures/all/0B.jpg' -r500 -transparent

end

function h = plot_uncertainty(X, Y, options)
% Input:
%   X - 1 x N vector
%   Y - M x N matrix or 1 x N cell array

%% Default options
if nargin < 3
    options = [];
end

handle = parse_param(options, 'figure', gcf);
%color_area = parse_param(options, 'color_area', [128 193 219]./255);    % Blue theme
%color_line = parse_param(options, 'color_line', [ 52 148 186]./255);
%options.color_area = [243 169 114]./255;    % Orange theme
%options.color_line = [236 112  22]./255;
color_area = parse_param(options, 'color_area', [169,169,169]./255);    % Black theme
color_line = parse_param(options, 'color_line', [0 0 0]./255);
alpha = parse_param(options, 'alpha', 0.5);
line_width = parse_param(options, 'line_width', 2);
line_shape = parse_param(options, 'line_shape', '-');
marker_size = parse_param(options, 'marker_size', 3);
center_type = parse_param(options, 'center_type', 'mean');
error_type = parse_param(options, 'error_type', 'std');

plot_type = parse_param(options, 'plot_type', 'plot');

%%
if size(X,2) == 1
    Y = Y';
    X = X';
end

if ~iscell(Y) % Y is M x N
    Y1 = {};
    for i = 1:size(Y,2)
        Y1{i} = Y(:,i);
    end
    Y = Y1;
end


%% Computing the mean and standard deviation of the data matrix
data_center = [];

if strcmp(center_type, 'mean')
    for i = 1:length(Y)
        data_center(:,i) = mean(Y{i},1);
    end
elseif strcmp(center_type, 'median')
    for i = 1:length(Y)
        data_center(:,i) = median(Y{i},1);
    end
end

for i = 1:length(Y)
    data_std(:,i)  = std(Y{i},0,1);
end

% Type of error plot
switch(error_type)
    case 'std'
        error = data_std;
        upper_bd = data_center + error;
        bottom_bd = data_center - error;
    case 'q1-q3'
        for i = 1:length(Y)
            upper_bd(:,i) = quantile(Y{i}, 0.75);
            bottom_bd(:,i) = quantile(Y{i}, 0.25);
        end
    otherwise
end

% Plotting the result
figure(handle);

if strcmp(plot_type, 'plot')
    % plot(X, data_center, 'color', color_line, 'LineWidth', line_width);
elseif strcmp(plot_type, 'semilogx')
    h = semilogx(X, data_center, line_shape, 'color', color_line, 'linewidth', line_width, 'MarkerSize', marker_size);
    grid on
elseif strcmp(plot_type, 'loglog')
    h = loglog(X, data_center, line_shape, 'color', color_line, 'linewidth', line_width, 'MarkerSize', marker_size);
    grid on
end
hold on;

x_vector = [X, fliplr(X)];
patch = fill(x_vector, [upper_bd,fliplr(bottom_bd)], color_area);
set(patch, 'edgecolor', 'none');
set(patch, 'FaceAlpha', alpha);
hold off;

% send patch object to background
set(gca,'children',flipud(get(gca,'children')))

end
