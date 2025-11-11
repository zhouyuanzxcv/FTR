function plot_ARI_runtime_0A()
close all
noises = [0.01, 0.1, 0.2];
x = {5:5:30, 5:5:30, 5:5:20};
generation_method = {'Sigmoid', 'SuStaIn'};
generation_label = {'Sigmoid function','Piecewise linear function'};
clustering_method = {'FTR_MCEM','FTR_kmeans', 'sustain'};
color_area = {[70 240 240]./255, [220 190 255]./255, [250 190 212]./255};
color_line = {[0 0.4470 0.7410], [145 30 180]./255, [0.8500 0.3250 0.0980]};
shape_line = {'-',':x',':|'};

%% ARI
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
            c_method = clustering_method{idxC};
            y = [];
            for b = x{idxC}
                y_b = [];
                result_dir = ['output/Synthetic_data/synthetic_data_',g_method,'/',g_method,'_b',...
                    int2str(b), '_n', num2str(n)];
                for f = 1:5
                    result_path = [result_dir, '_f', int2str(f), '/', c_method, '_result.csv'];
                    t = readmatrix(result_path);
                    y_b =  [y_b, t(1)];
                end
                y = [y; y_b];
            end

            plot_uncertainty(x{idxC}, y', options)
            hold on;
            h = plot(x{idxC},mean(y,2),shape_line{idxC},'Color',color_line{idxC},'LineWidth',1,'MarkerSize',3);
            hs = [hs, h];
            xlim([5 30])
            ylim([0 1])
            xlabel('Number of biomarkers');
            if idxN == 1
                ylabel({generation_label{idxG},'Adjusted Rand index'});
            else
                ylabel('Adjusted Rand index');
            end
            if idxG == 1 && idxN == 1 % 仅在第一个 noise 和 generation_method 时保存 legend handles
                hs_all = hs;
            end
            if idxG == 1
                title({['Noise = ',num2str(n * 5)],''})
            end
        end
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
export_fig './figures/all/0A.jpg' -r500 -transparent
