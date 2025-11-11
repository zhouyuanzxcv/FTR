function confusing_stages_as_subtypes_S13()
close all
dataset_names = {'ADNI_FSX_LM', 'ADNI_FSX_LS','ADNI_FSX_HM', 'ADNI_FSX_HS';
    'OASIS3_ROD1_LM', 'OASIS3_ROD1_LS','OASIS3_ROD1_HM', 'OASIS3_ROD1_HS';
    'NACC_LM', 'NACC_LS','NACC_HM', 'NACC_HS'};

method_names = {'FTR_MCEM', 'hierarch_clustering_baseline'};
method_labels = {'FTR (MCEM)', 'Hierarchical Clustering'};
nsubtype = 3;


% 定义总行数和列数
rows = 8;  % 总行数
cols = 6;  % 总列数

% 每个组内包含的行和列
group_rows = 2;  % 每组2行
group_cols = 2;  % 每组2列

% 设置组间距和组内间距
outer_gap = 0.04;  % 组间距
inner_gap = 0.01;  % 组内间距

% 计算每个组的宽度和高度
group_width = (1 - (cols/group_cols + 1) * outer_gap) / (cols / group_cols);
group_height = (1 - (rows/group_rows + 1) * outer_gap) / (rows / group_rows);

% 初始化figure窗口
f = figure;
% f.Position = [483 69 775 797];
f.Position = [483,1,816,865];
for idx_dataset = 1:size(dataset_names, 1)
    results = load_data_result(dataset_names(idx_dataset,:), method_names, nsubtype, 0);
    for idx_nbiom = 1:size(dataset_names,2)
        for idx_method = 1:length(method_names)
            group_col = idx_dataset;
            group_row = idx_nbiom;

            % 计算组的左下角位置
            group_left = (group_col - 1) * (group_width + outer_gap) + outer_gap;
            group_bottom = 1 - group_row * (group_height + outer_gap);

            % 计算子图在组内的位置
            sub_row = 0;
            sub_col = idx_method - 1;  % 组内列索引 (0 or 1)

            % 计算子图的宽度和高度
            sub_width = (group_width - inner_gap * (group_cols - 1)) / group_cols;
            sub_height = (group_height - inner_gap * (group_rows - 1)) / group_rows;

            % 计算子图左下角位置
            sub_left = group_left + sub_col * (sub_width + inner_gap);
            sub_bottom = group_bottom + (1 - sub_row) * (sub_height + inner_gap);

            % 创建axes
            ax = axes('Position', [sub_left, sub_bottom, sub_width, sub_height]);
            % colormap("sky");

            %% Count
            joindata = results{idx_nbiom, idx_method}.joindata;
            % [unique_RID, ia, ~] = unique(joindata.RID);
            % diagnosis = joindata.diagnosis(ia);
            % subtype = joindata.subtype(ia);
            diagnosis = joindata.diagnosis;
            subtype = joindata.subtype;
            edges1 = [0.5 1.5 2.5 3.5];
            edges2 = [-0.1 0.1 0.6 1.1];
            [counts,~,~] = histcounts2(subtype, diagnosis, edges1, edges2);
            imagesc(counts);
            colorbar;
            colormap("sky");
            clim([0 max(counts, [],"all")]);

            textStrings = num2str(counts(:), '%d');
            textStrings = strtrim(cellstr(textStrings));
            [xPos, yPos] = meshgrid(1:size(counts, 2), 1:size(counts, 1));
            text(xPos(:), yPos(:), textStrings, 'HorizontalAlignment', 'Center', 'Color', 'k');

            if group_row == 1 
                title('Diagnosis')
                % title(method_labels{idx_method});
            end
            % 

            ax_tmp = gca;
            ax_tmp.XAxisLocation = 'top';
            set(gca, 'XTick', 1:length(edges1)-1, 'XTickLabel', {'CN', 'MCI', 'AD'});

            
            if idx_method == 1
                set(gca, 'YTick', 1:length(edges1)-1, 'YTickLabel', {'S1', 'S2', 'S3'});
                if group_col == 1
                    ylabel('Count')
                end
            else
                set(gca, 'YTick', [])
            end
            % set(gca, 'XTick', 1:length(edges2)-1, 'XTickLabel', {'0', '0.5', '1'});
            % set(gca, 'YTick', 1:length(edges1)-1, 'YTickLabel', {'1', '2', '3'});
            % xlabel('第二列数据');
            % ylabel('第一列数据');
            % title('Subtype / Stage');
            % set(ax, 'XTick', [], 'YTick', []);
            % set(ax, 'XTick', [], 'YTick', []);

            %% Percentage
            sub_row = 1;
            sub_col = idx_method - 1;
            sub_width = (group_width - inner_gap * (group_cols - 1)) / group_cols;
            sub_height = (group_height - inner_gap * (group_rows - 1)) / group_rows;
            sub_left = group_left + sub_col * (sub_width + inner_gap);
            sub_bottom = group_bottom + (1 - sub_row) * (sub_height + inner_gap);
            ax = axes('Position', [sub_left, sub_bottom, sub_width, sub_height]);

            joindata = results{idx_nbiom, idx_method}.joindata;
            diagnosis = joindata.diagnosis;
            subtype = joindata.subtype;
            edges1 = [0.5 1.5 2.5 3.5];
            edges2 = [-0.1 0.1 0.6 1.1];
            [counts,~,~] = histcounts2(subtype, diagnosis, edges1, edges2);
            proportions = counts ./ sum(counts, 2);
            imagesc(proportions);
            colorbar;
            colormap("sky");
            clim([0 1]);

            textStrings = num2str(proportions(:), '%.2f');
            textStrings = strtrim(cellstr(textStrings));
            [xPos, yPos] = meshgrid(1:size(proportions, 2), 1:size(proportions, 1));
            text(xPos(:), yPos(:), textStrings, 'HorizontalAlignment', 'Center', 'Color', 'k');
            % set(gca, 'XTick', 1:length(edges2)-1, 'XTickLabel', {'0', '0.5', '1'});
            % set(gca, 'YTick', 1:length(edges1)-1, 'YTickLabel', {'1', '2', '3'});
            % xlabel('第二列数据');
            % ylabel('第一列数据');
            % title('Subtype / Stage');
            set(gca, 'XTick', [])
            % 
            if idx_method == 1
                set(gca, 'YTick', 1:length(edges1)-1, 'YTickLabel', {'S1', 'S2', 'S3'});
                if group_col == 1
                    ylabel('Proportion')
                end
            else
                set(gca, 'YTick', [])
            end
            % set(ax, 'XTick', [], 'YTick', []);
        end
    end
end

export_fig './figures/all/S13.jpg' -r500 -transparent
end