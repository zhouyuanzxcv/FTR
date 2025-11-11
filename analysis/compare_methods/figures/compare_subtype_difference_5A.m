function compare_subtype_difference_5ABC()

close all



num_of_region = {'2','13','26','41','82'};
num_of_region = repmat(num_of_region, 2, 1);
num_of_region(3,:) = {'2','7','14','32','64'};

data = {'ADNI_FSX_LM', 'ADNI_FSX_LS','ADNI_FSX_HM', 'ADNI_FSX_HS';
    'OASIS3_ROD1_LM', 'OASIS3_ROD1_LS','OASIS3_ROD1_HM', 'OASIS3_ROD1_HS';
    'NACC_LM', 'NACC_LS', 'NACC_HM', 'NACC_HS'};

ctv_hv_ratio_data_idx = 1;

method = {'FTR_MCEM','sustain','hierarch_clustering_baseline','ctv_hv_ratio_baseline'};

hcbl = sprintf('Hierarchical\nClustering\n');
hcad = sprintf('Hierarchical\nClustering\n');
method_labels = {'FTR (MCEM)','SuStaIn',hcbl,'CTV-HV-Ratio'};
% method_colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330]};
method_colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560]};
method_shapes = {'o','square','diamond','v'};



nsubtype = 3;

num_runs = 5;

limit_logps = [6,5,4];


%% Fit linear mixed model
% Define the variables

variables = {'MMSE'};
% variables = {'ADAS13'};
% variables = {'MMSE','Abeta_summary','TAU_Temporal'};
f = figure();
f.Position = get_figure_position('fig51');

dependentVariable = {'MMSE'}; % Dependent variable (MMSE measurements)
dataset_title = {'ADNI (AD+MCI)', 'ADNI (AD)', 'OASIS'};
for idx_dataset = 1:size(data,1)
    subplot(size(data,1), 1, idx_dataset)

    % B维度数据的间隔增大倍数
    b_spacing_multiplier = 40;

    % C维度数据的间隔增大倍数
    c_spacing_multiplier = 170;

    handles = cell(1, length(method));

    limit_logp = limit_logps(idx_dataset);

    for include_control = [0,2]
        hold on;

        results = cell(size(data,2), length(method), num_runs);
        for idx_run = 1:num_runs

            options = [];
            options.postfix = ['_run',num2str(idx_run)];
    
            tmp = load_data_result(data(idx_dataset,:), method, nsubtype, ...
                include_control, options);

            results(:,:,idx_run) = tmp;
        end

        x_ticks = [];
        plot_handles = [];
        
        for dataIdx = 1:size(data,2)
            for methodIdx = 1:length(method)
                ftr_or_sustain = strcmp(method{methodIdx}, 'FTR_MCEM') || ...
                        strcmp(method{methodIdx}, 'sustain');
                % only apply FTR and SuStaIn on the validation sets
                if include_control == 2 && ~ftr_or_sustain
                    continue;
                end

                if isempty(results{dataIdx, methodIdx})
                    continue;
                end

                if strcmp(method{methodIdx}, 'ctv_hv_ratio_baseline') 
                    if dataIdx == ctv_hv_ratio_data_idx
                        x_value = -c_spacing_multiplier + b_spacing_multiplier;
                    else
                        continue;
                    end
                else
%                     x_value = (dataIdx-1)*c_spacing_multiplier + (methodIdx-1)*b_spacing_multiplier;
                    x_value = (dataIdx-1)*c_spacing_multiplier + (2-1)*b_spacing_multiplier;
                end

                

                if ftr_or_sustain % if ftr or sustain, use 5 runs
                    num_runs1 = num_runs;
                else % if not ftr or sustain, only 1 run is executed
                    num_runs1 = 1;
                end

                y_values = nan(3, num_runs1);

                for idx_run = 1:num_runs1

                    joindata = results{dataIdx, methodIdx, idx_run}.joindata;
                    mmse = results{dataIdx, methodIdx}.mmse;
                    covs = {'AGE_baseline', 'PTGENDER', 'PTEDUCAT', 'diagnosis_bl'};
                    independentVariables = [{'RID','years', 'subtype'}, covs];
                    joindata = joindata(:,[{'RID','subtype'},covs]);
                    joindata = unique(joindata);
                    joindata1 = outerjoin(joindata, mmse,'Type','Left','Keys','RID','MergeKeys',true);
                    joindata = joindata1;
                    dx_bl = joindata.diagnosis_bl;
                    data_sel_ind = {dx_bl == 0.5 | dx_bl == 1};
                    data_sel_name = {'MCI/AD at baseline'};
                    markersizes = [5,4,3];

                    for id_sel = 1:length(data_sel_ind)
                        try
                            p_sel = p_all_lmem(joindata, dependentVariable, ...
                                independentVariables, nsubtype, data_sel_ind{id_sel});
                            y_value = -log10(p_sel);
                            % jitter = (rand(length(data_sel_ind),1) - 0.5) * jitter_range;

                            y_value(y_value > limit_logp) = limit_logp;

                            y_values(:,idx_run) = y_value;
                            
                        catch me
                            % sustain on NACC_LM validation set only identified
                            % 2 subtypes
                            fprintf('WARNING: exception in calculate p value at run %s on %s when %s.\n', ...
                                method{methodIdx},data{idx_dataset,dataIdx},data_sel_name{id_sel});
                            msgText = getReport(me);
                            disp(msgText);
                        end

                    end
                end

                for id_sel = 1:length(data_sel_ind)
                    try
                        
                        for idx_point = 1:size(y_values, 1)
                            if include_control == 0
                                facecolor = method_colors{methodIdx};
                                linestyle = '-';
                            elseif include_control == 2
                                facecolor = 'none';
                                linestyle = '--';
                            end

                            y_value_m = mean(y_values(idx_point, :), 'omitnan');
                            y_value_std = std(y_values(idx_point, :), 'omitnan');

%                             h = plot(x_value, y_value_m, method_shapes{methodIdx}, ...
%                                 'Color', method_colors{methodIdx}, 'MarkerFaceColor', ...
%                                 facecolor,'MarkerEdgeColor',method_colors{methodIdx}, ...
%                                 'MarkerSize', markersizes(id_sel), 'linewidth', 1);

                            h = errorbar(x_value, y_value_m, y_value_std, ...
                                [linestyle,method_shapes{methodIdx}], ...
                                "MarkerSize",markersizes(id_sel),...
                                'linewidth', 1, 'color', method_colors{methodIdx}, ...
                                "MarkerEdgeColor", method_colors{methodIdx}, ...
                                "MarkerFaceColor", facecolor, ...
                                'CapSize',0);

                            if dataIdx == 1 && id_sel == 1 && idx_point == 1
                                plot_handles = [plot_handles, h];
                            end
                            handles{methodIdx} = [handles{methodIdx}, h];
                        end
                    catch me
                        % sustain on NACC_LM validation set only identified
                        % 2 subtypes
                        fprintf('WARNING: exception in calculate p value at run %s on %s when %s.\n', ...
                            method{methodIdx},data{idx_dataset,dataIdx},data_sel_name{id_sel});
                        msgText = getReport(me);
                        disp(msgText);
                    end

                end

            end
            x_ticks = [x_ticks, (dataIdx-1)*c_spacing_multiplier + b_spacing_multiplier];
        end
        x_ticks = [- c_spacing_multiplier + b_spacing_multiplier, x_ticks];
        set(gca, 'xtick', x_ticks, 'xticklabel', num_of_region(idx_dataset,:));
        set(gca, 'TickLabelInterpreter', 'none');
        xlabel('Number of regions')
        %     if idx_dataset == 1
        ylabel('-log_{10} p');
        %     end
        % title(data_title(p))
        hLine = yline(-log10(0.05), '--');
        set(get(get(hLine, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
        xlim([-40-c_spacing_multiplier (size(data,2)-1)*c_spacing_multiplier + (length(method)-1)*b_spacing_multiplier + 40])
        ylim([0 limit_logp + 0.05])

        yticklabels = [arrayfun(@num2str,0:limit_logp-1,'UniformOutput',false), ...
            {['\geq ' num2str(limit_logp)]}];

        set(gca, 'ytick',[0:limit_logp], 'yticklabel', yticklabels);
        set(gca, 'TickLabelInterpreter', 'tex');
        % title(dataset_title{idx_dataset})
        hold off;
    end
    
    handles = [handles{:}];

    % For debug
    x_positions = arrayfun(@(h) h.XData, handles);
    % for i = 1:length(handles), handles(i).XData = x_positions(i); end
    handles = jitter_overlap1(handles, b_spacing_multiplier);
end

% legend(plot_handles, method_labels, 'Location', 'westoutside', 'Interpreter', 'none');

export_fig './figures/all/5A.jpg' -r500 -transparent



end

function p_j = p_all_lmem(joindata, dependentVariable, ...
    independentVariables, nsubtype, data_sel_ind)
p_j = [];
for ref_group = 1:nsubtype
    data = joindata(data_sel_ind, [dependentVariable, independentVariables]);
    data = rmmissing(data);
    [~,~,~,beta_p_j] = lmem(data, ...
        dependentVariable, nsubtype, ref_group);
    p_j = [p_j,beta_p_j];
end
p_j = sort(rmmissing(p_j));
p_j = p_j(2:2:end)';
end

function plot_p(data, x_labels, b_labels, data_title, p)

subplot(1, 3, p)

[A, B, C] = size(data);
colors = lines(B);
hold on;
x_ticks = [];
plot_handles = zeros(1, B);
rng('default'); % For reproducibility of random numbers

% 设置随机扰动的范围，可以根据需要调整
jitter_range = 10;

% B维度数据的间隔增大倍数
b_spacing_multiplier = 20;

% C维度数据的间隔增大倍数
c_spacing_multiplier = 150;

for c = 1:C
    matrix = data(:, :, c);
    for b = 1:B
        values = matrix(:, b);
        for a = 1:A
            % 添加随机扰动到x坐标
            jitter = (rand - 0.5) * jitter_range;
            x_value = (c-1)*c_spacing_multiplier + (b-1)*b_spacing_multiplier + jitter;
            h = plot(x_value, values(a), 'o', 'Color', colors(b, :), 'MarkerFaceColor', colors(b, :), 'MarkerSize', 6);
            if c == 1 && a == 1
                plot_handles(b) = h;
            end
        end
    end
    x_ticks = [x_ticks, (c-1)*c_spacing_multiplier + (B-1)*b_spacing_multiplier / 2];
end

% 设置x轴的ticks和labels
set(gca, 'xtick', x_ticks, 'xticklabel', x_labels);
set(gca, 'TickLabelInterpreter', 'none');

legend(plot_handles, b_labels, 'Location', 'best', 'Interpreter', 'none');

% title('Test');
xlabel('Number of regions')
ylabel('Subtype difference in MMSE slope (-log_{10} p)');
title(data_title(p))
hLine = yline(-log10(0.05), '--');
set(get(get(hLine, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

hold off;

end