function [outputArg1,outputArg2] = analyze_longitudinal_stability_6A()

close all

dataset_names = {'ADNI_FSX_LM', 'ADNI_FSX_LS','ADNI_FSX_HM', 'ADNI_FSX_HS';
    'OASIS3_ROD1_LM', 'OASIS3_ROD1_LS','OASIS3_ROD1_HM', 'OASIS3_ROD1_HS';
    'NACC_LM', 'NACC_LS', 'NACC_HM', 'NACC_HS'};

num_of_region = {'13','26','41','82'};
num_of_region = cat(1, repmat(num_of_region, 2, 1), {'7','14','32','64'});

method_names = {'FTR_MCEM','sustain'};


nsubtype_list = [3];
idx_subtype = 1; % use 3 subtypes

data_splits = {'last_point','baseline'};

save_path = ['longitudinal_stability'];

% dataset_linestyle = {'-','--',':','-.'};

% validate = 0;



%% plot using heatmaps

if 1
%     load('tmp_ys_all.mat');
    [metrics, ys_all, validate] = retrieve_longitudinal_stability(dataset_names, ...
        method_names, nsubtype_list, data_splits, save_path, idx_subtype);
    save('tmp_ys_all.mat', "ys_all");
else
    load('tmp_ys_all.mat');
end

% ys((dataset-1)*2+method, dim, run, split, metric, validate)

method_labels = {'FTR (MCEM)', 'SuStaIn'};
method_colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980]};

dataset_names1 = {'ADNI','OASIS','NACC'};

metric_names = {'Consistency','Adjusted Rand index'};

split_val = [1,1;2,1;1,2;2,2];

for i = 1:size(split_val,1)
    fh = figure;
    fh.Position = get_figure_position(10);
    plot_bar_error_comparison_3d(squeeze(ys_all(:,:,:,split_val(i,1),:,split_val(i,2))), method_labels,...
        method_colors, num_of_region, dataset_names1, metric_names);
    export_fig(['./figures/all/6A_',int2str(split_val(i,1)),'_', ...
        int2str(split_val(i,2)),'.jpg'],'-r500','-transparent');
end


end




function [metrics, ys_all, validate] = retrieve_longitudinal_stability(dataset_names, ...
    method_names, nsubtype_list, data_splits, save_path, idx_subtype, ys_all)

metrics = {'consistency','ARI'};

num_splits = length(data_splits);

data_inds_AD = 1:size(dataset_names, 2);

num_runs = 5;

if nargin < 7
    ys_all = nan(size(dataset_names,1)*length(method_names), ...
        size(dataset_names, 2), num_runs, num_splits, length(metrics),2);
end


for validate = [0,1]
    for idx_split = 1:num_splits
      
            num_subjs = [];
            
            for idx_dataset = 1:size(dataset_names,1)
                
                for run_idx = 1:num_runs

                    % the dimension of consistency is split_idx, data_idx, method_idx,
                    % subtype_idx
                    [confusion_matrices, percentage_consistency, run_time] = ...
                        run_algorithms_for_longitudinal_stability(...
                        dataset_names(idx_dataset,:), method_names, nsubtype_list, ...
                        data_splits(idx_split), save_path, validate, ['_run',num2str(run_idx)]);                   
                    
                    for idx_metric = 1:length(metrics)
                        metric = metrics{idx_metric};

                        curr_split_idx = 1;
                        
                        if strcmp(metric, 'consistency')
                            y = percentage_consistency(curr_split_idx, ...
                                data_inds_AD, :, idx_subtype);
                        elseif strcmp(metric, 'ARI')
                            y = [];
                            num_subj = [];
                            cms = squeeze(confusion_matrices(curr_split_idx, ...
                                data_inds_AD, :, idx_subtype));
                            for idx1 = 1:size(cms,1) % for each dimension
                                for idx2 = 1:size(cms,2) % for each method
                                    cm = cms{idx1,idx2};
                                    if isempty(cm) || any(isnan(cm),'all')
                                        y(idx1,idx2) = nan;
                                    else
                                        y(idx1,idx2) = rand_index_from_contingency(cm);
                                        num_subj(idx1,idx2) = sum(cm(:));
                                    end
                                end
                            end
                            
                            num_subj = squeeze(num_subj);
                            if run_idx == 1
                                num_subjs(end+1:end+2, :) = num_subj';
                            end
                        end
                        
                        y = squeeze(y);
                        rows_ind_ys = (idx_dataset-1)*length(method_names)+1: ...
                            idx_dataset*length(method_names);
%                         ys(rows_ind_ys,:,run_idx) = y';

                        ys_all(rows_ind_ys,:,run_idx,idx_split,idx_metric,validate+1) = y';
                    end
                end

            end

    end
    
end

% ys_all = cell2mat(ys_all);
end

