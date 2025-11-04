function [outputArg1,outputArg2] = calibrate_stability(inputArg1,inputArg2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
close all

% dataset_name = 'ADNI_FSX_HS';
% dataset_name = 'OASIS3_ROD1_HS';
dataset_name = 'NACC_HS';

dataset_folder = [dataset_name, '_FTR_MCEM'];

[data,mmse,time2event,dx_sel] = load_data_all(dataset_name);
PTID = data{data.group == 1, 'RID'};


for i = 1:5
    file_path = ['output/',dataset_folder,'/3_subtypes_run',num2str(i),'/init_runs.mat'];
    S = load(file_path);
    init_runs = S.init_runs;
    init_runs.RID = PTID(init_runs.inds_sel);
    results(i) = init_runs;
end



figure;
tiledlayout(2,3);

subtypes_sel = {};

metric_all_mode = [];
metric_all_meta = [];



for run_ind = 1:5
    nexttile;
    subtype = results(run_ind).subtype;
    sigma = results(run_ind).sigma;
    best_idx = results(run_ind).best_idx;
%     RID = results(run_ind).RID;
    
%     [~,ia] = unique(RID);
%     subtype1 = subtype(ia,:);
    subtype1 = subtype;

    meta_labels_old = ones(100,1);

    % run_ind
    for try_index = size(subtype,2)/100
    % for try_index = size(subtype1,2)/100
        subtype_try = subtype1(:, 1:try_index*100);
        thresh_same_partition = 0.98;
        [similarity, n_meta, meta_labels, proption, metric, bes_ep1, mean_ARIs] =  ...
            sel_init_by_meta_clustering(subtype_try, sigma, thresh_same_partition);
        % bes_ep1
        ari = rand_index(meta_labels(1:length(meta_labels_old)), meta_labels_old);
        meta_labels_old = meta_labels;
        if ari == 1
            break;
        end
    end
    fprintf('run_ind %d, first %d * 100, ARI %.2f \n', run_ind, try_index, ari);
    meta_labels1 = nan(size(subtype1,2),1);
    meta_labels1(1:length(meta_labels)) = meta_labels;
    meta_labels = meta_labels1;

    metric_all_meta(run_ind) = metric;




%     [best_partition, bes_ep1] = sel_init_by_clique(subtype, ...
%         thresh_same_partition);

    options = [];
    options.thresh_same_partition = 0.98;
    [bes_ep2, similarity, bes_subtype, eps, subtype1, bes_sgm, metric] = ...
        sel_init_by_partition_mode(subtype1, results(run_ind).sigma, options);
%     bes_ep1 = bes_ep2;
    metric_all_mode(run_ind) = metric;


    subtypes_sel{run_ind} = subtype(:, bes_ep1);
%     subtypes_sel{run_ind} = subtype_avg;

    
    [coeff, score] = pca(similarity);
%     score = tsne(similarity, 'NumDimensions', 2);

    colors = distinguishable_colors(n_meta);
    
    h = gscatter(score(:,1), score(:,2), meta_labels, colors, 'osd', 2, 'filled');
    h_scatter = findobj(h, 'Type', 'scatter');
    set(h_scatter, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
    hold on

    h1 = scatter(score(best_idx, 1), score(best_idx, 2), 'filled', ...
        'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1);
    h2 = scatter(score(bes_ep1, 1)+0.1, score(bes_ep1, 2)+0.1, 'filled', ...
        'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1);
    h3 = scatter(score(bes_ep2, 1)-0.1, score(bes_ep2, 2)-0.1, 'filled', ...
        'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1);


    legend_labels = cell(length(proption), 1);
    for k = 1:length(proption)
        legend_labels{k} = sprintf('%d (%.1f%%), Mean ARI %.2f', ...
            k, proption(k)*100, mean_ARIs(k));
    end
    legend_labels{n_meta+1} = 'previous';
    legend_labels{n_meta+2} = 'chosen';
    legend_labels{n_meta+3} = 'mode';

    legend([h;h1;h2;h3], legend_labels, 'Location', 'best', 'FontSize', 10);
    
    title('PCA of ARI Matrix');
    xlabel('PC1');
    ylabel('PC2');
end

disp('ARI or -consistency between choosen partitions')
ARI = [];
common_inds_len = [];
for try_ind1 = 1:5
    for try_ind2 = 1:5
        if isfield(results(try_ind1), 'inds_sel')
            common_inds_12 = results(try_ind1).inds_sel & results(try_ind2).inds_sel;
            common_inds_len(try_ind1, try_ind2) = length(find(common_inds_12));

            subtype_1 = sel_from_common(subtypes_sel{try_ind1}, ...
                results(try_ind1).inds_sel, common_inds_12);

            subtype_2 = sel_from_common(subtypes_sel{try_ind2}, ...
                results(try_ind2).inds_sel, common_inds_12);

            ARI(try_ind1, try_ind2) = rand_index(subtype_1, subtype_2);
        else

            if length(subtypes_sel{try_ind1}) == length(subtypes_sel{try_ind2})
                c1 = subtypes_sel{try_ind1};
                c2 = subtypes_sel{try_ind2};
                ARI(try_ind1, try_ind2) = rand_index(c1, c2);
            else
                ARI(try_ind1, try_ind2) = -calc_consistency_from_mismatched_partitions( ...
                    subtypes_sel{try_ind1}, subtypes_sel{try_ind2});
            end
        end
        
    end
end
ARI

metric_all_mode
metric_all_meta

figure, histogram(similarity(:));
       

end

function subtype_common = sel_from_common(subtype, inds_sel, common_inds)
subtype_original = nan(size(inds_sel));
subtype_original(inds_sel) = subtype;
subtype_common = subtype_original(common_inds);
end
