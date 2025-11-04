function [best_idx, is_stable, best_metric, options] = sel_from_multi_inits(subtype1, sgm1, ...
    loglik_all1, subtype_all1, options)
%SEL_FROM_MULTI_INITS Summary of this function goes here
%   Detailed explanation goes here

if size(subtype1, 2) == 1
    best_idx = 1;
    is_stable = 1;
    best_metric = 1;
    return;
end

% filter out those failed cases (one cluster becomes empty)
inds1 = find(~any(isnan(subtype1), 1));

num_runs = size(subtype1, 2);
if length(inds1) < num_runs
    fprintf('%d out of %d runs failed.\n', num_runs - length(inds1), num_runs);
end

subtype = subtype1(:,inds1);
sgm = sgm1(:,:,inds1);
loglik_all = loglik_all1(:,inds1);
subtype_all = subtype_all1(:,:,inds1);

% verify_stability = parse_param(options, 'verify_stability', 0);
% 
% if verify_stability && size(subtype, 2) > 1
%     % On ADNI, the strategy that selects the largest log-likelihood has a median RI of
%     % 0.52 while the strategy that selects the most common partition has a
%     % median RI of 0.97 over 5 runs, each with 20 random initializations
%     % Update on 12/27/2024: On OASIS, using the 50% most converged partitions
%     % as initializations leads to better ARI (median ARI 0.92) over 100
%     % inits
%     % Update on 1/8/2025: On OASIS, averaging gives more stable results than
%     % selecting the most common partition (median ARI 0.97 vs 0.79) over 100 inits.
%      
%     close all
%     
%     sims = [];
%     parfor run_idx = 1:30
%         [sim_mat1, sim_mat2] = compare_strategy_stability(subtype, sgm, ...
%             loglik_all, subtype_all, options);
%         sims(run_idx, :) = [sim_mat1(1,2), sim_mat2(1,2)];
%     end
%     figure, boxplot(sims);
%     disp(['Two strategies to select from multiple initializations ', ...
%         '(row 1: partition, row 2: average, cols: median, Q1, Q3)']);
%     [median(sims(:,1)), prctile(sims(:,1),25), prctile(sims(:,1),75); ...
%         median(sims(:,2)), prctile(sims(:,2),25), prctile(sims(:,2),75)]
% end


partition_selection_method = parse_param(options, 'partition_selection_method', ...
    'meta_clustering');

if strcmp(partition_selection_method, 'partition_mode')
    [best_idx, similarity_mat, bes_subtype, eps, subtype, bes_sgm, best_metric] = ...
        sel_init_by_partition_mode(subtype, sgm, options);
    num_same_partition = parse_param(options, 'num_same_partition', 20);
    if best_metric >= num_same_partition
        is_stable = 1;
    else
        is_stable = 0;
    end
    options.meta_proption = 1;
    options.similarity = similarity_mat;
    options.meta_metric = best_metric;
    
elseif strcmp(partition_selection_method, 'meta_clustering')
    thresh_same_partition = parse_param(options, 'thresh_same_partition', 0.98);
    [similarity_mat, n_meta, meta_labels, proption, best_metric, best_idx, mean_ARIs] = ...
        sel_init_by_meta_clustering(subtype, sgm, thresh_same_partition);
    
%     fprintf('Meta-clustering proportion: %s, number of same partitions: %d. \n', ...
%         mat2str(proption, 2), best_metric);

    % update meta_labels_old for checking stability convergence
    meta_labels_old = parse_param(options, 'meta_labels_old', ones(size(meta_labels)));
    ari = rand_index(meta_labels(1:length(meta_labels_old)), meta_labels_old);
    options.meta_labels_old = meta_labels;
    options.meta_proption = proption;
    options.meta_mean_ARIs = mean_ARIs;
    options.similarity = similarity_mat;
    options.meta_metric = best_metric;

    num_same_partition = parse_param(options, 'num_same_partition', 20);
    if 0
        is_stable = (ari == 1);
        best_metric = ari;
    else
        if best_metric >= num_same_partition
            is_stable = 1;
        else
            is_stable = 0;
        end
    end


%     if strcmp(options.methods, 'kmeans')
%         pre_subtype = subtype_avg;
%         pre_sigma = sgm_avg;
%     elseif strcmp(options.methods, 'MCEM')
%         pre_subtype = subtype_perms;
%     	pre_sigma = sgm_avg;
%     end
% 
%     options.pre_subtype = pre_subtype;
%     options.pre_sgm = pre_sigma;
end

best_idx = inds1(best_idx);

% [bes_subtype, bes_sgm, subtype_perms, sgm_perms] = ...
%     sel_based_on_partition_averaging(subtype, sgm);




% sel_strategy = parse_param(options, 'select_from_multi_inits', 'average');
% % use one strategy
% if strcmp(sel_strategy, 'loglikelihood')
%     bes_ep = sel_based_on_loglik(loglik_all);
% elseif strcmp(sel_strategy, 'partition')
%     bes_ep = sel_based_on_partition(subtype);
% elseif strcmp(sel_strategy, 'partition_hierarchical')
%     bes_ep = sel_based_on_partition_hierarchical(subtype);
% elseif strcmp(sel_strategy, 'part_w_loglik')
%     bes_ep = sel_based_on_weighted_loglik(subtype, loglik_all);
% end
% 
% pre_subtype = subtype(:,bes_ep);
% pre_sigma = sgm(:,:,bes_ep);

end





function [sim_mat1, sim_mat2] = compare_strategy_stability(subtype, sgm, ...
    loglik_all, subtype_all, options)
nsubtype = max(subtype(:));
dat = options.dat;
PTID = options.PTID;

[nsamp, nbioms] = size(dat);

[X, similarity_mat] = reduce_dimension(subtype, sgm);

fold_colors = [1,0,0;0,0,1];

method_colors = [0,1,0;0.5,0,0.5];

max_ep = size(subtype,2);
inds = crossvalind('Kfold', max_ep, 2);

show_fig = 0;

if show_fig
    figure;
end
% compare two strategies on selecting from multiple initializations
for k = 1:2
    inds_k = find(inds == k);
    subtype_k = subtype(:, inds_k);
    loglik_all_k = loglik_all(:, inds_k);
    subtype_all_k = subtype_all(:,:,inds_k);
    sgm_k = sgm(:, :, inds_k);
    
    if show_fig
        scatter(X(inds_k,1), X(inds_k,2), 10, fold_colors(k,:));
        hold on;
    end

     %% filter those converged results 
%     mismatch = double(subtype_all_k(:,2:end,:) ~= subtype_all_k(:,1:end-1,:));
%     mismatch_perc = sum(mismatch, 1) / size(mismatch, 1);
%     mismatch_perc = squeeze(mismatch_perc);
%     filter_inds = find(mismatch_perc(end,:) < 0.01);
%     subtype_k = subtype_k(:, filter_inds);
%     sgm_k = sgm_k(:,:,filter_inds);
    
    %% select the optimal initialization from the common partition
%     bes_ep = sel_based_on_loglik(loglik_k);
%     bes_subtype(:,k) = subtype_k(:,bes_ep);

%     keep_iter = 5;
%     subtype_all_k_keep = subtype_all_k(:,end-keep_iter+1:end,:);
%     subtype_k = subtype_all_k_keep;

    [bes_ep1,sim_mat,bes_subtype,~,~,bes_sgm] = sel_init_by_partition_mode(subtype_k, sgm_k);
    pre_subtype = bes_subtype;
    pre_sigma = bes_sgm;

    options1 = options;
    options1.max_iter = options1.max_iter_refine;
    
    [traj,re_traj,subtype_est,stage,sigma,loglik,proption,extra] = ...
        choose_version(dat,PTID,nsubtype,pre_subtype,pre_sigma,options1);
    
    bes_subtype1(:,k) = subtype_est;
    
    if show_fig
        scatter(X(inds_k(bes_ep1),1), X(inds_k(bes_ep1),2), 9, method_colors(1,:), 'filled');
    end
    
    %% use the similar partitions as initial input
    %     bes_ep2 = sel_based_on_partition_density(subtype_k, loglik_k);
%     bes_subtype2(:,k) = subtype_k(:,bes_ep2);
    
%     if 1
    [subtype_avg, sgm_avg, subtype_perms, sgm_perms] = ...
        sel_based_on_partition_averaging(subtype_k, sgm_k);
    pre_subtype = subtype_perms;
    pre_sigma = sgm_perms;

    options1 = options;

    options1.max_iter = options1.max_iter_refine;

%     pre_subtype = subtype_avg;
%     pre_sigma = sgm_avg;

%     options1.samp_multi = size(pre_subtype, 2);

    [bes_subtype, bes_sgm] = refine_for_similarity(subtype_perms, ...
        sgm_perms, dat, PTID, nsubtype, options1);

    subtype_est = bes_subtype;
    
    bes_subtype2(:,k) = subtype_est;
%     else
%         [bes_ep, bes_subtype] = sel_based_on_neighbor_loglik(subtype_k, loglik_all_k);
%         pre_subtype = bes_subtype;
%         pre_sigma = sgm_k(:,:,bes_ep);
% 
%         options.max_iter = options.max_iter_refine;
% 
%         [traj,re_traj,subtype_est,stage,sigma,loglik,proption,extra] = ...
%             choose_version(dat,PTID,nsubtype,pre_subtype,pre_sigma,options);
% 
%         bes_subtype2(:,k) = subtype_est;
%     end
    
%     scatter(X(inds_k(bes_ep2),1), X(inds_k(bes_ep2),2), 9, method_colors(2,:), 'filled');

%     bes_ep2 = sel_based_on_weighted_loglik(subtype_k, loglik_k);
%     bes_subtype2(:,k) = subtype_k(:, bes_ep2);
% 
%     bes_ep3 = sel_based_on_partition_hierarchical(subtype_k);
%     bes_subtype3(:,k) = subtype_k(:,bes_ep3);
end

% evaluate the 2 strategies using partition similarity
for ep1 = 1:2
    for ep2 = 1:2
        sim_mat1(ep1, ep2) = rand_index(bes_subtype1(:,ep1), bes_subtype1(:,ep2));
        sim_mat2(ep1, ep2) = rand_index(bes_subtype2(:,ep1), bes_subtype2(:,ep2));
%         similarity_mat2(ep1, ep2) = rand_index(bes_subtype2(:,ep1), bes_subtype2(:,ep2));
%         similarity_mat3(ep1, ep2) = rand_index(bes_subtype3(:,ep1), bes_subtype3(:,ep2));
    end
end

if show_fig
sim_mat1
sim_mat2
end

% similarity_mat2
% 
% similarity_mat3

end

function [bes_subtype, bes_sgm] = refine_for_similarity(subtype_perms, ...
    sgm_perms, dat, PTID, nsubtype, options1)

num_tries = 10;
num_folds = 10;

subtype_est_all = zeros(size(dat, 1), num_tries * num_folds);
sigma_est_all = zeros(size(sgm_perms,1), size(sgm_perms,2), num_tries * num_folds);
% sim_mat = zeros(num_tries * num_folds, num_tries * num_folds);

for try_idx = 1:num_tries
    inds = crossvalind('Kfold', size(subtype_perms, 2), num_folds);
    tmp = (try_idx-1)*num_folds;

    options1.samp_multi = 100;
    options1.max_iter = 10;

    for fold_idx = 1:num_folds
        pre_subtype = subtype_perms(:, inds == fold_idx);
        pre_sigma = sgm_perms(:,:, inds == fold_idx);

        pre_subtype = repmat(pre_subtype, [1, options1.samp_multi/size(pre_subtype,2)]);
        pre_sigma = mean(pre_sigma, 3);

        [traj,re_traj,subtype_est,stage,sigma_est,loglik,proption,extra,loglik_all,subtype_all] = ...
            choose_version(dat,PTID,nsubtype,pre_subtype,pre_sigma,options1);

%         check_convergence(loglik_all,subtype_all);

        subtype_est_all(:, tmp + fold_idx) = subtype_est;
        sigma_est_all(:,:, tmp + fold_idx) = sigma_est;
    end

    subtype_est_all1 = subtype_est_all(:,1:try_idx*num_folds);
    sigma_est_all1 = sigma_est_all(:,:,1:try_idx*num_folds);

    [bes_ep1, similarity_mat, bes_subtype, eps, subtype, bes_sgm] = ...
        sel_init_by_partition_mode(subtype_est_all1, sigma_est_all1, options1);

    sim_mat = similarity_mat;
    sim_mat = double(sim_mat > 0.98);
    if any(sum(sim_mat, 2) >= 10)
        disp(['break at ', num2str(try_idx)]);
        break
    end
end

end



function [bes_subtype, bes_sgm, subtype_perms, sgm_perms] = sel_based_on_partition_averaging(subtype, sgm)
[bes_ep1, similarity_mat,bes_subtype,eps,subtype] = sel_init_by_partition_mode(subtype, sgm);
sim1 = sort(similarity_mat, 2, 'descend');
sim_vec = sum(sim1(:,2:end), 2);

[vec_sorted, inds] = sort(sim_vec, 'descend');

% On OASIS, for 200 runs, split the runs into 100/100 to validate the
% reproducibility. Over 30 splitting, the competing strategy (sel_based_on_partition) has 
% median ARI 0.79. This strategy (sel_based_on_partition_averaging) has median 0.93 
% for percentage 0.3, median 0.97 for 0.5, median 0.96 for 1
keep_top_perc = 1;

inds_selected = inds(1:round(size(subtype,2)*keep_top_perc));

subtype_sel = subtype(:, inds_selected);
sgm_sel = sgm(:, :, eps(inds_selected));
[bes_subtype, bes_sgm, subtype_perms, sgm_perms] = average_subtypes(subtype_sel, sgm_sel);

end



function [X, similarity_mat] = reduce_dimension(subtype, sgm)
[bes_ep1, similarity_mat] = sel_init_by_partition_mode(subtype, sgm);
W = similarity_mat;
W = W - eye(size(W,1));
D = diag(sum(W,2));
L = D - W;
[V,D] = eig(L);
X = V(:,2:3);
end

function visualize_by_spectral_clustering(subtype, loglik_all)
[X, similarity_mat] = reduce_dimension(subtype);

% figure, scatter(X(:,1),X(:,2));
% 
% sz = loglik_all(end,:);
% sz = 1./(-sz);
% sz = sz / max(sz);
% sz = sz * 5;
figure, scatter(X(:,1),X(:,2),5);
hold on;
scatter(X(bes_ep1,1),X(bes_ep1,2),5,'r');
end

%% obsolete. Not stable compared to pick the one with the highest average similarity
function best_ep = sel_based_on_partition_density(subtype, loglik_all)
[bes_ep1, similarity_mat] = sel_init_by_partition_mode(subtype);
sim1 = sort(similarity_mat, 2, 'descend');
sim2 = sim1(:,2:11);
sim_vec = sum(sim2, 2);
[vec_max, best_ep] = max(sim_vec);
end

function best_ep = sel_based_on_partition_gmm(subtype, loglik_all)
[bes_ep1, similarity_mat] = sel_init_by_partition_mode(subtype);
sim1 = sort(similarity_mat, 2, 'descend');
sim2 = sim1(:,2:end);
sim3 = sum(sim2, 2);
gmm = fitgmdist(sim3, 2, 'Replicates', 100);
idx = cluster(gmm, sim3);
[max_mu, max_idx] = max(gmm.mu);
sim_gmm = similarity_mat(idx == max_idx, idx == max_idx);
sum_vec = sum(sim_gmm, 2);
[vec_max, idx_vec] = max(sum_vec);
inds = find(idx == max_idx);
best_ep = inds(idx_vec);
end

function bes_ep2 = sel_based_on_weighted_loglik(subtype, loglik_all)
max_ep = size(subtype, 2);
similarity_mat = zeros(max_ep, max_ep);
for ep1 = 1:max_ep
    for ep2 = 1:max_ep
        if any(isnan(subtype(:,ep1))) || any(isnan(subtype(:,ep2)))
            similarity_mat(ep1, ep2) = 0;
        else
            similarity_mat(ep1, ep2) = rand_index(subtype(:,ep1), subtype(:,ep2));
        end
    end
end

loglik_mat = loglik_all(end,:);
loglik = repmat(loglik_mat, [max_ep, 1]);

% weight = similarity_mat ./ repmat(sum(similarity_mat, 2), 1, max_ep);
weight = similarity_mat;

sim_vec = sum(loglik .* weight, 2, 'omitnan');
[~, bes_ep2] = max(sim_vec);
end

function [bes_ep, bes_subtype] = sel_based_on_neighbor_loglik(subtype, loglik_all)
max_ep = size(subtype, 2);
similarity_mat = zeros(max_ep, max_ep);
for ep1 = 1:max_ep
    for ep2 = 1:max_ep
        if any(isnan(subtype(:,ep1))) || any(isnan(subtype(:,ep2)))
            similarity_mat(ep1, ep2) = 0;
        else
            similarity_mat(ep1, ep2) = rand_index(subtype(:,ep1), subtype(:,ep2));
        end
    end
end

neighbor_sz = round(max_ep * 0.5);

[sim_sorted, inds] = sort(similarity_mat, 2, 'descend');

loglik_mat = mean(loglik_all(end-4:end,:), 1);

loglik_neighbor_avg = loglik_mat;

for i = 1:size(inds,1)
    inds_i = inds(i,1:neighbor_sz);
    weight_i = similarity_mat(i, inds_i);
    weight_i = weight_i / sum(weight_i);
    loglik_i = loglik_mat(inds_i);
    loglik_neighbor_avg(i) = mean(loglik_i .* weight_i, 2, 'omitnan');
end

[~, bes_ep] = max(loglik_neighbor_avg);
bes_subtype = subtype(:, bes_ep);

end

function bes_ep = sel_based_on_loglik(loglik_all)
loglik_mat = loglik_all(end,:);
[~,bes_ep] = max(loglik_mat);
end


function bes_ep3 = sel_based_on_partition_hierarchical(subtype)
max_ep = size(subtype, 2);
num_repeats = floor(sqrt(max_ep));
inds = crossvalind('Kfold', max_ep, num_repeats);

bes_eps = [];
for k = 1:num_repeats
    inds_k = find(inds == k);
    subtype_k = subtype(:, inds_k);
    [bes_ep1, similarity_mat] = sel_init_by_partition_mode(subtype_k);
    bes_eps(k) = inds_k(bes_ep1);
end

bes_ep2 = sel_init_by_partition_mode(subtype(:, bes_eps));
bes_ep3 = bes_eps(bes_ep2);

end
