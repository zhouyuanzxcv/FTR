function subtype = init_subtype(nsubj, nsubtype, dat, PTID, options)
%INIT_SUBTYPE Summary of this function goes here
%   Detailed explanation goes here
% [nsamp, nbiom] = size(dat);
% 
% [unique_PTIDs,ia,ic] = unique(PTID);
% [yy, xx] = ndgrid(ic, 1:nbiom);
% dat_grouped = accumarray([yy(:),xx(:)], dat(:), [], @mean);
% data1 = dat_grouped ./ repmat(sqrt(sum(dat_grouped.^2, 2)), [1, nbiom]);

if nsubtype == 1
    subtype = ones(size(dat,1), 1);
    
    return;
end


Ts = calc_Ts(PTID);

init_method = parse_param(options, 'init_method', 'random_permutation');

if strcmp(init_method, 'randi')
    subtype = randi(nsubtype, [nsubj, 1]);
    subtype = repelem(subtype, Ts, 1);
elseif strcmp(init_method, 'random_permutation')
    subtype = repmat((1:nsubtype)', [ceil(nsubj/nsubtype), 1]);
    inds = randperm(nsubj);
    subtype = subtype(inds);

    subtype = repelem(subtype, Ts, 1);
elseif strcmp(init_method, 'kmeans++')
    [nsamp, nbiom] = size(dat);
    start_pt = zeros(1, nbiom);
    end_pt = max(dat, [], 1);
    
    % generate centroids similar to kmeans++
    num_int = 101;
    re_traj = zeros(nbiom, num_int, nsubtype);

    sample_idx = randi(nsamp);

    D1s = [];
    traj_pts = [];

    for k = 1:nsubtype
        traj_pts = cat(1, traj_pts, dat(sample_idx, :));
        re_traj(:,:,k) = reparam_traj([start_pt; traj_pts(k,:); end_pt], num_int)';

        % square the distance to amplify the probability of distant points
        D = pdist2(re_traj(:,:,k)', dat, 'euclidean');
        D1 = min(D, [], 1);
        D1s = cat(1, D1s, D1);
        D1 = min(D1s, [], 1);
        D2 = D1.^2;
        D2_p = D2 ./ sum(D2);
        sample_idx = randsample(1:nsamp, 1, true, D2_p);
    end

    % assign subtypes
    D1s = D1s';
    [unique_PTIDs,ia,ic] = unique(PTID);
    [yy, xx] = ndgrid(ic, 1:nsubtype);
    D1s_grouped = accumarray([yy(:),xx(:)], D1s(:));

    [~, subtype] = min(D1s_grouped, [], 2);
    subtype = subtype(ic);
    
%     [stage,subtype,dist_samp_min] = cal_stage_subtype_kmeans(dat,PTID,re_traj,Ts);

    % figure; show_trajectory_in_PCA(dat, subtype, subtype, re_traj);
end


% centers = [];
% for k = 1:nsubtype
%     centers(k,:) = mean(data1(subtype == k, :), 1);
% end
% 
% subtype = kmeans(data1, nsubtype, 'Start', centers);

end

