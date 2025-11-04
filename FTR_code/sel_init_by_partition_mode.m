function [bes_ep1, similarity_mat, bes_subtype, eps, subtype, bes_sgm, metric] = ...
    sel_init_by_partition_mode(subtype, sgm, options)
if nargin < 3
    options = [];
end

thresh_same_partition1 = parse_param(options, 'thresh_same_partition', 0.98);

if ndims(subtype) == 2
    [nsamp, n_ep] = size(subtype);
    eps = (1:n_ep);
elseif ndims(subtype) == 3
    [nsamp, niters, n_ep] = size(subtype);
    subtype = reshape(subtype, [nsamp, niters*n_ep]);
    eps = repelem((1:n_ep), 1, niters);
end

max_ep = size(subtype, 2);
similarity_mat = zeros(max_ep, max_ep);
for ep1 = 1:max_ep
    for ep2 = 1:max_ep
        if any(isnan(subtype(:,ep1))) || any(isnan(subtype(:,ep2)))
            similarity_mat(ep1, ep2) = NaN;
        else
            similarity_mat(ep1, ep2) = rand_index(subtype(:,ep1), subtype(:,ep2));
        end
    end
end

thresh_same_partition = thresh_same_partition1 * max(similarity_mat(:));

sim1 = sort(similarity_mat, 2, 'descend');
sim2 = sim1(:,1:end);
sim2(sim2 <= thresh_same_partition) = 0;

sim_vec = sum(sim2, 2, 'omitnan');
[metric, bes_ep1] = max(sim_vec);
metric = sum(double(sim2(bes_ep1,:) > thresh_same_partition)); % number of same results

bes_subtype = subtype(:,bes_ep1);

% for ndims(subtype) == 3
bes_ep1 = eps(bes_ep1);

bes_sgm = sgm(:,:,bes_ep1);

% figure, histogram(similarity_mat(:));
% figure,plot(loglik_all);

end
