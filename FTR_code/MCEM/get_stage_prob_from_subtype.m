function stage_prob = get_stage_prob_from_subtype(loglik_stage_nor, subtype_samp)
%GET_STAGE_PROB_FROM_SUBTYPE Summary of this function goes here
%   Detailed explanation goes here
% Input:
%   loglik_stage_nor - nsamp x num_int x K tensor where the element at the ith row, kth page
%                       contains the stage probabilities from the ith data point on the kth
%                       trajectory
%   subtype_samp - L x nsamp matrix where the lth row contains the subtype
%                   labels of all the data points from the lth sample of z.
% Output:
%   stage_prob - nsamp x num_int x L matrix where the lth page contains the 
%                stage probabilities from the lth subtype labels

[nsamp, num_int, nsubtype] = size(loglik_stage_nor);
[samp_multi, nsamp] = size(subtype_samp);

stage_prob = zeros(nsamp, num_int, samp_multi);
for l = 1:samp_multi
    subtype_l = subtype_samp(l,:)';
    subtype_l = repmat(subtype_l, 1, num_int);
    stage_l = repmat((1:num_int), nsamp, 1);
    dat_l = repmat((1:nsamp)', 1, num_int);
    idx_l = sub2ind([nsamp, num_int, nsubtype], dat_l(:), stage_l(:), subtype_l(:));
    stage_prob(:,:,l) = reshape(loglik_stage_nor(idx_l), [nsamp, num_int]);
end

end

