function [ind_reject, ps_corrected] = multi_test_correction(ps, level, mode)
%MULTI_TEST_CORRECTION Summary of this function goes here
%   Detailed explanation goes here
m = length(ps);
if strcmp(mode, 'bonferroni')
    ps_corrected = ps * m;
    ind_reject = ps_corrected < level;
elseif strcmp(mode, 'BH') % Benjamini-Hochberg correction
    [ps_sorted, inds] = sort(ps);
    ps_sorted_corrected = ps_sorted(:) * m ./ (1:m)';
    ps_corrected = nan(size(ps));
    ps_corrected(inds) = ps_sorted_corrected;
    ind_reject = ps_corrected < level;

%     compare = ps_sorted(:) <= (1:m)'/m*level;
%     k = find(~compare, 1);
%     if isempty(k)
%         ind_reject = true(size(ps));
%     elseif k == 1
%         ind_reject = false(size(ps));
%     else
%         k = k-1;
%         ind_reject = false(size(ps));
%         ind_reject(inds(1:k)) = true;
%     end
end

end

