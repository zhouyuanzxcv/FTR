function [ord, reach_stages] = calc_order(traj, stages, thresh)
%CALC_ORDER Summary of this function goes here
%   Detailed explanation goes here
% Input:
%   traj - B x num_int matrix

comp = sum(traj < thresh,2);
[~,ord] = sortrows([comp, mean(traj,2)]);

[B,num_int] = size(traj);

reach_stages = [];
for j = 1:B
    idx = find(traj(j,:) >= thresh, 1);
    if isempty(idx)
        reach_stages(j) = stages(end) + (stages(end) - stages(end-1));
    else
        reach_stages(j) = stages(idx);
    end
end

reach_stages = reach_stages';

[~,ord1] = sortrows([reach_stages, mean(traj,2)]);
assert(all(ord1 == ord));

end

