function [stage,subtype,extra] = cal_stage_subtype(dat, PTID, mdl, options)

re_traj = mdl.re_traj;
proption = mdl.proption;
sgm = mdl.sigma;

method = parse_param(options, 'methods', 'kmeans');

% unique_PTIDs = unique(PTID);
% Ts = zeros(1, length(unique_PTIDs));
% for i = 1:length(Ts)
%     Ts(i) = length(find(unique_PTIDs(i) == PTID));
% end

Ts = calc_Ts(PTID);

if strcmp(method, 'kmeans')
    [stage,subtype,dist_samp_min] = cal_stage_subtype_kmeans(dat, PTID, re_traj, Ts);
    extra = [];
elseif strcmp(method, 'MCEM')
    [stage,subtype,extra] = cal_stage_subtype_MCEM(dat, PTID, re_traj, proption, sgm, Ts);
end

end


