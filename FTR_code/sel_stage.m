function [dat_sel,PTID_sel,stage_sel,inds_sel] = sel_stage(dat,PTID,options)

nsubtype = 1;

[mdl,subtype,stage,extra] = FTR_model(dat,PTID,nsubtype,options);

traj = mdl.re_traj;
[num_bioms,num_int] = size(traj);

% remove all points whose modeled point is the origin, i.e. f(s_i) = 0
% this is equivalent to removing the data points with stages being 0
traj = [zeros(num_bioms,1), traj];
dist = sqrt(sum((traj(:,2:end) - traj(:,1:end-1)).^2, 1));

dist2origin = cumsum(dist);
k = find(dist2origin, 1);

stage_line = linspace(0, 1, num_int);

if k > 1
    thresh = stage_line(k-1);
else
    % keep all the stages since the first point is already not at the
    % origin
    thresh = -Inf;
end

inds_sel = stage>thresh;
dat_sel = dat(inds_sel,:);
PTID_sel = PTID(inds_sel,:);
stage_sel = stage(inds_sel);

end
