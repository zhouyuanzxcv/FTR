function [stage,subtype,dist_samp_min] = cal_stage_subtype_kmeans(dat,PTID,re_traj,Ts)
[nsamp,nbiom] = size(dat);
nsubtype = size(re_traj,3);
dist_samp = zeros(nsamp,nbiom,nsubtype);
diff_samp = zeros(nsamp,nbiom,nsubtype);
stag_samp = zeros(nsamp,nsubtype);

% dist_samp is the squared distance of each biomarker
for k = 1:nsubtype
    [stag_samp(:,k),dist_samp(:,:,k),diff_samp(:,:,k)] = cal_dis(dat,re_traj(:,:,k));
end

% for each patient, an average of the squared distances of all the points
% is assigned to each point such that a patient has a consistent subtype
dist_samp1 = dist_samp;

if 0
    uni_id = unique(PTID);
    for i = 1:length(uni_id)
        ind_PTID_i = PTID==uni_id(i);
        rws = dist_samp(ind_PTID_i,:,:);
        dist_samp1(ind_PTID_i,:,:) = repmat(mean(rws,1),sum(ind_PTID_i),1,1);
    end
else
    % change the above code to a faster implementation
    group_ids = (1:length(Ts))';
    group_ids = repelem(group_ids, Ts, 1);
    [yy, xx, zz] = ndgrid(group_ids, 1:nbiom, 1:nsubtype);
    dist_samp1 = accumarray([yy(:),xx(:),zz(:)], dist_samp(:), [], @mean);
    dist_samp1 = repelem(dist_samp1, Ts, 1, 1);    
end

dist_samp = dist_samp1;

% find the subtype of each point by summing up the squared distances of all
% the biomarkers
[~,subtype] = min(sum(dist_samp,2),[],3);

p = 1:nsamp;
ind1 = sub2ind(size(stag_samp),p',subtype);
stage = stag_samp(ind1);


mat1 = repmat(1:nsamp,nbiom,1).';
mat2 = repmat(1:nbiom,nsamp,1);
mat3 = repmat(subtype,1,nbiom);
ind2 = sub2ind(size(dist_samp),mat1(:),mat2(:),mat3(:));
dist_samp_min = reshape(dist_samp(ind2),nsamp,nbiom);

% set nan values
nan_inds = any(isnan(dat), 2);
stage(nan_inds) = NaN;
subtype(nan_inds) = NaN;
dist_samp_min(nan_inds, :) = NaN; 

end