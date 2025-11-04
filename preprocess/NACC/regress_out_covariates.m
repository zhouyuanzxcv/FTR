function input = regress_out_covariates(input, idx_CN, biom_column, mri_column, cov_column)
num_cov = length(cov_column);

% diagnosis = input.diagnosis;
bm_all = input(:,biom_column);
bm_all_icv_ro = bm_all;

mat_all_ro = zeros(size(bm_all));
nbiom = size(bm_all,2);

cov = input{:,cov_column};
% abeta = input.ABETA;
ICV = input.ICV;

%idx_CN = cellfun(@(x) strcmp(x, 'CN'), input.diagnosis) & abeta>880;
% idx_CN = diagnosis == 0 & abeta>880;

cov_CN = input{idx_CN,cov_column};
% mat_CN = input{idx_CN,biom_column};

icv_CN = ICV(idx_CN);
mean_cov_CN = nanmean(cov_CN,1);
mean_icv_CN = nanmean(icv_CN,1);
one = ones(sum(idx_CN),1);

%% regress out ICV for MRI
b_icv = zeros(nbiom,2);
p_icv = zeros(nbiom,1);

for i = 1:length(mri_column)
    bm = mri_column{i};
    % for each MRI biomarker, the slope is calculated from the CN group and
    % propagated to the entire population
    [b_icv(i,:),~,~,~,stats] = regress(input{idx_CN, bm},[icv_CN,one]);
    p_icv(i) = stats(3);
    bm_all_icv_ro{:,bm} = bm_all{:,bm} - b_icv(i,1) * (ICV - mean_icv_CN);
end

mat_icv_ro_CN = bm_all_icv_ro{idx_CN,:};

%% determine if a covariate is significant
b_cov = zeros(nbiom,num_cov);
p_cov = zeros(nbiom,num_cov);
for i = 1:nbiom
    for j = 1:num_cov
        [b,~,~,~,stats] = regress(mat_icv_ro_CN(:,i),[cov_CN(:,j),one]);
        b_cov(i,j) = b(1);
        p_cov(i,j) = stats(3);
    end
end

%% regress out covariates
ro_tf = p_cov < 0.05 * ones(size(p_cov));
p_all = zeros(nbiom,1);
for i = 1:nbiom
    [b,~,~,~,stats] = regress(mat_icv_ro_CN(:,i),[cov_CN(:,ro_tf(i,:)),one]);
    p_all(i) = stats(3);
    % if none of the covariates is significant
    if sum(ro_tf(i,:)) == 0
        mat_all_ro(:,i) = bm_all_icv_ro{:,i};
    else
        mat_all_ro(:,i) = bm_all_icv_ro{:,i} - (cov(:,ro_tf(i,:))-mean_cov_CN(:,ro_tf(i,:))) * b(1:end-1);
    end
end

input{:, biom_column} = mat_all_ro;

end