clear

demo_num = 14;
cov_cols = 5:8;

hemisphere_status = {'HM','HS','LM','LS'};
num_volume_regions = [41,82,13,26];
group = {'ROD1'};

for i = 1:length(hemisphere_status)
    for j = 1:length(group)
        read_dir = ['./Result/OASIS3_input_', group{j}, '_', hemisphere_status{i}, '.csv'];
        save_dir = ['./Result/OASIS3_input_', group{j}, '_', hemisphere_status{i}, '_z.csv'];
        regress_z_table(read_dir, save_dir, demo_num, cov_cols, num_volume_regions(i), hemisphere_status{i});
    end
end

function regress_z_table(read_dir, save_dir, demo_num, cov_cols, num_volume, hemisphere_status)

input = readtable(read_dir,'VariableNamingRule',"preserve");
% input(input.FLDSTRENG == 1.5, :) = [];
mat = input{:,demo_num+1:end};
nbiom = size(mat,2);
mat_all_ro = regress_cov(input,mat,demo_num,cov_cols,num_volume);

%% get cortical & hippocampus

mapping = readtable('./FSX_biomarker_mapping.csv');

switch hemisphere_status
    case {'HM','HS'}
        CTV_column = mapping{strcmp(mapping.HV_CTV,'Cortical'), hemisphere_status};
        CTV_column = unique(cellfun(@(s) strrep(s, 'FSX-', ''), CTV_column, 'UniformOutput', false));
        HV_column = mapping{strcmp(mapping.HV_CTV,'Hippocampus'), hemisphere_status};
        HV_column = unique(cellfun(@(s) strrep(s, 'FSX-', ''), HV_column, 'UniformOutput', false));
        CTV = sum(input{:,CTV_column},2);
        HV = sum(input{:,HV_column},2);
        ctv_hv_need_regress = 1;
    case {'LM','LS'}
        % No Cortical data in LM and LS version.
        % Read previously processed data from HM and HS version.
        previous_result = readtable(strrep(save_dir, hemisphere_status, strrep(hemisphere_status, 'L', 'H')), 'VariableNamingRule','preserve');
        CTV = previous_result.CTV;
        HV = previous_result.HV;
        ctv_hv_need_regress = 0; % Thus no need to regress
    otherwise
end

if ctv_hv_need_regress
    ctv_hv = regress_cov(input,[CTV,HV],demo_num,cov_cols,2);
    CTV = ctv_hv(:,1);
    HV = ctv_hv(:,2);
end

%% zscore

% zscore only for mri
mat_z_score = z_score(input,mat_all_ro(:, 1:num_volume),demo_num,cov_cols,num_volume);


%% save table
for i = 1:num_volume
    input{:,i+demo_num} = mat_z_score(:,i);
end

input = input(:,[1,2,demo_num+1:end,3:demo_num]);

% add CTV & HV
input.CTV = CTV;
input.HV = HV;

writetable(input,save_dir)
end


function mat_all_ro = regress_cov(input,mat,demo_num,cov_cols,num_volume)

years = input.years;
diagnosis = input.diagnosis;
flag = input.group;
mat_icv_ro = zeros(size(mat));
mat_all_ro = zeros(size(mat));
mat_z_score = zeros(size(mat));
nbiom = size(mat,2);

ICV = input.ICV;
fldstr = input.FLDSTRENG;
cov = input{:,cov_cols};

idx_CN = flag==0;
cov_CN = input{idx_CN,cov_cols};
icv_CN = input{idx_CN,'ICV'};
mat_CN = input{idx_CN,demo_num+1:end};
one = ones(sum(idx_CN),1);

mean_icv_CN = mean(icv_CN);
mean_cov_CN = mean(cov_CN,1);

%% regress out ICV
b_icv = zeros(nbiom,2);
p_icv = zeros(nbiom,1);
for i = 1:nbiom
    if i <= num_volume
        [b_icv(i,:),~,~,~,stats] = regress(mat_CN(:,i),[icv_CN,one]); 
        p_icv(i) = stats(3);
        mat_icv_ro(:,i) = mat(:,i) - b_icv(i,1) * (ICV - mean_icv_CN);
    else
        mat_icv_ro(:,i) = mat(:,i);
    end
end

mat_all_ro(:,1:num_volume) = regress_out_if_significant(mat_icv_ro(:,1:num_volume), idx_CN, cov);
% mat_all_ro(:,1:num_volume) = regress_out_if_significant(mat_all_ro(:,1:num_volume), idx_CN, fldstr);

% %% covariant p-value
% mat_icv_ro_CN = mat_icv_ro(idx_CN,:);
% b_cov = zeros(nbiom,4);
% p_cov = zeros(nbiom,4);
% for i = 1:nbiom
%     for j = 1:4
%         [b,~,~,~,stats] = regress(mat_icv_ro_CN(:,i),[cov_CN(:,j),one]);
%         b_cov(i,j) = b(1);
%         p_cov(i,j) = stats(3);
%     end
% end
% 
% %% regress out age/sex/edu/APOE
% ro_tf = p_cov < 0.05 * ones(size(p_cov));
% p_all = zeros(nbiom,1);
% for i = 1:nbiom
%     [b,~,~,~,stats] = regress(mat_icv_ro_CN(:,i),[cov_CN(:,ro_tf(i,:)),one]);
%     p_all(i) = stats(3);
%     if sum(ro_tf(i,:)) == 0
%         mat_all_ro(:,i) = mat_icv_ro(:,i);
%     else
%         mat_all_ro(:,i) = mat_icv_ro(:,i) - (cov(:,ro_tf(i,:))-mean_cov_CN(:,ro_tf(i,:))) * b(1:end-1);
%     end
% end

end

function mat_z_score = z_score(input,mat_all_ro,demo_num,cov_cols,num_volume)

years = input.years;
diagnosis = input.diagnosis;
flag = input.group;
mat_z_score = zeros(size(mat_all_ro));
nbiom = size(mat_all_ro,2);

ICV = input.ICV;
cov = input{:,cov_cols};

idx_CN = flag==0;
cov_CN = input{idx_CN,cov_cols};
icv_CN = input{idx_CN,'ICV'};
mat_CN = input{idx_CN,demo_num+1:end};
one = ones(sum(idx_CN),1);

mean_icv_CN = mean(icv_CN);
mean_cov_CN = mean(cov_CN,1);

%% normalize z-score
for i = 1:nbiom
    pd = fitdist(mat_all_ro(idx_CN,i),'Normal');
    if i <= num_volume
        mat_z_score(:,i) = (mat_all_ro(:,i)-pd.mu)/pd.sigma * -1;
    else
        mat_z_score(:,i) = (mat_all_ro(:,i)-pd.mu)/pd.sigma;
    end
end

end

function mat_all_ro = regress_out_if_significant(bm_all_icv_ro, idx_CN, cov)
mat_icv_ro_CN = bm_all_icv_ro(idx_CN, :);
nbiom = size(mat_icv_ro_CN, 2);
num_cov = size(cov, 2);

cov_CN = cov(idx_CN, :);

one = ones(size(cov_CN, 1), 1);
mean_cov_CN = nanmean(cov_CN,1);

mat_all_ro = zeros(size(bm_all_icv_ro));

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
        mat_all_ro(:,i) = bm_all_icv_ro(:,i);
    else
        mat_all_ro(:,i) = bm_all_icv_ro(:,i) - (cov(:,ro_tf(i,:))-mean_cov_CN(:,ro_tf(i,:))) * b(1:end-1);
    end
end
end