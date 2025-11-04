function regress_normalize(freesurfer_ver, biomarker_ver)

input = readtable(['./Result/ADNI_demo_MRI_',freesurfer_ver,'_PET_', ...
    biomarker_ver,'_all.csv'],VariableNamingRule='preserve');
data = input;

%% Verify the csf abeta cutoff

csf_abeta_cutoff = get_csf_cutoff(data);

%% AD and CN criteria
% Iterate through each row of the table
group = NaN(size(data,1), 1);

unique_RIDs = unique(data.RID);

for curr_rid = unique_RIDs'
    inds = find(data.RID == curr_rid);
    
    curr_diagnosis = data.diagnosis(inds);
    curr_csf_abeta = data.ABETA(inds);
    curr_pet_abeta = data.("Abeta-positive")(inds);
    
    cond1 = any(curr_diagnosis == 0.5) || any(curr_diagnosis == 1);
    cond2 = any(curr_csf_abeta < csf_abeta_cutoff) || any(curr_pet_abeta == 1);
    
    % a subject belongs to the case group if at some time (i) the
    % diagnosis is MCI or AD and (ii) the CSF abeta is positive or the PET
    % abeta is positive
    if cond1 && cond2
        group(inds) = 1;
    end

    % a subject belongs to the validation group if at some time (i) the
    % diagnosis is AD and (ii) the CSF abeta and PET abeta is missing
    cond5 = any(curr_diagnosis == 1);
    cond6 = all(isnan(curr_csf_abeta)) && all(isnan(curr_pet_abeta));
    if cond5 && cond6
        group(inds) = 2;
    end
    
    % a subject belongs to the control group if at all times (i) the
    % diagnosis is CN and (ii) the CSF abeta is negative.
    curr_diagnosis1 = curr_diagnosis;
    curr_diagnosis1(isnan(curr_diagnosis1)) = [];

    curr_csf_abeta1 = curr_csf_abeta;
    curr_csf_abeta1(isnan(curr_csf_abeta1)) = [];

    cond3 = ~isempty(curr_diagnosis1) && all(curr_diagnosis1 == 0);
    cond4 = ~isempty(curr_csf_abeta1) && all(curr_csf_abeta1 > csf_abeta_cutoff);

    if cond3 && cond4
        group(inds) = 0;
    end
    
end

data = addvars(data,group,'After','diagnosis');

%% regress out age, gender, APOE4, education and then normalize to z-score
cov_column = ["AGE","PTGENDER","APOE4","PTEDUCAT"];

cov_nans = any(isnan(data{:,cov_column}), 2);
data(cov_nans, :) = [];


[all_names, mri_bio_names, abeta_bio_names, tau_bio_names] = ...
    get_biomarker_names(freesurfer_ver, biomarker_ver);

% only perform covariate regression and zscore normalization for MRI
% measures
biom_column = setdiff(all_names, union(abeta_bio_names, tau_bio_names));
mri_column = mri_bio_names;


% % impute missing scanner strength with freesurfer version
idx_fldstr_missing = isnan(data.FLDSTRENG) & ~isnan(data.FSVERSION);
fsver2fldstr = containers.Map([4.3, 5.1, 6], [1.5, 3, 3]);
fldstr_values = arrayfun(@(x) fsver2fldstr(x), data{idx_fldstr_missing, 'FSVERSION'});
data{idx_fldstr_missing, 'FLDSTRENG'} = fldstr_values;

% remove missing scanner strength afterwards
data(isnan(data.FLDSTRENG), :) = [];

data_all = data;
data_types = [0, 1.5, 3];

for data_type = data_types
    if data_type ~= 0
        data = data_all(data_all.FLDSTRENG == data_type, :);
    else
        data = data_all;
    end

    % regress and normalize
    idx_CN = data.group == 0;
    biom_add = {'Cortical','Hippocampus'};
    data = regress_out_covariates(data, idx_CN, [biom_column,biom_add], ...
        [mri_column,biom_add], cov_column, data_type);

    data = normalize_to_zscore(data, idx_CN, biom_column, mri_column);

    %% rename biomarkers

    curr_col_names = data.Properties.VariableNames;

    curr_col_names{strcmp(curr_col_names, 'Hippocampus')} = 'HV';
    curr_col_names{strcmp(curr_col_names, 'Cortical')} = 'CTV';

    if strcmp(freesurfer_ver, 'FSX')
        curr_col_names = strrep(curr_col_names, 'FSX-', '');
    elseif strcmp(freesurfer_ver, 'FSL')
        curr_col_names = strrep(curr_col_names, 'FSL-', '');
    end

    data.Properties.VariableNames = curr_col_names;

    if data_type == 0
        writetable(data, ['./Result/ADNI_demo_MRI_',freesurfer_ver,'_PET_',biomarker_ver,'_zscore.csv']);
    else
        writetable(data, ['./Result/ADNI_demo_MRI_',num2str(data_type),'_',freesurfer_ver,'_PET_',biomarker_ver,'_zscore.csv']);
    end

end

end

function input1 = normalize_to_zscore(input, idx_CN, biom_column, mri_column)
input1 = input;
nbiom = length(biom_column);

% normalize to z-score for all biomarkers
sigma = zeros(nbiom,1);
for i = 1:length(biom_column)
    bm = biom_column{i};
    pd = fitdist(input{idx_CN, bm},'Normal');
    input1{:, bm} = (input{:, bm} - pd.mu)/pd.sigma;
    sigma(i) = pd.sigma;
end

% add a minus sign for MRI
for bm = mri_column
    input1{:, bm} = -input1{:, bm};
end

end

function input = regress_out_covariates(input, idx_CN, biom_column, mri_column, cov_column, data_type)
% if data_type == 0
%     cov_column = cat(2, cov_column, 'FLDSTRENG');
% end

% num_cov = length(cov_column);

% diagnosis = input.diagnosis;
bm_all = input(:,biom_column);
bm_all_icv_ro = bm_all;

% mat_all_ro = zeros(size(bm_all));
% nbiom = size(bm_all,2);

cov = input{:,cov_column};
% abeta = input.ABETA;
ICV = input.ICV;
fldstr = input.FLDSTRENG;

%idx_CN = cellfun(@(x) strcmp(x, 'CN'), input.diagnosis) & abeta>880;
% idx_CN = diagnosis == 0 & abeta>880;

% cov_CN = input{idx_CN,cov_column};
% mat_CN = input{idx_CN,biom_column};

icv_CN = ICV(idx_CN);
% fldstr_CN = fldstr(idx_CN);

% mean_cov_CN = nanmean(cov_CN,1);
mean_icv_CN = nanmean(icv_CN,1);
% mean_fldstr_CN = nanmean(fldstr_CN, 1);

one = ones(sum(idx_CN),1);

%% regress out ICV and scanner strength for MRI
b_icv = [];
p_icv = [];

for i = 1:length(mri_column)
    bm = mri_column{i};
    % for each MRI biomarker, the slope is calculated from the CN group and
    % propagated to the entire population

    [b_icv(i,:),~,~,~,stats] = regress(input{idx_CN, bm}, [icv_CN, one]);
    p_icv(i,:) = stats(3);

    bm_all_icv_ro{:,bm} = bm_all{:,bm} - (ICV - mean_icv_CN) * b_icv(i,1)';

end


mat_all_ro = bm_all_icv_ro{:,:};
mat_all_ro = regress_out_if_significant(mat_all_ro, idx_CN, cov);

if data_type == 0
    mat_all_ro = regress_out_if_significant(mat_all_ro, idx_CN, fldstr);
end

input{:, biom_column} = mat_all_ro;

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