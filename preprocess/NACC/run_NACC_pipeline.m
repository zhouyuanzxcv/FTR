function [outputArg1,outputArg2] = run_NACC_pipeline(inputArg1,inputArg2)
%ADNI_PIPELINE Summary of this function goes here
%   Detailed explanation goes here

prepare();

population_selection();

biomarker_vers = {'HM','HS','LM','LS'};

mri = load_mri();

for j = 1:length(biomarker_vers)
    mri_ver = feature_selection(mri, biomarker_vers{j});
    regress_zcore(mri_ver, biomarker_vers{j});
end

copy_files_to_destination_NACC(biomarker_vers);

end

function prepare()
clinical_path = ['./data/investigator_ftldlbd_nacc69.csv'];
clinical_file = readtable(clinical_path, VariableNamingRule='preserve');

kept_columns = {'NACCID','NACCADC','FORMVER','NACCALZD','NACCUDSD','CDRGLOB',...
    'NACCMMSE','NACCMOCA','CDRSUM','HACHIN','NACCGDS','PACKET','NACCVNUM','NACCAVST',...
    'NACCNVST','NACCFDYS','VISITMO','VISITDAY','VISITYR',...
    'BIRTHMO','BIRTHYR','SEX','RACE','EDUC','NACCAGE','NACCAGEB','NACCBMI',...
    'NACCNE4S'};

inds = clinical_file.NACCALZD == 8 | clinical_file.NACCALZD == 1;
inds1 = ismember(clinical_file.NACCID, unique(clinical_file.NACCID(inds)));

file_reduced = clinical_file(inds1, kept_columns);
writetable(file_reduced, './data/nacc_clinical.csv');

end

function population_selection()

%% remove vascular dementia (HACHIN > 4) and depression (NACCGDS > 5)
clinical = readtable('./data/nacc_clinical.csv', VariableNamingRule='preserve');

inds = clinical.HACHIN > 4 | clinical.NACCGDS > 5;
vascular_dep_ids = unique(clinical{inds, 'NACCID'});
inds_remove = ismember(clinical.NACCID, vascular_dep_ids);
clinical(inds_remove, :) = [];


% remove other dementia
inds_remove = ismember(clinical.NACCID, unique(clinical{clinical.NACCALZD == 0, 'NACCID'}));
clinical(inds_remove, :) = [];

%% determine clinical category
%% if there is one visit with diagnosis as AD

% For example, subjects with an etiologic diagnosis of Alzheimer’s disease 
% (NACCALZD=1) can have a cognitive status of impaired-not-MCI, MCI, or dementia. 
% The only way to focus on those with AD dementia is to use both the cognitive 
% status variable (NACCUDSD=4) and the etiologic diagnosis variable (NACCALZD=1).
AD_ids = unique(clinical{clinical.NACCALZD == 1, 'NACCID'});

dementia_ids = unique(clinical{clinical.NACCUDSD > 2, 'NACCID'});
AD_34_ids = intersect(AD_ids, dementia_ids);

dementia_ids = unique(clinical{clinical.NACCUDSD == 4, 'NACCID'});
AD_4_ids = intersect(AD_ids, dementia_ids);

% nADD_ids = clinical{clinical.NACCALZD == 0, 'NACCID'};
% nADD_ids = unique(nADD_ids);
% 
% AD_ids = setdiff(AD_ids, nADD_ids);

%% if all visits have diagnosis as CN 

CN_ids = clinical{clinical.NACCUDSD == 1, 'NACCID'};
CN_ids = unique(CN_ids);

impaired_ids = unique(clinical{clinical.NACCUDSD > 1, 'NACCID'});
CN_ids = setdiff(CN_ids, impaired_ids);


%% determine Abeta positivity

% determine PET abeta positivity
pet_path = './data/UCBERKELEY_AMYLOID_MRIFREE_GAAIN_04Apr2025.csv';
pet = readtable(pet_path, VariableNamingRule='preserve');

% if there is one visit with positive abeta
amy_positive = pet{pet.AMYLOID_STATUS == 1, 'NACCID'};
amy_positive = unique(amy_positive);

% if all the visits have negative abeta
amy_negative = pet{pet.AMYLOID_STATUS == 0, 'NACCID'};
amy_negative = unique(amy_negative);
amy_negative = setdiff(amy_negative, amy_positive);


AD_abeta_missing_ids = setdiff(AD_4_ids, union(amy_positive, amy_negative));
AD_abeta_negative_ids = intersect(AD_4_ids, amy_negative);

AD_ids = intersect(AD_34_ids, amy_positive);
CN_ids = intersect(CN_ids, amy_negative);

convert_id = @(x) str2num(x(5:end));
CN_ids = cellfun(convert_id, CN_ids);
AD_ids = cellfun(convert_id, AD_ids);
AD_abeta_missing_ids = cellfun(convert_id, AD_abeta_missing_ids);
AD_abeta_negative_ids = cellfun(convert_id, AD_abeta_negative_ids);

save('population.mat', 'CN_ids', 'AD_ids', "AD_abeta_missing_ids", "AD_abeta_negative_ids");

end


function mri = load_mri()
mri_path = './data/UCDMRISBM_04Apr2025.csv';
mri = readtable(mri_path, VariableNamingRule='preserve');
mri = add_baseline_info(mri);

% remove records that are not available in the demographic table
mri(isnan(mri.years), :) = [];

mri = sortrows(mri, {'naccid','years'});
end


function mri2 = feature_selection(mri, biomarker_ver)
% combine regions
map = readtable('./FSX_biomarker_mapping.csv',VariableNamingRule='preserve');
modality = 'mri';
combined_data = combine_regions(mri, map, biomarker_ver, modality);

if ismember('HV_CTV', map.Properties.VariableNames)
    combined_data1 = combine_regions(mri, map, 'HV_CTV', modality);
    combined_data = [combined_data, combined_data1(:,{'Cortical','Hippocampus'})];
end

% add baseline info
mri1 = [mri(:,{'naccid','scandt','years'}), combined_data, mri(:,{'diagnosis_bl', ...
    'AGE_baseline','PTGENDER','APOE4','PTEDUCAT'})];

% add record age
AGE = mri1.AGE_baseline + mri1.years;
mri1 = addvars(mri1, AGE, 'Before','AGE_baseline'); 

mri2 = add_data_point_info(mri1);

end




function regress_zcore(mri, biomarker_ver)
load('population.mat');

%% concatenate the groups

mri_CN = mri(ismember(mri.RID, CN_ids), :);
mri_AD = mri(ismember(mri.RID, AD_ids), :);
mri_abeta_missing_AD = mri(ismember(mri.RID, AD_abeta_missing_ids), :);

% truncate the CN group to be age > 60 to ensures the means between the CN
% and AD groups have no significant difference
mri_CN = mri_CN(mri_CN.AGE_baseline > 60, :);


% concatenate all the subsets
data = cat(1, mri_CN, mri_AD, mri_abeta_missing_AD);
group = [zeros(size(mri_CN,1), 1); ones(size(mri_AD,1), 1); ...
    2*ones(size(mri_abeta_missing_AD,1), 1)];

data = addvars(data, group, 'before', 'AGE');


%% regress out covariates and normalize to z-score
cov_column = ["AGE","PTGENDER","APOE4","PTEDUCAT"];

[all_names, mri_bio_names, abeta_bio_names, tau_bio_names] = ...
    get_biomarker_names('FSX', biomarker_ver);
% biom_column = all_names;
mri_column = mri_bio_names(ismember(mri_bio_names, data.Properties.VariableNames));

% remove controls with missing covariates
cov_nans = any(isnan(data{:,cov_column}), 2) & data.group==0;
data(cov_nans, :) = [];

% remove data points with missing brain volumes
data(any(isnan(data{:, mri_column}), 2), :) = [];

idx_CN = data.group == 0;
idx_AD = data.group == 1;

% verify age distributions have the same mean
disp('age statistics: row 1 CN, row2 AD; col 1 mean, col 2 std');
[mean(data.AGE(idx_CN)), std(data.AGE(idx_CN)); mean(data.AGE(idx_AD)), std(data.AGE(idx_AD))]
[h,p,ci,stats] = ttest2(data.AGE(idx_CN), data.AGE(idx_AD), 'Vartype', 'unequal');
p


% there are 41 points in the case group (186 points) that have APOE4 missing. Fill 
% them with mean APOE4 in the control group such that the APOE4 slopes
% estimated from the control group will not be used in the missing points
missing_apoe_ind = data.group~=0 & isnan(data.APOE4);
data{missing_apoe_ind, 'APOE4'} = mean(data{data.group==0, 'APOE4'});

% regress out covariates
biom_add = {'Cortical','Hippocampus'};
data = regress_out_covariates(data, idx_CN, [mri_column,biom_add], ...
    [mri_column,biom_add], cov_column);

data{missing_apoe_ind, 'APOE4'} = NaN;


% normalize to zscores
data = normalize_to_zscore(data, idx_CN, mri_column, mri_column);

%% rename variables
curr_col_names = data.Properties.VariableNames;

curr_col_names{strcmp(curr_col_names, 'Hippocampus')} = 'HV';
curr_col_names{strcmp(curr_col_names, 'Cortical')} = 'CTV';

curr_col_names = strrep(curr_col_names, 'FSX-', '');

data.Properties.VariableNames = curr_col_names;

data = sortrows(data, {'RID','years'});

if  ~exist('./Result/', 'dir')
    mkdir('./Result/');
end

writetable(data, ['./Result/NACC_',biomarker_ver,'_z.csv']);

end
