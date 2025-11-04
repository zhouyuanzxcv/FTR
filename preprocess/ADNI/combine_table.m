function combine_table(freesurfer_ver, biomarker_ver)
%% Input
mri_fsx = readtable(sprintf('./Result/UCSFFSX_%s_all.csv',biomarker_ver),VariableNamingRule='preserve');
% mri_fsx.ICV = [];
% mri_fsl = readtable(sprintf('./Result/UCSFFSL_%s_all.csv',biomarker_ver),VariableNamingRule='preserve');
% mri_fsl.ICV = [];
pet_abeta = readtable(sprintf('./Result/UCBERKELEYAV45_%s_all.csv',biomarker_ver),VariableNamingRule='preserve');
% pet_tau = readtable(sprintf('./Result/UCBERKELEYAV1451_%s_all.csv',biomarker_ver),VariableNamingRule='preserve');
data = readtable('./data/ADNIMERGE_16Aug2023.csv',VariableNamingRule='preserve');
old_csf =  readtable('./data/UPENNBIOMK_MASTER_FINAL_04Oct2023.csv',VariableNamingRule='preserve');

demo_sel = data(:,["RID","VISCODE","EXAMDATE","DX","AGE","PTGENDER","APOE4",...
    "PTEDUCAT","ADAS11","ADAS13","MMSE","FLDSTRENG","FSVERSION"]);
old_csf_sel = old_csf(:,["RID","VISCODE","VISCODE2","EXAMDATE","ABETA42","TAU","PTAU","RUNDATE"]);

old_csf_sel = add_years(old_csf_sel);
demo_sel = add_years(demo_sel);
demo_sel = add_diagnosis_bl(demo_sel);

demo_sel = removevars(demo_sel, {'VISCODE', 'EXAMDATE'});
old_csf_sel = removevars(old_csf_sel, ["VISCODE", "VISCODE2", "EXAMDATE"]);

AGE_baseline = demo_sel.AGE;
demo_sel = addvars(demo_sel, AGE_baseline,'After','AGE');
demo_sel.AGE = demo_sel.AGE + demo_sel.years;

old_csf_sel = renamevars(old_csf_sel, 'ABETA42', 'ABETA');

%% remove same RID & VISCODE in csf file
% Sort the table by 'VERSION' in descending order
sortedTable = sortrows(old_csf_sel, 'RUNDATE', 'descend');

% Find the indices of the first occurrence of each combination of 'RID' and 'VISCODE'
[~, uniqueIndices, ~] = unique(sortedTable(:, {'RID', 'years'}));

% Extract the rows with the latest 'VERSION' for each combination of 'RID' and 'VISCODE'
old_csf_sel = sortedTable(uniqueIndices, :);
old_csf_sel = removevars(old_csf_sel, 'RUNDATE');

%% Combine data
if strcmp(freesurfer_ver, 'FSX')
    joindata1 = mri_fsx;
elseif strcmp(freesurfer_ver, 'FSL')
    joindata1 = mri_fsl;    
end

nanRows = any(isnan(table2array(joindata1)), 2);
joindata1(nanRows, :) = [];

% demo_sel = renamevars(demo_sel, 'FLDSTRENG', 'ADNIMERGE_FLDSTR');
% demo_sel = renamevars(demo_sel, 'FSVERSION', 'ADNIMERGE_FSVERSION');
joindata1 = join_by_RID_years(joindata1, pet_abeta, 'insert_right_into_left');
% joindata1 = join_by_RID_years(joindata1, pet_tau, 'insert_right_into_left');
joindata1 = join_by_RID_years(joindata1, demo_sel, 'find_right_for_each_left');
joindata1 = join_by_RID_years(joindata1, old_csf_sel, 'insert_right_into_left');

data = joindata1;

% Replace "Male" with 1 and "Female" with 0 in the 'PTGENDER' column
% data.PTGENDER = cellfun(@(x) strcmp(x, 'Male'), data.PTGENDER);

new_genders = NaN(size(data,1), 1);
for i = 1:size(data, 1)
    gender_i = data{i, "PTGENDER"}{:};
    if ~isempty(gender_i)
        new_genders(i) = strcmp(gender_i, 'Male');
    end
end
data.PTGENDER = new_genders;

% Change DX to diagnosis
% data.DX = replace(data.DX, {'Dementia', 'MCI', 'CN'}, {'1', '0.5', '0'});
% data.DX = str2double(data.DX);

DX_new = convert_diagnosis_to_number(data.DX);
data.DX = DX_new;
data = renamevars(data, 'DX', 'diagnosis');

data.diagnosis_bl = convert_diagnosis_to_number(data.diagnosis_bl);

% replace scanner strength
data.FLDSTRENG(cellfun(@isempty, data.FLDSTRENG)) = {''};
replaceMap = containers.Map({'1.5 Tesla MRI', '3 Tesla MRI', ''}, [1.5, 3, NaN]);
data.FLDSTRENG = cellfun(@(x) replaceMap(x), data.FLDSTRENG);

% replace freesurfer version
data.FSVERSION(cellfun(@isempty, data.FSVERSION)) = {''};
replaceMap = containers.Map({'Cross-Sectional FreeSurfer (FreeSurfer Version 4.3)',...
    'Cross-Sectional FreeSurfer (5.1)', 'Cross-Sectional FreeSurfer (6.0)', ''}, ...
    [4.3, 5.1, 6.0, NaN]);
data.FSVERSION = cellfun(@(x) replaceMap(x), data.FSVERSION);

%% remove nan
% rowsToDelete1 = all(isnan(data{:, 16:181}), 2);
% data(rowsToDelete1, :) = [];
% rowsToDelete2 = any(isnan(data{:, 3:8}), 2);
% data(rowsToDelete2, :) = [];
% data = unique(data);

writetable(data,['./Result/ADNI_demo_MRI_', freesurfer_ver, '_PET_',biomarker_ver,'_all.csv']);
end

function DX_new = convert_diagnosis_to_number(DX)
DX_new = NaN(size(DX,1), 1);
DX_new(strcmp(DX, 'Dementia')) = 1;
DX_new(strcmp(DX, 'MCI')) = 0.5;
DX_new(strcmp(DX, 'CN')) = 0;
end