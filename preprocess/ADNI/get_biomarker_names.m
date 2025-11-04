function [all_names, mri_bio_names, abeta_bio_names, tau_bio_names] = ...
    get_biomarker_names(freesurfer_ver, biomarker_ver)
%GET_BIOMARKER_NAMES Summary of this function goes here
%   Detailed explanation goes here
if strcmp(freesurfer_ver, 'FSX')
    map_mri = readtable('./FSX_biomarker_mapping.csv',VariableNamingRule='preserve');
elseif strcmp(freesurfer_ver, 'FSL')
    map_mri = readtable('./FSL_biomarker_mapping.csv',VariableNamingRule='preserve');
end

map_abeta = readtable('./Abeta_biomarker_mapping.csv',VariableNamingRule='preserve');
map_tau = readtable('./TAU_biomarker_mapping.csv',VariableNamingRule='preserve');

mri_bio_names = map_mri.(biomarker_ver)(logical(map_mri.('REQUIRE_NORMALIZE')));
mri_bio_names = unique(mri_bio_names);

abeta_bio_names = map_abeta.(biomarker_ver)(logical(map_abeta.('REQUIRE_NORMALIZE')));
abeta_bio_names = unique(abeta_bio_names);

tau_bio_names = map_tau.(biomarker_ver)(logical(map_tau.('REQUIRE_NORMALIZE')));
tau_bio_names = unique(tau_bio_names);

mri_bio_names = mri_bio_names';
abeta_bio_names = abeta_bio_names';
tau_bio_names = tau_bio_names';

all_names = [mri_bio_names, abeta_bio_names, tau_bio_names];

end

