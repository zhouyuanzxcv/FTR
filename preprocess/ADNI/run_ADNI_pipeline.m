function [outputArg1,outputArg2] = run_ADNI_pipeline(inputArg1,inputArg2)
%ADNI_PIPELINE Summary of this function goes here
%   Detailed explanation goes here

freesurfer_vers = {'FSX'};
biomarker_vers = {'HM','HS','LM','LS'};

if 1
for j = 1:length(biomarker_vers)
    biomarker_ver = biomarker_vers{j};
    data_filter('mri_fsx', biomarker_ver);
    %data_filter('mri_fsl', biomarker_ver);
    data_filter('pet_abeta', biomarker_ver);
    %data_filter('pet_tau', biomarker_ver);
end
end


if 1
for i = 1:length(freesurfer_vers)
    for j = 1:length(biomarker_vers)       
        combine_table(freesurfer_vers{i}, biomarker_vers{j});
    end
end
end

for i = 1:length(freesurfer_vers)
    for j = 1:length(biomarker_vers)        
        regress_normalize(freesurfer_vers{i}, biomarker_vers{j});
    end
end

% copy the zscore files to input
copy_files_to_destination(freesurfer_vers, biomarker_vers);

end


