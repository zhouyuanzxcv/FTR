function copy_files_to_destination(freesurfer_vers, biomarker_vers)
output_dir = '../../input/ADNI/';
if  ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

for i = 1:length(freesurfer_vers)
    for j = 1:length(biomarker_vers)   
        filename = ['ADNI_demo_MRI_',freesurfer_vers{i},'_PET_',biomarker_vers{j},'_zscore.csv'];
        source_file = ['./Result/', filename];
        des_file = [output_dir, filename];
        [status, msg] = copyfile(source_file, des_file);
    end
end

adnimerge_file = 'ADNIMERGE_16Aug2023.csv';
copyfile(['./data/',adnimerge_file], [output_dir,adnimerge_file]);

if exist([output_dir,'ADNI_time2event.csv'], 'file')
    delete([output_dir,'ADNI_time2event.csv']);
end

end

