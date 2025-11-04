function copy_files_to_destination_NACC(biomarker_vers)
output_dir = '../../input/NACC/';
if  ~exist(output_dir, 'dir')
    mkdir(output_dir);
end


for j = 1:length(biomarker_vers)
    filename = ['NACC_',biomarker_vers{j},'_z.csv'];
    source_file = ['./Result/', filename];
    des_file = [output_dir, filename];
    [status, msg] = copyfile(source_file, des_file);
end


clinical_file = 'nacc_clinical.csv';
copyfile(['./data/',clinical_file], [output_dir,clinical_file]);

if exist([output_dir,'NACC_time2event.csv'], 'file')
    delete([output_dir,'NACC_time2event.csv']);
end

end

