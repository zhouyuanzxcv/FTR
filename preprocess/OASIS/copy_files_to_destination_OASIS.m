function copy_files_to_destination_OASIS()
biomarker_vers = {'HM','HS','LM','LS'};

output_dir = '../../input/OASIS3/';
if  ~exist(output_dir, 'dir')
    mkdir(output_dir);
end


for j = 1:length(biomarker_vers)   
    filename = ['OASIS3_input_ROD1_',biomarker_vers{j},'_z.csv'];
    source_file = ['./Result/', filename];
    des_file = [output_dir, filename];
    [status, msg] = copyfile(source_file, des_file);
end


cdr_file = 'OASIS3_UDSb4_cdr.csv';
copyfile(['./data/',cdr_file], [output_dir,cdr_file]);

if exist([output_dir,'OASIS3_time2event.csv'], 'file')
    delete([output_dir,'OASIS3_time2event.csv']);
end

end

