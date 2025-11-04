function data = get_ADNI_time2event_data()
%GET_ADNI_TIME2EVENT_DATA Summary of this function goes here
%   Detailed explanation goes here
save_file = 'input/ADNI/ADNI_time2event.csv';
if ~exist(save_file, 'file')
    clinical_path = ['input/ADNI/ADNIMERGE_16Aug2023.csv'];
    clinical_table = readtable(clinical_path,'VariableNamingRule','preserve');
    
    data = clinical_table(:,["RID","Years_bl","DX"]);
    
    data = renamevars(data, 'Years_bl', 'years');
    
    % to run conversion, DX should be renamed to diagnosis
    data.DX = convert_diagnosis_to_number(data.DX);
    data = renamevars(data, 'DX', 'diagnosis');
    data = convert_timeseries_to_timetoevent(data);
    
%     assert(all(data_bl.years == 0));
    data = renamevars(data, 'diagnosis', 'dx_time2event');
    data = renamevars(data, 'years', 'years_time2event');
        
    writetable(data, save_file);
else
    data = readtable(save_file,'VariableNamingRule','preserve');
end

end

function DX_new = convert_diagnosis_to_number(DX)
DX_new = NaN(size(DX,1), 1);
DX_new(strcmp(DX, 'Dementia')) = 1;
DX_new(strcmp(DX, 'MCI')) = 0.5;
DX_new(strcmp(DX, 'CN')) = 0;
end

