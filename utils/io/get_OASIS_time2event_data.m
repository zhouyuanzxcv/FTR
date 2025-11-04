function data = get_OASIS_time2event_data()
%GET_ADNI_TIME2EVENT_DATA Summary of this function goes here
%   Detailed explanation goes here
save_file = 'input/OASIS3/OASIS3_time2event.csv';
if ~exist(save_file, 'file')
    clinical_path = ['input/OASIS3/OASIS3_UDSb4_cdr.csv'];
    mmse = readtable(clinical_path,'VariableNamingRule','preserve');

    mmse = mmse(:,["OASISID","days_to_visit","CDRTOT"]);
    rid = cellfun(@(x) str2num(x(5:end)), mmse.OASISID, 'UniformOutput', 1);
    mmse.OASISID = rid;
    mmse.days_to_visit = mmse.days_to_visit / 365;
    mmse = renamevars(mmse, 'OASISID', 'RID');
    mmse = renamevars(mmse, 'days_to_visit', 'years');

    mmse{mmse.CDRTOT>1, 'CDRTOT'} = 1;

    mmse = renamevars(mmse, 'CDRTOT', 'diagnosis');
    
    data = convert_timeseries_to_timetoevent(mmse);
    
%     assert(all(data_bl.years == 0));
    data = renamevars(data, 'diagnosis', 'dx_time2event');
    data = renamevars(data, 'years', 'years_time2event');
        
    writetable(data, save_file);
else
    data = readtable(save_file,'VariableNamingRule','preserve');
end

end

