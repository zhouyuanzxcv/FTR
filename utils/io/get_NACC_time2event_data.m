function data = get_NACC_time2event_data()
%GET_NACC_TIME2EVENT_DATA Summary of this function goes here
%   Detailed explanation goes here
save_file = 'input/NACC/NACC_time2event.csv';
if ~exist(save_file, 'file')
    clinical_path = ['input/NACC/nacc_clinical.csv'];
    mmse = readtable(clinical_path,'VariableNamingRule','preserve');

    mmse = mmse(:,["NACCID","NACCFDYS","NACCUDSD"]);
    rid = cellfun(@(x) str2num(x(5:end)), mmse.NACCID, 'UniformOutput', 1);
    mmse.NACCID = rid;
    mmse.NACCFDYS = mmse.NACCFDYS / 365;
    mmse = renamevars(mmse, 'NACCID', 'RID');
    mmse = renamevars(mmse, 'NACCFDYS', 'years');

    CN_inds = mmse.NACCUDSD == 1 | mmse.NACCUDSD == 2;
    MCI_inds = mmse.NACCUDSD == 3;
    AD_inds = mmse.NACCUDSD == 4;
    mmse{CN_inds, 'NACCUDSD'} = 0;
    mmse{MCI_inds, 'NACCUDSD'} = 0.5;
    mmse{AD_inds, 'NACCUDSD'} = 1;

    mmse = renamevars(mmse, 'NACCUDSD', 'diagnosis');
    
    data = convert_timeseries_to_timetoevent(mmse);
    
%     assert(all(data_bl.years == 0));
    data = renamevars(data, 'diagnosis', 'dx_time2event');
    data = renamevars(data, 'years', 'years_time2event');
        
    writetable(data, save_file);
else
    data = readtable(save_file,'VariableNamingRule','preserve');
end

end

