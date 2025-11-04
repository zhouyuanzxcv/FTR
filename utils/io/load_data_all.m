function [data,mmse,time2event,dx_sel] = load_data_all(dataset_name)
%LOAD_DATA_PATH_MMSE_TIME2EVENT Summary of this function goes here
%   Detailed explanation goes here
dx_sel = false;
switch dataset_name
    case {'ADNI_FSX_LM','ADNI_FSX_LS','ADNI_FSX_HM','ADNI_FSX_HS',...
            'ADNI_FSX_LM_DS','ADNI_FSX_LS_DS','ADNI_FSX_HM_DS','ADNI_FSX_HS_DS',...
            'ADNI_FSX_1.5_HM', 'ADNI_FSX_1.5_HS', 'ADNI_FSX_1.5_LM', 'ADNI_FSX_1.5_LS', ...
            'ADNI_FSX_3_HM', 'ADNI_FSX_3_HS', 'ADNI_FSX_3_LM', 'ADNI_FSX_3_LS'}

        fs_ver = dataset_name(6:8);
        dx_sel = false;
        fldstr = '';

        if length(dataset_name) == 11 % 'ADNI_FSX_HM'
            biom_ver = dataset_name(10:11);
        elseif length(dataset_name) == 15 % 'ADNI_FSX_1.5_HM'
            biom_ver = dataset_name(end-1:end);
            fldstr = '1.5_';
        elseif length(dataset_name) == 13 % 'ADNI_FSX_3_HM'
            biom_ver = dataset_name(end-1:end);
            fldstr = '3_';
        elseif length(dataset_name) == 14 % 'ADNI_FSX_HM_DS'
            biom_ver = dataset_name(end-4:end-3);
            dx_sel = true;
        end

        data_path = ['input/ADNI/ADNI_demo_MRI_',fldstr,fs_ver, '_PET_', ...
            biom_ver,'_zscore.csv'];


        % load MMSE for ADNI
        clinical_path = ['input/ADNI/ADNIMERGE_16Aug2023.csv'];
        clinical_table = readtable(clinical_path,'VariableNamingRule','preserve');
        use_metric = 'MMSE';
        if strcmp(use_metric, 'MMSE') % use MMSE
            mmse = clinical_table(:,["RID","Years_bl","MMSE","MMSE_bl", ...
                "MOCA","MOCA_bl","CDRSB","CDRSB_bl"]);
            mmse = renamevars(mmse, 'Years_bl', 'years');
        elseif strcmp(use_metric, 'MOCA') % check MOCA
            mmse = clinical_table(:,["RID","Years_bl","MOCA","MOCA_bl"]);
            mmse = renamevars(mmse, 'Years_bl', 'years');
            mmse = renamevars(mmse, 'MOCA', "MMSE");
            mmse = renamevars(mmse, 'MOCA_bl', 'MMSE_bl');
        elseif strcmp(use_metric, 'CDRSUM')
            mmse = clinical_table(:,["RID","Years_bl","CDRSB","CDRSB_bl"]);
            mmse = renamevars(mmse, 'Years_bl', 'years');
            mmse = renamevars(mmse, 'CDRSB', "MMSE");
            mmse = renamevars(mmse, 'CDRSB_bl', 'MMSE_bl');
            mmse.MMSE = -mmse.MMSE;
            mmse.MMSE_bl = -mmse.MMSE_bl;
        end

        time2event = get_ADNI_time2event_data();

    case {'OASIS3_ROD0_HM','OASIS3_ROD0_HS', ...
            'OASIS3_ROD1_HM','OASIS3_ROD1_HS','OASIS3_ROD1_LM','OASIS3_ROD1_LS',...
            'OASIS3_ROD2_HM','OASIS3_ROD2_HS','OASIS3_ROD2_LM','OASIS3_ROD2_LS',...
            'OASIS3_ROD3_HM','OASIS3_ROD3_HS','OASIS3_ROD3_LM','OASIS3_ROD3_LS'}
        data_path = ['input/OASIS3/OASIS3_input_',dataset_name(end-6:end),'_z.csv'];

        % load MMSE for OASIS (it has been ordered by id and years)
        mmse_path = ['input/OASIS3/OASIS3_UDSb4_cdr.csv'];
        mmse = readtable(mmse_path,'VariableNamingRule','preserve');
        mmse = mmse(:,["OASISID","days_to_visit","MMSE","CDRSUM"]);
        rid = cellfun(@(x) str2num(x(5:end)), mmse.OASISID, 'UniformOutput', 1);
        mmse.OASISID = rid;
        mmse.days_to_visit = mmse.days_to_visit / 365;
        mmse = renamevars(mmse, 'OASISID', 'RID');
        mmse = renamevars(mmse, 'days_to_visit', 'years');
        mmse = renamevars(mmse, 'CDRSUM', 'CDRSB');
        % add MMSE bl
        [~,ia,ic] = unique(mmse.RID);
        mmse_bl = mmse(ia, {'RID', 'MMSE', 'CDRSB', 'years'});
        ind_bl_missing = mmse_bl.years >= 0.5;
        mmse_bl(ind_bl_missing, :) = [];
        mmse_bl = renamevars(mmse_bl, 'MMSE', 'MMSE_bl');
        mmse_bl = renamevars(mmse_bl, 'CDRSB', 'CDRSB_bl');
        mmse_bl = removevars(mmse_bl, 'years');
        mmse = outerjoin(mmse, mmse_bl,'Type','Left','Keys','RID','MergeKeys',true);

        time2event = get_OASIS_time2event_data();
    case {'NACC_HM','NACC_HS','NACC_LM','NACC_LS'}
        data_path = ['input/NACC/',dataset_name,'_z.csv'];

        % load MMSE
        mmse_path = ['input/NACC/nacc_clinical.csv'];
        mmse = readtable(mmse_path,'VariableNamingRule','preserve');
        mmse = sortrows(mmse, {'NACCID','NACCFDYS'});

        % for NACC, use MOCA, MMSE, or CDRSUM
        
        mmse = mmse(:,{'NACCID','NACCFDYS','CDRSUM','NACCMMSE','NACCMOCA'});
%         if strcmp(use_metric, 'NACCMMSE')
        mmse_nan_inds = mmse.NACCMMSE == -4 | mmse.NACCMMSE == 88 | ...
            (mmse.NACCMMSE >= 95 & mmse.NACCMMSE <= 98);
        mmse.NACCMMSE(mmse_nan_inds) = NaN;
%         elseif strcmp(use_metric, 'NACCMOCA')
        mmse_nan_inds = mmse.NACCMOCA == 88 | mmse.NACCMOCA == 99 | ...
            mmse.NACCMOCA == -4;
        mmse.NACCMOCA(mmse_nan_inds) = NaN;
%         elseif strcmp(use_metric, 'CDRSUM')
        mmse.CDRSUM = -mmse.CDRSUM;
%         end

        rid = cellfun(@(x) str2num(x(5:end)), mmse.NACCID, 'UniformOutput', 1);
        mmse.NACCID = rid;
        mmse.NACCFDYS = mmse.NACCFDYS / 365;
        mmse = renamevars(mmse, 'NACCID', 'RID');
        mmse = renamevars(mmse, 'NACCFDYS', 'years');

        use_metric = 'CDRSUM';
        CDRSB = -mmse.(use_metric);

        mmse = renamevars(mmse, use_metric, 'MMSE');  
        
        mmse = addvars(mmse, CDRSB, 'after', 'MMSE');
        mmse = renamevars(mmse, 'NACCMOCA', 'MOCA');

        % add MMSE bl
        [~,ia,ic] = unique(mmse.RID);
        mmse_bl = mmse(ia, {'RID', 'MMSE', 'CDRSB', 'MOCA', 'NACCMMSE', 'years'});
        ind_bl_missing = mmse_bl.years >= 0.5;
        mmse_bl(ind_bl_missing, :) = [];
        mmse_bl = renamevars(mmse_bl, {'MMSE','CDRSB','MOCA','NACCMMSE'}, ...
            {'MMSE_bl','CDRSB_bl','MOCA_bl','NACCMMSE_bl'});
        mmse_bl = removevars(mmse_bl, 'years');
        mmse = outerjoin(mmse, mmse_bl,'Type','Left','Keys','RID','MergeKeys',true);

        time2event = get_NACC_time2event_data();
    otherwise
end


data = readtable(data_path,'VariableNamingRule','preserve');
if strcmp(dataset_name(1:4), 'ADNI')
    data.Properties.VariableNames = strrep(data.Properties.VariableNames, 'TAU', 'Tau');
end
if strcmp(dataset_name(1:5), 'OASIS')
    data = renamevars(data, 'APOE', 'APOE4');
end
if strcmp(dataset_name(1:4), 'NACC')
    % For NACC, change CDRSUM to MMSE and make it negative
    data = removevars(data, 'MMSE');
    data = renamevars(data, 'CDRSUM', 'MMSE');
    data.MMSE = -data.MMSE;
end

end

