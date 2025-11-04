function data_filter(data_type, biomarker_ver)

mri_fsx = 0;
mri_fsl = 0;
pet_abeta = 0;
pet_tau = 0;

%% Choose data type
if strcmp(data_type, 'mri_fsx')
    mri_fsx = 1;
end
if strcmp(data_type, 'mri_fsl')
    mri_fsl = 1;
end
if strcmp(data_type, 'pet_abeta')
    pet_abeta = 1;
end
if strcmp(data_type, 'pet_tau')
    pet_tau = 1;
end

%% Load data
if mri_fsx
    data1 = readtable('./data/UCSFFSX6_08_17_22_16Aug2023.csv',VariableNamingRule='preserve');
    data2 = readtable('./data/UCSFFSX51_11_08_19_16Aug2023.csv',VariableNamingRule='preserve');
    data3 = readtable('./data/UCSFFSX_11_02_15_24Sep2023.csv',VariableNamingRule='preserve');

%     data1.FLDSTRENG = 3 * ones(size(data1, 1), 1);
%     data1.FSVERSION = 6 * ones(size(data1, 1), 1);
% 
%     data2.FLDSTRENG = 3 * ones(size(data2, 1), 1);
%     data2.FSVERSION = 5.1 * ones(size(data2, 1), 1);
% 
%     data3.FLDSTRENG = 1.5 * ones(size(data3, 1), 1);
%     data3.FSVERSION = 4.3 * ones(size(data3, 1), 1);

    joindata = outerjoin(data1, data2,'MergeKeys',true);
    data = outerjoin(joindata, data3,'MergeKeys',true);
    map = readtable('./FSX_biomarker_mapping.csv',VariableNamingRule='preserve');
elseif mri_fsl
    data1 = readtable('./data/UCSFFSL51_03_01_22_10Oct2023.csv',VariableNamingRule='preserve');
    data2 = readtable('./data/UCSFFSL51ALL_08_01_16_10Oct2023.csv',VariableNamingRule='preserve');
    data3 = readtable('./data/UCSFFSL51Y1_08_01_16_10Oct2023.csv',VariableNamingRule='preserve');
    data4 = readtable('./data/UCSFFSL_02_01_16_10Oct2023.csv',VariableNamingRule='preserve');
    joindata1 = outerjoin(data1, data2,'MergeKeys',true);
    joindata2 = outerjoin(data3, data4,'MergeKeys',true);
    data = outerjoin(joindata1, joindata2,'MergeKeys',true);
    map = readtable('./FSL_biomarker_mapping.csv',VariableNamingRule='preserve');
elseif pet_abeta
    data = readtable('./data/UCBERKELEYAV45_8mm_02_17_23_16Aug2023.csv',VariableNamingRule='preserve');
    map = readtable('./Abeta_biomarker_mapping.csv',VariableNamingRule='preserve');
elseif pet_tau 
    data = readtable('./data/UCBERKELEYAV1451_8mm_02_17_23_16Aug2023.csv',VariableNamingRule='preserve');
    map = readtable('./TAU_biomarker_mapping.csv',VariableNamingRule='preserve');
end

%% date to years
data = add_years(data);

%% Quality control
if mri_fsx || mri_fsl
%     data = data(data.OVERALLQC == "Pass" & data.STATUS == "complete",:);
    data = data(data.OVERALLQC == "Pass",:);
    list_regions1 = ["TEMPQC", "FRONTQC", "PARQC", "INSULAQC", "OCCQC", ...
        "BGQC", "CWMQC", "VENTQC", "LHIPQC", "RHIPQC"];
    list_regions2 = ["TEMPQC", "FRONTQC", "PARQC", "INSULAQC", "OCCQC", ...
        "BGQC", "CWMQC", "VENTQC"];
    list_regions3 = ["TEMPQC", "FRONTQC", "PARQC", "INSULAQC", "OCCQC", ...
        "BGQC", "CWMQC", "VENTQC", "HIPPOQC"];
    list_regions = [list_regions1, list_regions2, list_regions3];
    list_regions = unique(list_regions);
    for region = list_regions
        if ismember(region, data.Properties.VariableNames)
            data(string(data.(region)) == "Fail", :) = [];
        end
    end

    %% Remove same data
    % Sort the table by 'VERSION' in descending order
    sortedTable = sortrows(data, 'VERSION', 'descend');

    % Find the indices of the first occurrence of each combination of 'RID' and 'VISCODE'
    [~, uniqueIndices, ~] = unique(sortedTable(:, {'RID', 'years'}));

    % Extract the rows with the latest 'VERSION' for each combination of 'RID' and 'VISCODE'
    data = sortedTable(uniqueIndices, :);
end

%% Biomarkers Mapping

if mri_fsx || mri_fsl
    modality = 'mri';
elseif pet_abeta || pet_tau
    modality = 'pet';
end
combined_data = combine_regions(data, map, biomarker_ver, modality);

if ismember('HV_CTV', map.Properties.VariableNames)
    combined_data1 = combine_regions(data, map, 'HV_CTV', modality);
    combined_data = [combined_data, combined_data1(:,{'Cortical','Hippocampus'})];
end


% normalize
% if pet_abeta
%     l = size(combined_data,2);
%     combined_data{:,3:l} = combined_data{:,3:l} ./ data.WHOLECEREBELLUM_SUVR;
% end

myTable = [num2cell(data.RID), num2cell(data.years), combined_data];
% if ismember('IMAGEUID', data.Properties.VariableNames)
%     myTable = [myTable,num2cell(data.IMAGEUID)];
%     myTable = renamevars(myTable, myTable.Properties.VariableNames(end), 'IMAGEUID');
% end
myTable = renamevars(myTable, myTable.Properties.VariableNames(1:2), {'RID','years'});

% FLDSTRENG = data.FLDSTRENG;
% FSVERSION = data.FSVERSION;
% myTable = addvars(myTable, FLDSTRENG, FSVERSION);

% myTable(isnan(years), :) = [];
% myTable = unique(myTable);


%% Output
if mri_fsx
    writetable(myTable, sprintf('./Result/UCSFFSX_%s_all.csv',biomarker_ver));
elseif mri_fsl
    writetable(myTable, sprintf('./Result/UCSFFSL_%s_all.csv',biomarker_ver));
elseif pet_abeta
    writetable(myTable, sprintf('./Result/UCBERKELEYAV45_%s_all.csv',biomarker_ver));
elseif pet_tau
    writetable(myTable, sprintf('./Result/UCBERKELEYAV1451_%s_all.csv',biomarker_ver));
end

end

