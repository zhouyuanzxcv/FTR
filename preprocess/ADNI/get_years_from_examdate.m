function [years, diagnosis_bl] = get_years_from_examdate(data)
%GET_YEARS_FROM_EXAMDATE Summary of this function goes here
%   Detailed explanation goes here

data_adni_merge = readtable('./data/ADNIMERGE_16Aug2023.csv',VariableNamingRule='preserve');

bl_dates = data_adni_merge(data_adni_merge.VISCODE == "bl", ["RID","EXAMDATE","DX"]);

converted_dates = zeros(size(data,1), 1);
diagnosis_bl = cell(size(data,1), 1);

for i = 1:size(data, 1)
    ind = data{i, "RID"} == bl_dates{:, "RID"};
    bl_date_i = bl_dates{ind, "EXAMDATE"};
    bl_dx = bl_dates{ind, "DX"};
    
    if isempty(bl_dx)
        diagnosis_bl(i) = {''};
    else
        diagnosis_bl(i) = bl_dx;
    end
    
    if isempty(bl_date_i)
        converted_dates(i) = nan;
    else
        % if a subject has 2 screening dates, choose the earlies one
        bl_date_i = bl_date_i(1);
    
        curr_date_i = data{i, "EXAMDATE"};
        years_diff = yearfrac(bl_date_i, curr_date_i, 0);
        
        % round it to half a year
%         converted_dates(i) = round(years_diff * 2) / 2;
        converted_dates(i) = years_diff;
    end
end

% there are 53 subjects whose screening date is more than 3 months before
% baseline. Their years become -0.5.
years = converted_dates;
end

