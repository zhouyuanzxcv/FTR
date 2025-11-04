function [outputArg1,outputArg2] = change_excel_contents(inputArg1,inputArg2)
%CHANGE_EXCEL_CONTENTS Summary of this function goes here
%   Detailed explanation goes here
map = readtable('./TAU/TAU_biomarker_mapping.csv',VariableNamingRule='preserve');
biom = map.LS;
biom1 = biom; 
for i = 1:length(biom)
    if i >= 1 && i <= 41
        biom1{i} = strrep(biom{i}, 'TAU-', 'TAU-Left-'); 
    elseif i >= 42 && i <= 82
        biom1{i} = strrep(biom{i}, 'TAU-', 'TAU-Right-'); 
    end
end
LS = biom1; 
map.LS = LS;
writetable(map, './TAU/TAU_biomarker_mapping1.csv');
end

