function [outputArg1,outputArg2] = change_excel_contents(inputArg1,inputArg2)
%CHANGE_EXCEL_CONTENTS Summary of this function goes here
%   Detailed explanation goes here
map = readtable('./FSX_biomarker_mapping.csv',VariableNamingRule='preserve');
biom = map.FLDNAME;
biom1 = biom; 
for i = 1:length(biom)
    if contains(biom{i}, 'FSX-Left-')
        biom1{i} = strrep(biom{i}, 'FSX-Left-', 'lh-'); 
    elseif contains(biom{i}, 'FSX-Right-')
        biom1{i} = strrep(biom{i}, 'FSX-Right-', 'rh-'); 
    end
end

for i = 1:length(biom1)
    if ~contains(biom1{i}, '-gvol')
        biom1{i} = [biom1{i},'-gvol'];
    end
end

map.FLDNAME = biom1;
writetable(map, './FSX_biomarker_mapping.csv');
end

