function colors = get_method_color(method_names)
%GET_METHOD_COLOR Summary of this function goes here
%   Detailed explanation goes here
myMap = containers.Map;

hcbl = sprintf('Hierarchical\nClustering\n');
method_labels = {'FTR','FTR (MCEM)','SuStaIn','CTV-HV-Ratio',hcbl};
% method_colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330]};
method_colors = {[0 0.4470 0.7410], [0 0.4470 0.7410], [0.8500 0.3250 0.0980], ...
    [0.9290 0.6940 0.1250], [0.4660 0.6740 0.1880]};

method_labels{end+1} = 'Hippocampus volume';
method_colors{end+1} = [67,217,22]/255;

method_labels{end+1} = 'CSF Abeta';
method_colors{end+1} = [217,174,22]/255;

method_labels{end+1} = 'CSF pTau';
method_colors{end+1} = [200,217,22]/255;

method_labels{end+1} = 'MMSE';
method_colors{end+1} = [47,217,22]/255;

for i = 1:length(method_labels)
    myMap(method_labels{i}) = method_colors{i};
end

colors = [];
for i = 1:length(method_names)
    color = myMap(method_names{i});
    colors = [colors; color];
end

end

