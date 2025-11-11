function [outputArg1,outputArg2] = analyze_reproducibility(inputArg1,inputArg2)
%ANALYZE_REPRODUCIBILITY Summary of this function goes here
%   Detailed explanation goes here
close all

data_names = {{'ADNI_FSX_HS','ADNI_FSX_HM'},{'ADNI_FSX_HS','ADNI_FSX_LS'},{'ADNI_FSX_HS','ADNI_FSX_LM'}};

method = {'FTR_MCEM'};

nsubtypes = [3;3];


f = figure;
f.Position = get_figure_position(1);

tiledlayout(2,2, "TileSpacing", "compact")

for plotIdx = 1:length(data_names)
    
    data = data_names{plotIdx};
    
    results = load_data_result(data, method, nsubtypes, 0, 1);
    
    biomarker_names = results{1,1}.biomarker_names;
    
    [nbiom,num_int,~] = size(results{1,1}.mdl.re_traj);
    stage = (0:num_int-1)/(num_int-1);
    
    num_datasets = length(data);
    
    joindata1 = results{1,1}.joindata;
    joindata2 = results{2,1}.joindata;


    conf_matrix = get_confusion_matrix(joindata1, joindata2, 3);
    

    if plotIdx == 1
        nexttile([1,2])
    else
        nexttile([1,1])
    end

    show_confusion_matrix(conf_matrix, plotIdx)
    

end

export_fig './figures/all/4F.jpg' -r500 -transparent
end

function conf_matrix = get_confusion_matrix(train_subtype_stage, test_subtype_stage, nsubtype)

conf_matrix = zeros(nsubtype);
for i = 1:nsubtype
    for j = 1:nsubtype
        % Find PTIDs of train_data with subtype i and test_data with subtype j
        train_PTID_subtype_i = unique(train_subtype_stage.RID(train_subtype_stage.subtype ==  i));
        test_PTID_subtype_j = unique(test_subtype_stage.RID(test_subtype_stage.subtype ==  j));

        % Count the number of common PTIDs between train and test subtypes
        common_PTIDs = intersect(train_PTID_subtype_i, test_PTID_subtype_j);
        conf_matrix(i, j) = numel(common_PTIDs);
    end
end

end

function [outputArg1,outputArg2] = show_confusion_matrix(conf_matrix, plotIdx)
%SHOW_CONFUSION_MAT Summary of this function goes here
%   Detailed explanation goes here
 % Display confusion matrix

perc_consistency = sum(diag(conf_matrix)) / sum(sum(conf_matrix));

% Example subtypes 
num_subtypes = size(conf_matrix, 1);
train_subtypes = {'S1', 'S2', 'S3'};
test_subtypes = {'S1', 'S2', 'S3'};
labels = {'41', '26', '13'};

% Plotting the confusion matrix
% colormap sky
colormap(get_colormap_sky());
imagesc(conf_matrix);
axis equal;
xlim([0.5 3.5])
ylim([0.5 3.5])

% Adjusting axis labels and ticks
set(gca, 'XTick', 1:length(train_subtypes), 'XTickLabel', train_subtypes, 'YTick', 1:length(test_subtypes), 'YTickLabel', test_subtypes);
xlabel([labels{plotIdx},' regions']);
ylabel('82 regions');

% title({['ADNI 82 vs. ', labels{plotIdx}], ['Consistency: ',num2str(perc_consistency*100,3),'%'], ''});
title({['Consistency: ',num2str(perc_consistency*100,3),'%']});

% Displaying values in the cells (optional)
textStrings = num2str(conf_matrix(:), '%d');
textStrings = strtrim(cellstr(textStrings));
[x, y] = meshgrid(1:length(train_subtypes), 1:length(test_subtypes));
hStrings = text(x(:), y(:), textStrings(:), 'HorizontalAlignment', 'center');
midValue = mean(get(gca, 'CLim'));

end