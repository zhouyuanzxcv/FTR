function [outputArg1,outputArg2] = analyze_reproducibility(inputArg1,inputArg2)
%ANALYZE_REPRODUCIBILITY Summary of this function goes here
%   Detailed explanation goes here
close all

atlases = {{'HS','HM'},{'HS','LS'},{'HS','LM'},{'HM','LS'},{'HM','LM'},{'LM','LS'}};

datasets = {'ADNI_FSX', 'OASIS3_ROD1', 'NACC'};

f = figure;
f.Position = get_figure_position();
xlabel('Comparison between atlases')
ylabel('Consistency (%)')
hold on;

xlim([0.5 length(atlases)+0.5])
ylim([40 100])

xticks(1:length(atlases))
xticklabels({'82 vs. 41', '82 vs. 26', '82 vs. 13', '41 vs. 26', '41 vs. 13', '26 vs. 13'})


colors = {[18,53,85] / 255, [230,155,3] / 255, [25,227,66] / 255};
marker_shape = {'o', 'diamond', 's'};

marker_size = 80;

ps = [];

for i = 1:length(atlases)
    at = atlases{i};
    for j = 1:length(datasets)

        data1 = [datasets{j}, '_', at{1}];
        data2 = [datasets{j}, '_', at{2}];

        data = {data1, data2};
        method = {'FTR_MCEM'};
        nsubtypes = [3;3];
        results = load_data_result(data, method, nsubtypes, 0, 1);
        joindata1 = results{1,1}.joindata;
        joindata2 = results{2,1}.joindata;
        subtype1 = dummyvar(joindata1.subtype);
        subtype2 = dummyvar(joindata2.subtype);
    
        [~, subtype2] = permute_endmembers(subtype1', subtype2');
        
        [~, subtype2] = max(subtype2', [], 2);
    
        joindata2.subtype = subtype2;
    
        conf_matrix = get_confusion_matrix(joindata1, joindata2, 3);
        consistency = sum(diag(conf_matrix)) / sum(sum(conf_matrix)) * 100;
        
        p = scatter(i, consistency, marker_size, colors{j}, "filled", marker_shape{j});
%         p = scatter(i, consistency, marker_size, "filled", marker_shape{j});
        p.MarkerFaceAlpha = 0.8;
        if i == 1
            ps = [ps, p];
        end

    end
end


legend(ps, {'ADNI', 'OASIS', 'NACC'}, 'Location', 'best')
export_fig './figures/all/4G.jpg' -r500 -transparent
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
