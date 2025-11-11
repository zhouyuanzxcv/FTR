function [rs, ris] = analyze_cv()
%ANALYZE_CV
close all

output_file_name = 'ADNI_FSX_HS_FTR_MCEM';


cv_loglikelihood = readmatrix(['output/', output_file_name, '/loglikelihood_cross_validation.csv']);
% S11A
f = figure();
f.Position = get_figure_position();
plot(sum(cv_loglikelihood,1));
xlabel('Number of subtypes');
ylabel('Loglikelihood from cross validation');
xlim([1 6])

export_fig './figures/all/S11A.jpg' -r500 -transparent
% close all

max_nsubtype = size(cv_loglikelihood, 2);

num_fold = size(cv_loglikelihood, 1);

check_inds = 2:max_nsubtype;


max_len_cef = 0;
max_len_rand = 0;
for j = check_inds
    [cef_mat, rand_mat] = cal_stability(output_file_name, j, 0, num_fold);
    max_len_cef = max(max_len_cef, numel(cef_mat));  % 记录cef_mat的最大长度
    max_len_rand = max(max_len_rand, numel(rand_mat));
end

cef = [];
rand = [];
for j = check_inds
    [cef_mat, rand_mat] = cal_stability(output_file_name, j, 0, num_fold);
    % 将cef_mat填充为max_len长度
    padded_cef_mat = padarray(cef_mat(:), max_len_cef - numel(cef_mat), NaN, 'post');
    padded_rand_mat = padarray(rand_mat(:), max_len_rand - numel(rand_mat), NaN, 'post');
    cef = horzcat(cef, padded_cef_mat);
    rand = horzcat(rand, padded_rand_mat);
    rs(j) = mean(cef_mat(:),'omitnan');
    ris(j) = mean(rand_mat(:),'omitnan');
end
f = figure;
f.Position = get_figure_position(2);
subplot(1,2,1)
hold on;
% 绘制箱线图，check_inds将作为分类标签
%     boxplot(cef, check_inds);
x = repmat(check_inds,size(cef,1),1);
swarmchart(x, cef, 10, "red", 'filled', 'MarkerFaceAlpha', 0.4);
boxchart(reshape(x,[],1), cef(:), 'markerstyle', 'none');

xlabel('Number of subtypes')
ylabel('Trajectory similarity between training sets');
xticks(1:7)

subplot(1,2,2)
hold on;
%boxplot(rand, check_inds);
x = repmat(check_inds,size(rand,1),1);
swarmchart(x, rand, 10, "red", 'filled', 'MarkerFaceAlpha', 0.4);
boxchart(reshape(x,[],1), rand(:), 'markerstyle', 'none');

xlabel('Number of subtypes')
ylabel('Adjusted Rand index between training sets');
xticks(1:7)

export_fig './figures/all/S11B.jpg' -r500 -transparent
% close all

rs(1) = nan;
ris(1) = nan;

f = figure;
f.Position = get_figure_position();
plot(rs)
xlabel('Number of subtypes');
ylabel('Mean trajectory similarity between training sets');
export_fig './figures/all/S11C.jpg' -r500 -transparent
% close all
f = figure;
f.Position = get_figure_position();
plot(ris)
xlabel('Number of subtypes');
ylabel('Mean adjusted Rand index between training sets');
export_fig './figures/all/S11D.jpg' -r500 -transparent
% close all
end

