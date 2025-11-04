function AUCs = calc_AUCs(data, nboot)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%   AUCs - 3 x 3 matrix. Each row is a classification task (e.g. CN vs. AD)
%           The first column is the mean AUC. The second column and the
%           third column is the 95% confidence interval of AUC from 1000
%           bootstraps.
if nargin < 2
    nboot = 1000;
end

dx_sel = [0.5,0,1];
pos_class = [1,1,0.5];
class_labels = {'CN vs. AD', 'MCI vs. AD', 'CN vs. MCI'};

stage = data.stage;

f = figure;
f.Position = get_figure_position();
hs = [];
legend_labels = {};

colors = [0.8 0.2 0.2; 0.2 0.8 0.2; 0.2 0.2 0.8;];


AUCs = [];
for i = 1:length(dx_sel)
    pred_stage_sel = stage(data.diagnosis~=dx_sel(i));
    true_label_sel = data.diagnosis(data.diagnosis~=dx_sel(i));
    [X,Y,T,AUC, OPTROCPT] = perfcurve(true_label_sel,pred_stage_sel, pos_class(i), ...
        'NBoot', nboot);
    AUCs(i,:) = AUC;


    % Extract confidence intervals (stored in X and Y)
    fpr = X(:,1); % Mean FPR
    tpr_mean = Y(:,1); % Mean TPR

    if nboot > 0 
        tpr_lower = Y(:,2); % Lower bound of TPR (95% CI)
        tpr_upper = Y(:,3); % Upper bound of TPR (95% CI)
    end

    h = plot(fpr, tpr_mean, '-', 'color',colors(i,:),'LineWidth', 2); % Mean ROC
    hold on;

    if nboot > 0
        fill([fpr; flipud(fpr)], [tpr_lower; flipud(tpr_upper)], ...
            colors(i,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Confidence band
    end

    hs(i) = h;

    legend_labels{i} = ['AUC (',class_labels{i},') = ',num2str(AUC(1),2) ...
        ];
    
end

xlabel('False positive rate')
ylabel('True positive rate')
legend(hs, legend_labels, 'Location', 'southeast');
title('ROC curve')
end

