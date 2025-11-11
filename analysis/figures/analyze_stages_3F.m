function [outputArg1,outputArg2] = analyze_stages(inputArg1,inputArg2)
%ANALYZE_STAGES Summary of this function goes here
%   Detailed explanation goes here

close all

data_names = {'ADNI_FSX_HS'};
method_names = {'FTR_MCEM'};

nsubtype = 3;

% 1 for discovery set, 20 for validation set
include_control = 20;
results = load_data_result(data_names, method_names, nsubtype, include_control);
joindata = results{1,1}.joindata;
joindata.Properties.VariableNames = strrep(joindata.Properties.VariableNames, '-', '_');


stage = joindata.stage;
data = joindata;

%% calculate AUC

AUCs = calc_AUCs(data);

disp(['AUC for CN vs. AD (row 1), MCI vs. AD (row 2), CN vs. MCI (row 3). ', ...
    'The columns are mean (col. 1) and 95% CI (col. 2 and 3)']);
disp(AUCs);

% show cross validation
if 0
biomarker_names = results{1,1}.biomarker_names;
control_data = data(data.group == 0, :);
id_inds = crossvalind('Kfold', size(control_data,1), 10);
options = [];
options.method = 'MCEM';
output_file_name = [data_names{1}, '_', method_names{1}];

for i = 1:10

    save_dir_nsubtype1 = ['output/', output_file_name, ... 
        '/cross_validation_nsubtype3_fold', int2str(i)];

    test_dir = [save_dir_nsubtype1, '/test_subtype_stage.csv'];

    test_subtype_stage = readtable(test_dir);
    test_PTID = unique(test_subtype_stage.PTID);

    test_data = data(ismember(data.RID, test_PTID),:);
    
    input_data = [control_data(id_inds == i, :); test_data];
    
    % [sigma, re_traj, proption]
    mdl = load_parameters_theta(save_dir_nsubtype1);
    
    stage_all = cal_stage_subtype(input_data{:, biomarker_names}, input_data.RID, mdl, options);

    pred_label_sel = stage_all(input_data.diagnosis~=0.5);

    true_label_sel = input_data.diagnosis(input_data.diagnosis~=0.5);
    
    [X,Y,T,~] = perfcurve(true_label_sel, pred_label_sel, 1);
    
    plot(X,Y,':','Color',[0,0,0],'LineWidth',0.5)
end

end

if include_control == 1
    export_fig './figures/all/3F.jpg' -r500 -transparent
elseif include_control == 20
    export_fig './figures/all/3F1.jpg' -r500 -transparent
end

end
