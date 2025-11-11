function [outputArg1,outputArg2] = analyze_stages_3E(inputArg1,inputArg2)
%ANALYZE_STAGES Summary of this function goes here
%   Detailed explanation goes here

close all

data_names = {'ADNI_FSX_HS'};
method_names = {'FTR_MCEM'};


nsubtype = 3;

include_control = 1;
results = load_data_result(data_names, method_names, nsubtype, include_control, 1);
joindata = results{1,1}.joindata;
joindata.Properties.VariableNames = strrep(joindata.Properties.VariableNames, '-', '_');


stage = joindata.stage;
data = joindata;

%% Plot stage

f = figure;
f.Position = get_figure_position();

plot_stage_distribution_by_diagnosis(data);

export_fig './figures/all/3E.jpg' -r500 -transparent

end


