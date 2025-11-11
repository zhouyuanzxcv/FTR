function [outputArg1,outputArg2] = analyze_stages_3C(inputArg1,inputArg2)
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

data = joindata;


%% calculate r (vs. MMSE) and p
f = figure;
f.Position = get_figure_position();

plot_stage_vs_MMSE(data, nsubtype);

export_fig './figures/all/3C.jpg' -r500 -transparent
% legend(hs, 'Subtype 1', 'Subtype 2', 'Subtype 3')



end

