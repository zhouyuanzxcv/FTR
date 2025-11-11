function [outputArg1,outputArg2] = analyze_stages_4E(inputArg1,inputArg2)
%ANALYZE_STAGES Summary of this function goes here
%   Detailed explanation goes here

close all

data_names = {'ADNI_FSX_HS','OASIS3_ROD1_HS','NACC_HS'};
method_names = {'FTR_MCEM'};


nsubtype = 3;
data_sel = 2;

include_control = 1;
results = load_data_result(data_names, method_names, nsubtype, include_control, 1);
joindata = results{data_sel,1}.joindata;
joindata.Properties.VariableNames = strrep(joindata.Properties.VariableNames, '-', '_');


data = joindata;

if data_sel == 3
    data.MMSE = -data.MMSE;
end

%% calculate r (vs. MMSE) and p
f = figure;
f.Position = get_figure_position(5);

plot_stage_vs_MMSE(data, nsubtype);

if data_sel == 3
    ylim([0,18]);
    set(gca, 'YDir', 'reverse');
    ylabel('CDRSB');
    % yticklabels(arrayfun(@(x) num2str(-x), yticks, 'UniformOutput',false));
end

if data_sel == 2
    export_fig './figures/all/4E1.jpg' -r500 -transparent
elseif data_sel == 3
    export_fig './figures/all/4E2.jpg' -r500 -transparent
end
% legend(hs, 'Subtype 1', 'Subtype 2', 'Subtype 3')



end


