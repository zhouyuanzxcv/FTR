function [sim_mats,elapsed_ts] = compare_algorithmic_stability_1G(dataset_name)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin < 1
    dataset_name = {'ADNI_FSX_LM', 'OASIS3_ROD1_LM', 'NACC_LM', ...
        'ADNI_FSX_HS', 'OASIS3_ROD1_HS', 'NACC_HS'};
end

method_name = {'FTR_MCEM','sustain'};

[elapsed_ts, sim_mats] = get_time_stability(dataset_name, method_name);

if 1
fh = figure;
fh.Position = get_figure_position(9);

dataset_name1 = {'ADNI-13', 'OASIS-13', 'NACC-7', ...
    'ADNI-82', 'OASIS-82', 'NACC-64'};
method_name1 = {'FTR','SuStaIn'};
crossed_boxplot(sim_mats, elapsed_ts, dataset_name1, method_name1);

export_fig './figures/all/1G.jpg' -r500 -transparent
end

end



