function [outputArg1,outputArg2] = compare_subtype_progression_3B(input,subtype_stage,subtype_name)

close all

data = {'ADNI_FSX_HS'};

method = {'FTR_MCEM'}; % FTR_MCEM, sustain


nsubtype = 3;
data_sel = 1;

ref_group = 3;

% 0 for discovery set, 2 for validation set
include_control = 0;

results = load_data_result(data, method, nsubtype, include_control, 1);
joindata = results{data_sel,1}.joindata;
mmse = results{data_sel,1}.mmse;

joindata.Properties.VariableNames = strrep(joindata.Properties.VariableNames, '-', '_');


%% Fit linear mixed model
% Define the variables

variables = {'MMSE'};


for varIdx = 1:length(variables)
    dependentVariable = variables(varIdx); % Dependent variable (MMSE measurements)
%     independentVariables = {'RID','years', 'subtype', 'AGE', 'PTGENDER', 'PTEDUCAT', 'diagnosis', 'stage'}; 
%     age_bl = joindata.AGE - joindata.years;
%     joindata = addvars(joindata, age_bl, 'After','AGE');
    
    covs = {'AGE_baseline', 'PTGENDER', 'PTEDUCAT', 'diagnosis_bl'};
    independentVariables = [{'RID','years', 'subtype'}, covs]; 
    joindata = joindata(:,[{'RID','subtype'},covs]);
    joindata = unique(joindata);
    mmse1 = mmse(:,['RID','years',dependentVariable]);
    joindata1 = outerjoin(joindata, mmse1,'Type','Left','Keys','RID','MergeKeys',true);
    joindata = joindata1;

    % dx_bl1 = get_baseline_diagnosis(joindata);
    dx_bl = joindata.diagnosis_bl;
    
    data_sel_inds = {dx_bl == 0.5 | dx_bl == 1, dx_bl == 0.5, dx_bl == 1};
    data_sel_names = {'MCI/AD at BL', 'MCI at BL', 'AD at BL'};
        
    slopes_all_stages = compare_lmem(joindata, dependentVariable, ...
        independentVariables, nsubtype, data_sel_inds, ref_group);
    
    %% plot progression rates distribution
    f = figure;
    f.Position = get_figure_position();
    hs = boxplot_by_subtype_stage(slopes_all_stages, data_sel_names, 'MMSE decrease per year');
    legend(hs, {'Subtype 1', 'Subtype 2', 'Subtype 3'}, "Location", "best")
    ylim([-10,2]);

    if include_control == 0        
        export_fig './figures/all/3B.jpg' -r500 -transparent
    else
        export_fig './figures/all/3B1.jpg' -r500 -transparent
    end
    
end

end


function slopes_all_stages = compare_lmem(joindata, dependentVariable, ...
    independentVariables, nsubtype, data_sel_inds, ref_group)


num_data = length(data_sel_inds);


slopes_all_stages = cell(num_data, nsubtype);

for j = 1:num_data
    data = joindata(data_sel_inds{j}, [dependentVariable, independentVariables]);
    
    data = rmmissing(data);
    

    [slopes_all, beta_all, intercept_all, beta_p_all] = lmem(data, ...
        dependentVariable, nsubtype, ref_group);
    
    slopes_all_stages(j,:) = slopes_all;
    

end
end




