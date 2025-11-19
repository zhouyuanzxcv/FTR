function [outputArg1,outputArg2] = compare_subtype_progression_3A(input,subtype_stage,subtype_name)

close all

data = {'ADNI_FSX_HS'};

method = {'FTR_MCEM'}; % FTR_MCEM, sustain


nsubtype = 3;
data_sel = 1;

ref_group = 1;

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
    
    data_sel_inds = {dx_bl == 0.5 | dx_bl == 1};
    data_sel_names = {'MCI/AD at baseline'};
        
    slopes_all_stages = compare_lmem(joindata, dependentVariable, ...
        independentVariables, nsubtype, data_sel_inds, ref_group);
    
end

end


function slopes_all_stages = compare_lmem(joindata, dependentVariable, ...
    independentVariables, nsubtype, data_sel_inds, ref_group)


num_data = length(data_sel_inds);

f = figure;
f.Position = get_figure_position();

slopes_all_stages = cell(num_data, nsubtype);

for j = 1:num_data
    data = joindata(data_sel_inds{j}, [dependentVariable, independentVariables]);
    
    data_original = data;
    data = rmmissing(data);
    
    subtype = data.subtype;
    years = data.years;
    

    [slopes_all, beta_all, intercept_all, beta_p_all] = lmem(data, ...
        dependentVariable, nsubtype, ref_group);
    
    slopes_all_stages(j,:) = slopes_all;

    subplot(num_data, 1, j);
    hold on;

    max_years = max(years);

    for k = nsubtype:-1:1
        data_k = data(subtype == k, :);
        years_k = years(subtype == k);
        
        RID_k = data_k{:,'RID'};
        MMSE_k = data_k{:,dependentVariable};
        
        MMSE_k1 = MMSE_k + rand(size(MMSE_k)) - 0.5;
        
        %% plot progression line
        
        markerSize = 20;
        dot_color = [0,0,0];
        
        [unique_RID, ia, ic] = unique(RID_k);
        
        for each_RID = unique_RID'
            inds = each_RID == RID_k;
            x = years_k(inds)+rand(size(years_k(inds)));
            y = MMSE_k1(inds);
            s = scatter(x, y, markerSize, get_subtype_color(k), "filled", "MarkerFaceAlpha", 0.25);
            % s.AlphaData = sqrt(((x+8)/16*30.5).^2 + (y-45).^2);
            % s.MarkerFaceAlpha = 'flat';
        end
        
        % scatter(years_k,  MMSE_k1, markerSize, ...
        %     'MarkerFaceColor', dot_color, 'MarkerEdgeColor', 'none');
        xlabel('Years from baseline');
        ylabel(dependentVariable);
        
        beta_k = beta_all(k);
        intercept_k = intercept_all(k);
        beta_p_k = beta_p_all(k);
        
        line_x = linspace(0, max_years, 101)';
        plot(line_x, line_x * beta_k + intercept_k, 'color', get_subtype_color(k), 'linewidth', 2);
        
        % legend(['Slope: ',num2str(beta_k), ', p value: ',num2str(beta_p_k)]);
             
    end
    xlim([0 max_years])
    ylim([0 30.5])
    export_fig './figures/all/3A.jpg' -r500 -transparent
end
    
% slopes_all_stages = slopes_all_stages;

end

