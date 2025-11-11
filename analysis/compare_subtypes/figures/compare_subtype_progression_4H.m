function [outputArg1,outputArg2] = compare_subtype_progression_4H(input,subtype_stage,subtype_name)

close all

data = {'ADNI_FSX_HS', 'ADNI_FSX_HM', 'ADNI_FSX_LS', 'ADNI_FSX_LM'};

method = {'FTR_MCEM'}; % FTR_MCEM, sustain

nsubtype = 3;
results = load_data_result(data, method, nsubtype, 0);
ref_group = 1;

%% Fit linear mixed model
% Define the variables
r = {};

f = figure;
f.Position = get_figure_position();

ps = [];
for data_sel = 2:4

%     ref_data = results{1,1}.joindata;
    
%     results = load_data_result(data, method, nsubtype, 0, 1);
    joindata = results{data_sel,1}.joindata;
    mmse = results{data_sel,1}.mmse;

    % subtype1 = dummyvar(ref_data.subtype);
    % subtype2 = dummyvar(joindata.subtype);
    % 
    % [P, subtype2] = permute_endmembers(subtype1', subtype2');
    % 
    % [~, subtype2] = max(subtype2', [], 2);
    % 
    % joindata.subtype = subtype2;
    
    joindata.Properties.VariableNames = strrep(joindata.Properties.VariableNames, '-', '_');
    dependentVariable = {'MMSE'}; % Dependent variable (MMSE measurements)

    covs = {'AGE_baseline', 'PTGENDER', 'PTEDUCAT', 'diagnosis_bl'};
    independentVariables = [{'RID','years','subtype'}, covs]; 
    joindata = joindata(:,[{'RID','subtype'},covs]);
    joindata = unique(joindata);
    mmse1 = mmse(:,['RID','years',dependentVariable]);
    joindata1 = outerjoin(joindata, mmse1,'Type','Left','Keys','RID','MergeKeys',true);
    joindata = joindata1;

    % dx_bl1 = get_baseline_diagnosis(joindata);
    dx_bl = joindata.diagnosis_bl;
    
    data_sel_inds = {dx_bl == 0.5 | dx_bl == 1};
        
    [slopes_all_stages, beta_p_all] = compare_lmem(joindata, dependentVariable, ...
        independentVariables, nsubtype, data_sel_inds, ref_group);

    ps(end+1,:) = beta_p_all;

    r = [r; slopes_all_stages];

end

disp('ps for difference (row: data, col: subtype)');
disp(ps);

atlases = {'41','26','13'};
boxplot_by_subtype_stage(r, atlases, 'MMSE decrease per year');
ylim([-10,2]);
xlabel('Number of regions')
export_fig './figures/all/4H.jpg' -r500 -transparent

end


function [slopes_all_stages, beta_p_all] = compare_lmem(joindata, dependentVariable, ...
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


function boxplot_by_subtype_stage(slopes_all,stage_names,yname)

nsubtype = size(slopes_all,2);

subtype_offsets = [-floor(nsubtype/2):floor(nsubtype/2)];
subtype_name = cellfun(@(x) ['Subtype ', num2str(x)], num2cell(1:nsubtype), 'UniformOutput', 0);

hold on;

for j = 1:length(stage_names)
    for k = 1:nsubtype
        y = slopes_all{j,k};
        
        stage_width = (nsubtype+1);
        x = (j-1)*stage_width + 1 + subtype_offsets(k);
        x = repmat(x, size(y,1), 1);
        
        swarmchart(x, y, 5, get_subtype_color(k), 'filled', 'MarkerFaceAlpha', 0.3);
        hs(k) = boxchart(x, y, 'boxfacecolor', ...
            get_subtype_color(k)/2, 'whiskerlinecolor', get_subtype_color(k)/2, ...
            'markerstyle', 'none');
        
        
    end
end

xticks((0:nsubtype-1)*stage_width + 1);
xticklabels(stage_names);
ylabel(yname);
end

