function compare_subtype_abeta_tau_LR()
close all

data_names = {'ADNI_FSX_HS','OASIS3_ROD1_HS','NACC_HS'};

method = {'FTR_MCEM'};

nsubtype = 3;
data_sel = 1;

% include_control = 0;
results = load_data_result(data_names, method, nsubtype, 0);
joindata = results{data_sel,1}.joindata;
traj = results{data_sel,1}.mdl.re_traj;
biomarker_names = [results{data_sel,1}.biomarker_names];

joindata.Properties.VariableNames = strrep(joindata.Properties.VariableNames, '-', '_');
biomarker_names = strrep(biomarker_names, '-', '_');

Abeta = 0;
Tau = 0;

if Abeta == 1 % abeta
    biomarker_names_pet = cellfun(@(x) ['Abeta_',x], biomarker_names, 'UniformOutput', 0);
elseif Tau == 1 % tau
    biomarker_names_pet = cellfun(@(x) ['Tau_',x], biomarker_names, 'UniformOutput', 0);
else % for MRI
    biomarker_names_pet = biomarker_names;
end

% remove nan entries in the pet data
nan_inds = any(isnan(joindata{:,biomarker_names_pet}), 2);
data = joindata;
data(nan_inds, :) = [];

fprintf('%d subjects, %d points\n', length(unique(data.RID)), size(data,1));

%% compare one subtype with all the other subtypes in terms of mean pattern difference

data_for_mean = data;


disp_cells = compare_mean_pattern(data_for_mean, nsubtype, biomarker_names_pet);

for i = 1:size(disp_cells, 1)
    for j = 1:size(disp_cells, 2)
        if isempty(disp_cells{i,j})
            disp_cells{i,j} = 0;
        end
    end
end

disp_cells

image_t = cell2mat(disp_cells(:,2:end));
M = max(abs(image_t), [], "all");
image_t = image_t' / M + 1;

image_t = array2table(image_t);
image_t.Properties.VariableNames = cellfun(@(x) strrep(x, '_', '-'), disp_cells(:,1), 'UniformOutput', false);
names = arrayfun(@(x) sprintf('Image%d', x), 1:3, 'UniformOutput', false);
image_t = addvars(image_t, names', 'Before', 1, 'NewVariableNames',{ 'Image-name-unique' });

writetable(image_t, ['analysis/brainpainting/input/4B_', data_names{data_sel}, '.csv']);



end




