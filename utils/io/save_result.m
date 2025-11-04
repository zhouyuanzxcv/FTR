function save_result(all_data, stage_all, subtype_all, mdl, ...
    biomarker_names, thresh, extra, save_dir1)

model_name = '';
if isfield(mdl, 're_traj')
    model_name = 'FTR';
elseif isfield(mdl, 'samples_sequence')
    model_name = 'sustain';
end

if strcmp(model_name, 'FTR')
    num_subtype = size(mdl.re_traj, 3);
elseif strcmp(model_name, 'sustain')
    num_subtype = size(mdl.samples_sequence, 1);
end

writetable(array2table([all_data.RID, all_data.years, all_data.labels, stage_all, subtype_all], ...
    'VariableNames', {'PTID' ,'years_from_baseline','Diagnosis','stage','subtype'}), ...
    [save_dir1,'/subtype_stage.csv'])

if strcmp(model_name, 'FTR')
    save_parameters_theta(save_dir1, mdl, biomarker_names);
else
    save([save_dir1,'/mdl.mat'], 'mdl', 'biomarker_names');
end

if ~isempty(extra)    
    if isfield(extra, 'subtype_prob')
        for k = 1:num_subtype
            subtype_names{k} = ['subtype', int2str(k)];
        end

        for k = 1:num_subtype
            subtype_names_visit{k} = ['subtype', int2str(k),'_visit'];
        end

        subtype_prob = array2table([all_data.RID, all_data.years, ...
            all_data.labels, extra.subtype_prob, extra.subtype_prob_per_visit], ...
            'VariableNames', [{'PTID' ,'years_from_baseline','Diagnosis'},...
            subtype_names, subtype_names_visit]);
        writetable(subtype_prob, [save_dir1,'/subtype_prob.csv']);
    end

    if isfield(extra, 'stage_prob')
        % do not store stage probability for now
    end
    
    if isfield(extra, 'elapsed_t') && ~extra.re_subtype_staging
        writematrix(extra.elapsed_t, [save_dir1,'/elapsed_t.csv']);
    end

%     extra.stage_prob;

    if isfield(extra, 'loglik_all')
        writematrix(extra.loglik_all, [save_dir1,'/loglik_all.csv']);
    end

    if isfield(extra, 'init_runs') && ~isempty(extra.init_runs)
        init_runs = extra.init_runs;
        save([save_dir1,'/init_runs.mat'], 'init_runs');
    end
end

if strcmp(model_name, 'FTR')
    ord = cal_order(mdl.re_traj,thresh,biomarker_names,save_dir1);
end

end

