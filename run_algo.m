function [subtype_stage, mdl, options, train_data, test_data] = ...
    run_algo(dataset_name, method_name, options)
%RUN_ALGO Summary of this function goes here
%   Detailed explanation goes here


addpath(genpath("FTR_code/"));
addpath(genpath("utils/"));
addpath(genpath('analysis'));
addpath(genpath('competing_methods'));
addpath('figures');



if nargin < 3
    options = [];
end

if nargin < 2
    method_name = 'FTR_MCEM';
end

if nargin < 1
    dataset_name = 'ADNI_FSX_HS';
end

% if validate on other datasets, add rows of other data to the current data
% in csv
re_subtype_staging = parse_param(options, 're_subtype_staging', 0);

thresh_order = parse_param(options, 'thresh_order', 1);

% method_name1 = method_name;

switch method_name
    case {'FTR','FTR_kmeans','FTR_MCEM'}     
        % if nsubtype is a list (e.g. (1:6)), run cross validation
        options = insert_param_when_absent(options, 'nsubtype', 3);

        % number of folds in cross validation
        options = insert_param_when_absent(options, 'num_folds', 10);
        
        % number of samples (L) to approximate the Q function in MCEM
        options = insert_param_when_absent(options, 'samp_multi', 100);

        % number of iterations in MCEM
        options = insert_param_when_absent(options, 'max_iter', 30);

        % set number of initializations
        options = insert_param_when_absent(options, 'max_ep', 100);
        options = insert_param_when_absent(options, 'max_tries', 10);
        options = insert_param_when_absent(options, 'num_same_partition', 20);

        % set the method to select from multiple random initializations,
        % 'meta_clustering' or 'partition_mode'
        options = insert_param_when_absent(options, 'partition_selection_method', ...
            'meta_clustering');

        % set the method to initialize the partition (z)
        % 'kmeans++' or 'random_permutation' or 'randi'
        options = insert_param_when_absent(options, 'init_method', 'kmeans++');

%         options.num_int = 101;
        
        % set the method to perform deconvolution in FTR
        % 'search_max', 'min_max', 'approx'
        options = insert_param_when_absent(options, 'deconv', 'search_max'); 
%         options = insert_param_when_absent(options, 'biom_max', 5);
%         options = insert_param_when_absent(options, 'stage_thresh', 1);
        
        if length(method_name) > 3
            if strcmp(method_name, 'FTR_kmeans')
                options.methods = 'kmeans';
            elseif strcmp(method_name, 'FTR_MCEM')
                options.methods = 'MCEM';
            end
            method_name = 'FTR';
        else
            % This is the method proposed in the paper.
%             options.methods = 'kmeans';
            options.methods = 'MCEM';
        end
        
        % sigma_type can be '1x1' or '1xK' or 'BxK'
        options = insert_param_when_absent(options, 'sigma_type', '1xK'); 

        % for zscore data, initial sigma is set to 1
        options = insert_param_when_absent(options, 'sigma_init', 1); 
        
        % This parfor is for multiple random intializations given a 
        % single number of subtypes
        options = insert_param_when_absent(options, 'parfor', true);
        
        method_name1 = [method_name,'_',options.methods];
        
        % whether to first remove the data points whose stages are 0 in the 
        % single trajectory setting
        options.stagsel = 1;
        
        default_output = sprintf('%s_%s', dataset_name, method_name1);        
%         default_output = sprintf('%s_%s_%s_%s_%s', dataset_name, method_name1, ...
%             options.sigma_type, options.deconv, options.select_from_multi_inits);    

        options = insert_param_when_absent(options, 'output_file_name', default_output);

    case 'sustain'  
        options = insert_param_when_absent(options, 'nsubtype', 3);
        
        output_folder = sprintf('%s_%s', dataset_name, method_name);
        options = insert_param_when_absent(options, 'output_file_name', output_folder);
        
        % use the default setting in SuStaIn
        options.N_startpoints = 25;
        options.N_iterations_MCMC = int32(1e5);
        options = insert_param_when_absent(options, 'Z_vals', [1,2,3]);
        options = insert_param_when_absent(options, 'Z_max', 5);
        options = insert_param_when_absent(options, 'std_biomarker_zscore', 1);
    case {'hierarch_clustering_baseline'}
        options = insert_param_when_absent(options, 'nsubtype', 3);
        if strcmp(method_name(end-7:end),'baseline')
            options.data_split = 'baseline';
        elseif strcmp(method_name(end-7:end),'first_AD')
            options.data_split = 'first_AD';
        end

        method_name1 = method_name;
        method_name = 'hierarch_clustering';
        default_output = sprintf('%s_%s', dataset_name, method_name1);        
        options = insert_param_when_absent(options, 'output_file_name', default_output);

    case {'ctv_hv_ratio_baseline'}
        options = insert_param_when_absent(options, 'nsubtype', 3);
        if strcmp(method_name(end-7:end),'baseline')
            options.data_split = 'baseline';
        elseif strcmp(method_name(end-7:end),'first_AD')
            options.data_split = 'first_AD';
        end
        options.biomarker_type = 'ctv_hv';

        method_name1 = method_name;
        method_name = 'ctv_hv_ratio';
        default_output = sprintf('%s_%s', dataset_name, method_name1);        
        options = insert_param_when_absent(options, 'output_file_name', default_output);
    case {'debm','ebm'}
        default_output = sprintf('%s_%s', dataset_name, method_name);
        options = insert_param_when_absent(options, 'output_file_name', default_output);

        options.group_sel = 10;

    otherwise
end

   
options.dataset_name = dataset_name;

switch dataset_name
    case {'ADNI_FSX_HM', 'ADNI_FSX_HS', 'ADNI_FSX_LM', 'ADNI_FSX_LS'}
        fs_ver = dataset_name(6:8);
        default_dx_sel = false;
        fldstr = '';

        if length(dataset_name) == 11 % 'ADNI_FSX_HM'
            biom_ver = dataset_name(10:11);
        elseif length(dataset_name) == 15 % 'ADNI_FSX_1.5_HM'
            biom_ver = dataset_name(end-1:end);
            fldstr = '1.5_';
        elseif length(dataset_name) == 13 % 'ADNI_FSX_3_HM'
            biom_ver = dataset_name(end-1:end);
            fldstr = '3_';
        elseif length(dataset_name) == 14 % 'ADNI_FSX_HM_DS'
            biom_ver = dataset_name(end-4:end-3);
            default_dx_sel = true;
        end

        options.input_file_name = ['ADNI/ADNI_demo_MRI_',fldstr,fs_ver, '_PET_', ...
            biom_ver,'_zscore.csv'];

        num_bioms = get_num_bioms_from_ver(biom_ver, 'ADNI');
        options = insert_param_when_absent(options, 'group_sel', true);
        options = insert_param_when_absent(options, 'diagnosis_sel', default_dx_sel);
        options = insert_param_when_absent(options, ...
            'biomarker_column_range', [3, 3 + num_bioms - 1]);

    case {'OASIS3_ROD1_HM','OASIS3_ROD1_HS','OASIS3_ROD1_LM','OASIS3_ROD1_LS'}
        oasis_ver = dataset_name(end-1:end);
        oasis_rod = dataset_name(end-3);
        options.input_file_name = ['OASIS3/OASIS3_input_ROD',oasis_rod,'_',oasis_ver,'_z'];
        num_bioms = get_num_bioms_from_ver(oasis_ver, 'OASIS');
        options = insert_param_when_absent(options, 'group_sel', true);
        options = insert_param_when_absent(options, ...
            'biomarker_column_range', [3,3 + num_bioms - 1]);

    case {'NACC_HM','NACC_HS','NACC_LM','NACC_LS'}
        nacc_ver = dataset_name(end-1:end);
        
        options.input_file_name = ['NACC/NACC_',nacc_ver,'_z.csv'];
        num_bioms = get_num_bioms_from_ver(nacc_ver, 'NACC');
        options = insert_param_when_absent(options, 'group_sel', true);
        options = insert_param_when_absent(options, ...
            'biomarker_column_range', [3,3 + num_bioms - 1]);

    case 'custom'
        options.group_sel = false;
        options.diagnosis_sel = false;
        options.biomarker_column_range = [4,0];
    otherwise
end

options = insert_param_when_absent(options, 'biomarker_column_range', [4,0]);
options = insert_param_when_absent(options, 'diagnosis_sel', false);
options = insert_param_when_absent(options, 'group_sel', false);

% if options.diagnosis_sel
%     options.output_file_name = [options.output_file_name, '_dxsel=1'];
% end

%% Read data

train_inds = parse_param(options,'train_inds',[]);
test_inds = parse_param(options,'test_inds',[]);

[train_data,test_data,all_data,biomarker_name,options] = load_dataset( ...
    dataset_name, train_inds, test_inds, options);
     
%% Run FTR
start_t = tic;

nsubtype = options.nsubtype;

[save_dir1, save_dir] = get_save_result_filepath(options, nsubtype(1));
save([save_dir, '/all_data.mat'], 'biomarker_name');


switch method_name
    case 'FTR'
        if options.stagsel && ~re_subtype_staging
            [train_vol_sel,train_PTID_sel,stage_sel,inds_sel] = sel_stage(...
                train_data.vols,train_data.RID,options);
            
            fprintf('Removing stage 0 points from %d points leads to %d remaining.\n', ...
                size(train_data.vols,1), size(train_vol_sel,1));

            options.inds_sel = inds_sel;
            
            % for debug
%             options.subtypes_gt = options.subtypes_gt(inds_sel);
%             options.stages_gt = options.stages_gt(inds_sel);
        else
            train_vol_sel = train_data.vols;
            train_PTID_sel = train_data.RID;            
        end

        if length(nsubtype) > 1 && ~re_subtype_staging
            % cross validation
%             options.parfor = false;

            FTR_model_selection(train_vol_sel,train_PTID_sel,options);
            
            for nsubtype1 = options.nsubtype
                save_dir_nsubtype1 = ['output/', options.output_file_name, ... 
                    '/cross_validation_nsubtype', int2str(nsubtype1), '_fold0'];
                
                % [sigma, re_traj, proption] 
                mdl = load_parameters_theta(save_dir_nsubtype1);
                
                [stage_all,subtype_all,extra] = cal_stage_subtype(all_data.vols, ...
                    all_data.RID, mdl, options);
                
                [save_dir1, save_dir] = get_save_result_filepath(options, nsubtype1);

                save_result(all_data, stage_all, subtype_all, ...
                    mdl, biomarker_name, thresh_order, extra, save_dir1);
            end   
            
            return;
        else
            if re_subtype_staging
                mdl = load_parameters_theta(save_dir1);
                init_runs = [];
            else
                [mdl,pre_subtype,stage,extra] = FTR_model( ...
                    train_vol_sel,train_PTID_sel,nsubtype,options);
                init_runs = extra.init_runs;
            end
            
            % predict test
            [stage_all,subtype_all,extra] = cal_stage_subtype(all_data.vols, ...
                all_data.RID, mdl, options);
            extra.init_runs = init_runs;
        end
    case 'sustain'
        if re_subtype_staging
            S = load([save_dir1, '/mdl.mat']);
            mdl = S.mdl;
        else
            mdl = sustain_wrapper(train_data.vols, options.nsubtype, options);
            %         mdl = load('tpc7be50d5_7866_47ec_8497_868921b9d255_result.mat');
        end

        % In sustain, an average of the probabilities of all the time
        % points of a subject should be calculated to derive the subtype
        % label
        [stage_all,subtype_all,extra] = cal_stage_subtype_sustain( ...
            all_data.vols, all_data.RID, mdl, options);
    case 'hierarch_clustering'
        mdl = hierarchical_clustering(train_data.vols, train_data.RID, options.nsubtype, options);
        stage_all = NaN(size(all_data.RID));
        subtype_all = NaN(size(all_data.RID));
        for k = 1:options.nsubtype
            ind = ismember(all_data.RID, mdl.subtypes{k});
            subtype_all(ind) = k;
        end
        extra = [];

    case 'ctv_hv_ratio'
        mdl = hv_ctv_ratio(train_data.vols, train_data.RID);
        stage_all = NaN(size(all_data.RID));
        subtype_all = NaN(size(all_data.RID));
        for k = 1:options.nsubtype
            ind = ismember(all_data.RID, mdl.subtypes{k});
            subtype_all(ind) = k;
        end
        extra = [];
    case 'debm'
        mdl = ebm(train_data.vols, train_data.RID, train_data.labels, ...
            train_data.years, biomarker_name, 'debm');
        stage_all = NaN(size(all_data.RID));
        subtype_all = NaN(size(all_data.RID));
        stage_all(options.train_inds) = mdl.stage;
        extra = [];
    case 'ebm'
        mdl = ebm(train_data.vols, train_data.RID, train_data.labels, ...
            train_data.years, biomarker_name, 'ebm');
        stage_all = NaN(size(all_data.RID));
        subtype_all = NaN(size(all_data.RID));
        stage_all(options.train_inds) = mdl.stage;
        extra = [];

    otherwise
end

elapsed_t = toc(start_t)/60;
fprintf('The method takes %.1f minutes\n', elapsed_t);
extra.elapsed_t = elapsed_t;
extra.re_subtype_staging = re_subtype_staging;

subtype_stage = array2table( ...
    [all_data.RID,all_data.years,all_data.labels,stage_all,subtype_all], ...
    'VariableNames', {'PTID' ,'years_from_baseline','Diagnosis','stage','subtype'});

% %% permute the subtypes
% [data,mmse,time2event,dx_sel] = load_data_all(dataset_name);
% P = determine_trajectory_order(joindata, mmse, nsubtype);
% [mdl, joindata] = permute_subtypes(mdl, joindata, P);

%% save output

save_results = parse_param(options, 'save_results', 1);

if ~isempty(options.output_file_name) && save_results
    [save_dir1, save_dir] = get_save_result_filepath(options, nsubtype);
    save([save_dir, '/all_data.mat'], 'biomarker_name');
    
    save_result(all_data, stage_all, subtype_all, mdl, ...
        biomarker_name, thresh_order, extra, save_dir1);
end


end

function [save_dir1, save_dir] = get_save_result_filepath(options, nsubtype)
postfix = parse_param(options, 'postfix', '');
save_dir = ['./output/',options.output_file_name];
save_dir1 = [save_dir, '/',num2str(nsubtype),'_subtypes',postfix];
if ~exist(save_dir1, 'dir')
    mkdir(save_dir1)
end

end


