function results = load_data_result(datasets, methods, nsubtypes, ...
    include_control, options)
%LOAD_DATA_AND_RESULT Summary of this function goes here
%   Detailed explanation goes here
% Input:
%   include_control - 0 for only including the case group, 
%                     1 for including both the case and control group,
%                     2 for including only the validation group,
%                     20 for including both the validation and control
%                        group

if nargin < 4
    include_control = 0;
end

if nargin < 5
    options = [];
end

postfix = parse_param(options, 'postfix', ''); % e.g. _run1
% if nargin < 5
%     stagsel = 1;
% end
% 
% if nargin < 6
%     postfix = repmat({''}, 1, length(datasets));
% end

if isscalar(nsubtypes)
    nsubtypes = repmat(nsubtypes, [length(datasets), length(methods)]);
end

results = {};

Ps = cell(length(datasets), length(methods));
for idx_data = 1:length(datasets)
    dataset_name = datasets{idx_data};

    [data,mmse,time2event,dx_sel] = load_data_all(dataset_name);
    data.years = round(data.years, 4);

    for idx_method = 1:length(methods)   
        nsubtype = nsubtypes(idx_data, idx_method);

%         appendix = '_1x1';
        appendix = '';
        result_dir = ['output/',datasets{idx_data},'_',methods{idx_method},appendix];

        % if the postfix path (3_subtypes_run1) does not exist while the 
        % original path (3_subtypes) exist, use the original path
        postfix_exist = isfolder([result_dir, '/',num2str(nsubtype),'_subtypes',postfix]);
        original_exist = isfolder([result_dir, '/',num2str(nsubtype),'_subtypes']);
        if ~postfix_exist && ~original_exist
            results{idx_data, idx_method} = [];
            continue
        elseif postfix_exist
            postfix1 = postfix;
        elseif original_exist
            postfix1 = '';
        end

        elapsed_t = [];
        switch methods{idx_method}
            case {'FTR_kmeans', 'FTR_MCEM'}
                %% Change here to choose specific FTR version
%                 ftr_result_dir = [result_dir, '_min_max_0.97'];
                ftr_result_dir = result_dir;
                [subtype_stage, mdl, biomarker_names, P, elapsed_t] = load_ftr_result(...
                    ftr_result_dir, nsubtype, postfix1);
                Ps{idx_data, idx_method} = P;
            case {'sustain'}
                [subtype_stage, mdl, biomarker_names, elapsed_t] = load_sustain_result(...
                    result_dir, nsubtype, postfix1);
            otherwise
                [subtype_stage, mdl, biomarker_names] = load_other_result(...
                    result_dir, nsubtype, postfix1);
        end
        
        % round 'years' to 1e-4 precision before joining
        subtype_stage.years = round(subtype_stage.years, 4);
        
        joindata = outerjoin(subtype_stage,data,'Keys',{'RID','years'},'MergeKeys',true);

        joindata = sortrows(joindata, {'RID','years'});

        % select subjects
        case_indicator = joindata.group == 1;
        if dx_sel
            dx_ind = get_diagnosis_sel_indicators(joindata.RID, joindata.diagnosis);
            case_indicator = case_indicator & dx_ind;
        end
        control_indicator = joindata.group == 0;
        validation_indicator = joindata.group == 2;
        
        % % get permutation from HV/CTV
        % P = HV_CTV_perm(joindata{case_indicator, 'HV'}, ...
        %     joindata{case_indicator, 'CTV'}, ...
        %     joindata{case_indicator, 'subtype'});
        % 
        % Ps{idx_data, idx_method} = P;
        

        if include_control == 1
            joindata1 = joindata(case_indicator | control_indicator, :);
        elseif include_control == 0
            % select only the case subjects
            joindata1 = joindata(case_indicator, :);
        elseif include_control == 2
            joindata1 = joindata(validation_indicator, :);
        elseif include_control == 20
            joindata1 = joindata(validation_indicator | control_indicator, :);
        elseif include_control == 21
            joindata1 = joindata(validation_indicator | case_indicator, :);
        end
        
       
        result = [];
        result.elapsed_t = elapsed_t;
        result.joindata = joindata1;
        result.mdl = mdl;
        result.biomarker_names = biomarker_names;
        result.mmse = mmse;
        result.time2event = time2event;
        result.joindata_case = joindata(case_indicator, :);
        
        result = reorder_biomarkers(result);
        
        results{idx_data, idx_method} = result;
    end
end

% permute trajectories such that the trajectories from different datasets
% are comparable
for idx_method = 1:length(methods)
    for idx_data = 1:length(datasets)
        if ~strcmp(methods{idx_method}(1:3),'FTR')
            continue;
        end
        if isempty(results{idx_data, idx_method})
            continue
        end
        
        % If there is a permutation file specified, use it
        if ~isempty(Ps{idx_data, idx_method})
            P = Ps{idx_data, idx_method};
        % use the first dataset as reference and permute the remaining
        % datasets
        elseif idx_data > 1 && length(results{idx_data,idx_method}.biomarker_names)...
                == length(results{1,idx_method}.biomarker_names)
            % use the trajectories to determine the permutation
            traj1 = results{1,idx_method}.mdl.re_traj;
            traj_idx = results{idx_data,idx_method}.mdl.re_traj;
            
            if 0
                % use the mean pattern for permutation calculation
                [P,B] = permute_endmembers(squeeze(mean(traj1,2))', squeeze(mean(traj_idx,2))');
            else
                % use the flattened trajectory for permutation calculation
                reshape_size1 = [size(traj1,1)*size(traj1,2), size(traj1,3)];
                traj1_1 = reshape(traj1, reshape_size1)';
                reshape_size2 = [size(traj_idx,1)*size(traj_idx,2), size(traj_idx,3)];
                traj_idx_1 = reshape(traj_idx, reshape_size2)';
                
                [P,B] = permute_endmembers(traj1_1, traj_idx_1);
            end
        else % for the first dataset, reorder based on MMSE decrease rate
            joindata_case = results{idx_data,idx_method}.joindata_case;
            mmse1 = results{idx_data,idx_method}.mmse;
            % The order is determined by the MMSE rate of the joined big table instead of the
            % complete MMSE table. For OASIS, there is a discrepancy
            % between the order from the joined big table and the order
            % from the MMSE table. Hence, it is recommended to determine
            % the order of OASIS by ADNI instead of by its MMSE rate from
            % the joined big table
            % Update on Dec 18, 2024: Change it to order based on the
            % complete MMSE. It cannot be saved since if it is saved at
            % some run, running it again will use the previously save
            % permutation.
            nsubtype = nsubtypes(idx_data, idx_method);
            P = determine_trajectory_order(joindata_case, mmse1, nsubtype);
            result_dir = ['output/',datasets{idx_data},'_',methods{idx_method}];
            ind1 = P * (1:nsubtype)';
            %writematrix(ind1', [result_dir, '/',num2str(nsubtype),'_subtypes/permutation.csv']);
        end
        
        mdl = results{idx_data,idx_method}.mdl;


        joindata = results{idx_data,idx_method}.joindata;
        joindata_case = results{idx_data,idx_method}.joindata_case;

        [mdl1, joindata1] = permute_subtypes(mdl, joindata, P);
        [mdl1, joindata_case1] = permute_subtypes(mdl, joindata_case, P);
        
        results{idx_data,idx_method}.mdl = mdl1;
        
        results{idx_data,idx_method}.joindata = joindata1;
        results{idx_data,idx_method}.joindata_case = joindata_case1;
        results{idx_data,idx_method}.P = P;
        

    end
end

   
end



function result = reorder_biomarkers(result)
% reorder biomarkers in alphabetical order
biomarker_names = result.biomarker_names;
biomarker_names = strrep(biomarker_names, '-', '_');

[biomarker_names1, idx] = sort(biomarker_names);
% biomarker_names1 = cellfun(@lower, biomarker_names1, 'UniformOutput', false);
result.biomarker_names = biomarker_names1;

joindata = result.joindata;
joindata.Properties.VariableNames = strrep(joindata.Properties.VariableNames, '-', '_');

joindata1 = joindata;
joindata1(:,biomarker_names) = joindata(:,biomarker_names1);
joindata1 = renamevars(joindata1, biomarker_names, biomarker_names1);
result.joindata = joindata1;

if isfield(result.mdl, 're_traj') % FTR
    result.mdl.re_traj = result.mdl.re_traj(idx,:,:);
end

end

%% load functions for different methods
function [subtype_stage, mdl, biomarker_names, P, elapsed_t] = load_ftr_result(...
    result_ftr_path, nsubtype, postfix)


if exist([result_ftr_path, '/',num2str(nsubtype),'_subtypes',postfix,'/permutation.csv'],'file')
    order = readmatrix([result_ftr_path, '/',num2str(nsubtype), ...
        '_subtypes',postfix,'/permutation.csv']);
    P = zeros(max(order));
    for i = 1:max(order)
        P(i, order(i)) = 1;
    end
else
    P = [];
end
% load the model parameters
mdl = load_parameters_theta([result_ftr_path,'/',num2str(nsubtype),'_subtypes',postfix]);

S = load([result_ftr_path, '/all_data.mat']);
biomarker_names = S.biomarker_name;

subtype_stage = load_subtype_stage(result_ftr_path, nsubtype, postfix);

elapsed_t = readmatrix([result_ftr_path, '/',num2str(nsubtype),'_subtypes',postfix,'/elapsed_t.csv']);

end

function subtype_stage = load_subtype_stage(result_ftr_path, nsubtype, postfix)

result_subtype_stage_path = [result_ftr_path, '/',num2str(nsubtype), ...
    '_subtypes',postfix,'/subtype_stage.csv'];

% load the subtypes and stages
subtype_stage = readtable(result_subtype_stage_path,'VariableNamingRule','preserve');
subtype_stage = renamevars(subtype_stage, 'PTID', 'RID');
subtype_stage = renamevars(subtype_stage, 'years_from_baseline', 'years');
subtype_stage.Diagnosis = [];

% load the subtype probabilities
result_subtype_prob_path = [result_ftr_path, '/',num2str(nsubtype),'_subtypes', ...
    postfix, '/subtype_prob.csv'];
subtype_prob = readtable(result_subtype_prob_path,'VariableNamingRule','preserve');

assert(all(all(subtype_stage{:,1:2} == subtype_prob{:,1:2})));
subtype_names = cellfun(@(x) ['subtype', num2str(x)], num2cell(1:nsubtype), 'UniformOutput', 0); 
subtype_names_visit = cellfun(@(x) ['subtype', num2str(x),'_visit'], ...
    num2cell(1:nsubtype), 'UniformOutput', 0); 
subtype_stage = [subtype_stage, subtype_prob(:, [subtype_names,subtype_names_visit])];
end


function [subtype_stage, mdl, biomarker_names, elapsed_t] = load_sustain_result(...
    result_path, nsubtype, postfix)

% load the model parameters
S = load([result_path, '/',num2str(nsubtype),'_subtypes',postfix,'/mdl.mat']);
mdl = S.mdl;

S = load([result_path, '/all_data.mat']);
biomarker_names = S.biomarker_name;

subtype_stage = load_subtype_stage(result_path, nsubtype, postfix);

elapsed_t = readmatrix([result_path, '/',num2str(nsubtype),'_subtypes',postfix,'/elapsed_t.csv']);

end

function [subtype_stage, mdl, biomarker_names] = load_other_result(result_path, nsubtype, postfix)
result_subtype_stage_path = [result_path, '/',num2str(nsubtype),'_subtypes',...
    postfix, '/subtype_stage.csv'];

% load the model parameters
% S = load([result_path, '/',num2str(nsubtype),'_subtypes/mdl.mat']);
mdl = [];

S = load([result_path, '/all_data.mat']);
biomarker_names = S.biomarker_name;

% load the subtypes and stages
subtype_stage = readtable(result_subtype_stage_path,'VariableNamingRule','preserve');
subtype_stage = renamevars(subtype_stage, 'PTID', 'RID');
subtype_stage = renamevars(subtype_stage, 'years_from_baseline', 'years');
subtype_stage.Diagnosis = [];
end

function P = HV_CTV_perm(HV, CTV, subtype)

ratio = HV ./ CTV;
nsubtype = max(subtype);
ms = [];
for k = 1:nsubtype
    m = median(ratio(subtype == k));
    ms = [ms, m];
end
[~, I] = sort(ms);
P = zeros(k,k);
for i = 1:k
    P(i, I(i)) = 1;
end
end

