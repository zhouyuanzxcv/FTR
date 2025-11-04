function [traj_c,re_traj_c,subtype_c,stage_c,sigma_c,loglik_c,proption_c,extra_c] = ...
    init_wof(dat,PTID,nsubtype,options)

[nsamp, nbioms] = size(dat);

% options.filter = false;
max_ep = parse_param(options, 'max_ep', 100);
max_tries = parse_param(options, 'max_tries', 10);
% num_same_partition = parse_param(options, 'num_same_partition', 20);


% if set a single trajectory, no need to randomly initialize multiple times
if nsubtype == 1
    max_ep = 1;
    max_tries = 1;
end

max_runs = max_tries*max_ep;

subtype = cell(1, max_runs);
stage = cell(1, max_runs);

traj = cell(1, 1, 1, max_runs);
re_traj = cell(1, 1, 1, max_runs);
sgm = cell(1,1,max_runs);
loglik_mat = cell(1, max_runs);
proption = cell(1,1,max_runs);
extra = cell(1, max_runs);
loglik_all = cell(1, max_runs);
subtype_all = cell(1, 1, max_runs);

meet_stability_metric = 0;

% guarantees that the first time meta_labels_old is accessed in sel_init_by
% _meta_clustering, it is initialized to be all ones
if isfield(options, 'meta_labels_old')
    options = rmfield(options, 'meta_labels_old');
end


    
for try_idx = 1:max_tries
    s1 = (try_idx - 1) * max_ep;

    start_t = tic;

    if options.parfor
        % if 0
        parfor ep = 1:max_ep
            try
                [traj{:,:,:,s1+ep},re_traj{:,:,:,s1+ep},subtype{s1+ep},stage{s1+ep},...
                    sgm{:,:,s1+ep},loglik_mat{s1+ep},proption{:,:,s1+ep},extra{s1+ep},...
                    loglik_all{s1+ep},subtype_all{:,:,s1+ep}] = ...
                    choose_version(dat,PTID,nsubtype,[],[],options);
            catch me
                fprintf('WARNING: exception in init_wof at run %d/%d with %d subtypes.\n', ...
                    ep, max_ep, nsubtype);
                loglik_mat{s1+ep} = -Inf;
                subtype{s1+ep} = NaN(nsamp, 1);
                loglik_all{s1+ep} = NaN(options.max_iter, 1);

                msgText = getReport(me);
                disp(msgText);
            end

%             progress = ep+(try_idx-1)*max_ep;
%             if nsubtype > 1 && mod(progress, max_runs/10) == 0
%                 fprintf('Progressed to run %d out of total %d runs \n', progress, max_runs);
%             end
        end
    else
        for ep = 1:max_ep
%             try
            [traj{:,:,:,s1+ep},re_traj{:,:,:,s1+ep},subtype{s1+ep},stage{s1+ep},...
                sgm{:,:,s1+ep},loglik_mat{s1+ep},proption{:,:,s1+ep},extra{s1+ep},...
                loglik_all{s1+ep},subtype_all{:,:,s1+ep}] = ...
                choose_version(dat,PTID,nsubtype,[],[],options);
%             catch me
%                 fprintf('WARNING: exception in init_wof at run %d/%d with %d subtypes.\n', ...
%                     ep, max_ep, nsubtype);
%                 loglik_mat{s1+ep} = -Inf;
%                 subtype{s1+ep} = NaN(nsamp, 1);
%                 loglik_all{s1+ep} = NaN(options.max_iter, 1);
% 
%                 msgText = getReport(me);
%                 disp(msgText);
%             end

%             progress = ep+(try_idx-1)*max_ep;
%             if nsubtype > 1 && mod(progress, max_runs/10) == 0
%                 fprintf('Progressed to run %d out of total %d runs \n', progress, max_runs);
%             end
        end
    end

    % options.dat = dat;
    % options.PTID = PTID;

    s2 = try_idx * max_ep;

%     [~,ia] = unique(PTID);
    subtype_mat = cell2mat(subtype(1:s2));
    subtype_all_mat = cell2mat(subtype_all(:,:,1:s2));
%     subtype_mat = subtype_mat(ia,:);
%     subtype_all_mat = subtype_all_mat(ia,:,:);

    [best_idx, is_stable, best_metric, options] = sel_from_multi_inits(subtype_mat, ...
        cell2mat(sgm(:,:,1:s2)), cell2mat(loglik_all(1:s2)), ...
        subtype_all_mat, options);

    if nsubtype > 1
        tmp = options.similarity(triu(true(s2), 1));
        fprintf(['Run %d/%d takes %.2g mins, pairwise ARI IQR/max: %.2f-%.2f/%.2f, ', ...
            'meta-cluster prop.: %s, ', ...
            'meta mode: %d. \n'], ...
            s2, max_runs, toc(start_t)/60, prctile(tmp, 25), prctile(tmp,75), max(tmp), ...
            mat2str(options.meta_proption, 2), options.meta_metric);
    end

    % if nsubtype > 1 && best_metric >= num_same_partition && max_tries > 1
    if nsubtype > 1 && is_stable && max_tries > 1
        disp(['Meet stability metric ', num2str(best_metric), ...
            ' at ',num2str(s2), ' initializations']);
        meet_stability_metric = 1;

        break;
    end
end



if nsubtype > 1 && ~meet_stability_metric && max_tries > 1
    % disp(['Did not meet stability requirement (required ',  ...
    %     num2str(num_same_partition), ', last is ', num2str(best_metric), ...
    %     '). May need to rerun to check result stability or increase max_tries']);
    disp(['Did not meet stability requirement. ', ...
        'May need to rerun to check result stability or increase max_tries']);
end

traj_c = traj{:,:,:,best_idx};
re_traj_c = re_traj{:,:,:,best_idx};
subtype_c = subtype{best_idx};
stage_c = stage{best_idx};
sigma_c = sgm{:,:,best_idx};
loglik_c = loglik_mat{best_idx};
proption_c = proption{:,:,best_idx};
extra_c = extra{best_idx};

init_runs = [];
init_runs.subtype = cell2mat(subtype(1:s2));
init_runs.sigma = cell2mat(sgm(:,:,1:s2));
init_runs.loglik_all = cell2mat(loglik_all(1:s2));
init_runs.proption = cell2mat(proption(:,:,1:s2));
init_runs.best_idx = best_idx;
init_runs.best_metric = best_metric;

if isfield(options, 'inds_sel')
    init_runs.inds_sel = options.inds_sel;
end

extra_c.init_runs = init_runs;

if isfield(options, 'pre_subtype') && isfield(options, 'pre_sgm')
    extra_c.pre_subtype = options.pre_subtype;
    extra_c.pre_sgm = options.pre_sgm;
end

% pre_subtype = subtype(:,bes_ep);
% pre_sigma = sgm(:,:,bes_ep);

% extra = [];
% extra.loglik_all = loglik_all;


end

