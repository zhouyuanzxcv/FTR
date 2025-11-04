function bes_nsubtype = FTR_model_selection(dat,PTID,options)

num_folds = options.num_folds;
output_file_name = options.output_file_name;

nsubtype_list = options.nsubtype;

save_path = ['./output/',output_file_name];
if ~exist(save_path, 'dir')
    mkdir(save_path)
end

train_all_subtypes_folds(dat, PTID, nsubtype_list, num_folds, save_path, options);

nnsubtype = length(nsubtype_list);
loglik_mat = zeros(nnsubtype,num_folds);
BIC_mat = zeros(nnsubtype,1);
% rss_mat = zeros(nnsubtype,num_folds);

% Read result
for i = 1:nnsubtype
    nsubtype = nsubtype_list(i);
    for w = 1:num_folds
        result_save_dir = join([save_path,'/cross_validation_nsubtype',int2str(nsubtype),'_fold',int2str(w)]);
        result = readmatrix([result_save_dir,'/test_loglik.csv']);
        loglik_mat(i,w) =  result(1);
        % BIC_mat(i,w) =  result(2);
        % rss_mat(i,w) =  result(2);
    end
end

% BIC should be taken from the whole data (fold 0)
for i = 1:nnsubtype
    nsubtype = nsubtype_list(i);

    result_save_dir = join([save_path,'/cross_validation_nsubtype',int2str(nsubtype),'_fold0']);
    result = readmatrix([result_save_dir,'/train_loglik_BIC.csv']);

    BIC_mat(i) =  result(2);
end

% Save result    
loglik_cv = sum(loglik_mat,2);

[~,idx] = max(loglik_cv);
bes_nsubtype = nsubtype_list(idx);
header_cell = cell(1,nnsubtype);
for j = 1:nnsubtype
    header_cell(j) = cellstr(['nsubtype = ',int2str(nsubtype_list(j))]);
end

writetable(array2table(loglik_mat', 'VariableNames', header_cell), [save_path,'/loglikelihood_cross_validation','.csv']);
writetable(array2table(BIC_mat', 'VariableNames', header_cell), [save_path,'/BIC_cross_validation','.csv']);

% data point-wise model residues will be plotted and compared. No need
% to store the summary
% writetable(array2table(rss_mat', 'VariableNames', header_cell), [save_path,'/model_residue_cross_validation','.csv']);

fprintf('Best number of subtypes from model selection is %d.\n', bes_nsubtype);

end

function train_all_subtypes_folds(dat, PTID, nsubtype_list, num_folds, save_path, options)

% num_int = options.num_int;

[nsamp,nbiom] = size(dat);
[uni_id, ia, ic] = unique(PTID);
cv_inds_place = [save_path,'/cv_inds.csv'];

% solve the problem that a specific (nsubtype, fold) met an
% error in model selection
specify_nsubtype_fold_to_run = parse_param(options, 'specify_nsubtype_fold_to_run', []);

if isempty(specify_nsubtype_fold_to_run)
    if 1 % split based on subjects
        if 0 % use fixed splitting
            id_inds = repmat((1:num_folds)', 1, ceil(length(uni_id)/num_folds)); 
            id_inds = id_inds(:); 
            id_inds = id_inds(1:length(uni_id));
        else % use a random splitting strategy
            id_inds = crossvalind('Kfold',length(uni_id),num_folds);
        end
%         inds = zeros(nsamp,1);
%         for i = 1:nsamp
%             inds(i) = id_inds(uni_id==PTID(i));
%         end

        inds = id_inds(ic);
    else
        % in OASIS3, splitting based on subjects is unstable in cross
        % validation. So change it to split based on data points
        inds = crossvalind('Kfold', nsamp, num_folds);
    end
    writematrix(inds,cv_inds_place);
else 
    inds = readmatrix(join([save_path,'/cv_inds.csv']));
end


% select specified (nsubtype,fold) to run in case an error was
% met when running model selection
if isempty(specify_nsubtype_fold_to_run)
    % added a specific fold 0 to indicate using the whole data
    [m,n] = ndgrid(nsubtype_list, 0:num_folds);
    Z = [m(:),n(:)];
else
    Z = specify_nsubtype_fold_to_run;
end

if isempty(gcp('nocreate')), parpool; end

disp('Begin cross validation for model selection:');
% fprintf('%s\n', repmat('.',1,size(Z,1)));

start_t = tic;

for i = 1:size(Z,1)
    nsubtype = Z(i,1);
    fold = Z(i,2);
    save_dir = join([save_path,'/cross_validation_nsubtype',int2str(nsubtype),'_fold',int2str(fold)]);
    if ~exist(save_dir, 'dir')
        mkdir(save_dir)
    end
    try
        fprintf('Run cross validation with %d subtypes, fold %d\n', nsubtype, fold);
        % FTR algorithm
        if fold == 0
            test_inds = [];
            train_inds = (1:nsamp)';
        else
            test_inds = find(inds == fold);
            train_inds = find(inds ~= fold);
        end

        [mdl,train_subtype,train_stage,extra] = ...
            FTR_model(dat(train_inds,:),PTID(train_inds,:),nsubtype,options);
        % [traj,re_traj,train_subtype,train_stage,sigma,loglik_train,proption] = ...
            
        num_int = size(mdl.traj, 2);
        N_train = length(unique(PTID(train_inds,:)));
        D = (nsubtype-1) + (nsubtype*(num_int-1)*nbiom);
        if strcmp(options.sigma_type, '1xK')
            D = D + nsubtype;
        elseif strcmp(options.sigma_type, 'BxK')
            D = D + nsubtype*nbiom;
        end
        BIC = mdl.loglik - D/2 * log(N_train);

        % Save data
        save_parameters_theta(save_dir, mdl);
        
        writetable(array2table([train_inds,PTID(train_inds,:),train_subtype,train_stage], ...
            'VariableNames', {'Index' ,'PTID' ,'subtype','stage'}), [save_dir,'/train_subtype_stage.csv']);

        writetable(array2table([mdl.loglik, BIC], 'VariableNames', ...
            {'loglikelihood','BIC'}), strcat(save_dir,'/train_loglik_BIC.csv'));

        if ~isempty(test_inds)
            [test_stage,test_subtype,extra] = cal_stage_subtype(dat(test_inds,:), ...
                PTID(test_inds,:), mdl, options);
    
            loglik_test = cal_loglik(dat(test_inds,:), PTID(test_inds,:), mdl.traj, mdl.proption, mdl.sigma);

            % changed it to data point-wise mean absolute difference
%             rss = mean(sqrt(dist_samp_min),2);

            writetable(array2table([test_inds,PTID(test_inds,:),test_subtype,test_stage], ...
                'VariableNames', {'Index' ,'PTID' ,'subtype','stage'}), [save_dir,'/test_subtype_stage.csv']);

%             writetable(array2table([test_inds,PTID(test_inds,:),rss], ...
%                 'VariableNames', {'Index' ,'PTID' , 'model residue'}), [save_dir,'/test_model_residue.csv'])

            writetable(array2table([loglik_test], 'VariableNames', ...
                {'loglikelihood'}), strcat(save_dir,'/test_loglik.csv'));
        end

%         fprintf('.');
        if nsubtype == nsubtype_list(end)
            fprintf('Cross validation %.1f%% completed.\n', (fold+1)/(num_folds+1)*100);
            pause(0.01);
        end

    catch me
        fprintf('ERROR: Failed to run %d subtypes, fold %d in model selection.\n', nsubtype, fold);
        disp(me.message)
    end
end

fprintf('\nEnd cross validation (duration: %.2f hours)\n', toc(start_t)/3600);

end