function [mdl,subtype,stage,extra] = FTR_model(dat,PTID,nsubtype,options)

if 0 % search for lambda
    params = [];
    params.dat = dat;
    params.PTID = PTID;
    params.nsubtype = nsubtype;
    params.options = options;


    lambdas = [0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000];
    start_lambda = parse_param(options, 'lambda_start', 10);
    fprintf('Begin searching lambda. Start with lambda = %g. \n', start_lambda);
        
    eval_fun = @eval_obj_fun;
    [lambda_best, extra_best, obj_vals, best_ind] = search_for_lambda_with_step(...
        params, lambdas, start_lambda, eval_fun);

    fprintf('Lambdas: %s \n', mat2str(lambdas));
    fprintf('Obj vals: %s \n', mat2str(obj_vals,5));
    fprintf('Finish search. Best lambda: %g, best objective value: %.5g \n', ...
        lambda_best, obj_vals(best_ind));

    mdl = extra_best.mdl;
    subtype = extra_best.subtype;
    stage = extra_best.stage;
    extra = extra_best.extra;
else % use predefined lambda
    [mdl,subtype,stage,extra] = FTR_model_impl(dat,PTID,nsubtype,options);
end

end

function [mdl,subtype,stage,extra] = FTR_model_impl(dat,PTID,nsubtype,options)

verify_PTID_order = 1;
for i = 2:length(PTID)
    if ~(PTID(i) == PTID(i-1) || PTID(i) > PTID(i-1))
        verify_PTID_order = 0;
    end
end
assert(verify_PTID_order, ['PTID is not monotonically increasing in FTR_model. ',...
    'Possible error may occur in cal_loglik.m and MCEM_subtype.m']);

options.num_int = 101;

[traj,re_traj,subtype,stage,sigma,loglik,proption,extra] ...
    = init_wof(dat,PTID,nsubtype,options);

% if nsubtype > 1
%     % for kmeans, pre_subtype is a column vector consisting of averaged
%     % subtypes. For MCEM, pre_subtype is a matrix consisting of permuted
%     % subtypes from the largest meta-cluster. For both cases, pre_sigma is
%     % the averaged sigma
%     pre_subtype = extra.pre_subtype;
%     pre_sigma = extra.pre_sgm;
%     
%     if strcmp(options.methods,'MCEM')   
%         % if pre_subtype has less partitions than samp_multi, select the
%         % partitions in a circular way
%         sel_inds = mod(0:options.samp_multi-1, size(pre_subtype,2)) + 1;
%         pre_subtype = pre_subtype(:, sel_inds);
%     end
% 
%     options.max_iter = parse_param(options, 'max_iter_refine', 10);
% 
%     init_runs = extra.init_runs;
%     % methods
%     [traj,re_traj,subtype,stage,sigma,loglik,proption,extra,loglik_all,subtype_all] = ...
%         choose_version(dat,PTID,nsubtype,pre_subtype,pre_sigma,options);
%     extra.init_runs = init_runs;
% end

mdl = [];
mdl.traj = traj;
mdl.re_traj = re_traj;
mdl.sigma = sigma;
mdl.proption = proption;
mdl.loglik = loglik;

% extra.loglik_all_init = extra_init.loglik_all;

end    



function [obj_val, extra1] = eval_obj_fun(params, lambda)
% Need to change it to cross validation. But the current version is already
% too slow.
% solve the problem 
dat = params.dat;
PTID = params.PTID;
nsubtype = params.nsubtype;
options = params.options;

options.lambda = lambda;
[mdl,subtype,stage,extra] = FTR_model_impl(dat,PTID,nsubtype,options);


extra1 = [];
extra1.mdl = mdl;
extra1.subtype = subtype;
extra1.stage = stage;
extra1.extra = extra;

% evaluate the solution

[loglik, log_pdf, log_int, log_int_grouped] = cal_loglik(dat, PTID, ...
    mdl.re_traj, mdl.proption, mdl.sigma);
obj_val = -loglik;

end