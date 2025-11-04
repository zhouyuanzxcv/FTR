function [traj,re_traj,subtype,stage,sgm,loglik,proption,extra] = ...
    kmeans_subtype(dat,PTID,nsubtype,pre_subtype,pre_sgm,options)

max_iter = options.max_iter;
num_int = options.num_int;
% caloglik = false;

options.filter = parse_param(options, 'filter', 1);

Ts = calc_Ts(PTID);

nsubj = length(Ts);

biom_max = parse_param(options, 'biom_max', []);
if isempty(biom_max)
    biom_max = max(dat, [], 1);
end

if numel(biom_max) == 1
    biom_max = repmat(biom_max, 1, size(dat, 2));
end

[nsamp, nbiom] = size(dat);

traj = zeros(nbiom,num_int,nsubtype);
re_traj = zeros(nbiom,num_int,nsubtype);

%% initialize 
if length(pre_subtype) == nsamp
    subtype = pre_subtype;
else
%     subtype = randi(nsubtype,[nsamp 1]);
    subtype = init_subtype(nsubj, nsubtype, dat, PTID, options);
%     subtype = init_subtype(nsubj, nsubtype);
%     subtype = repelem(subtype, Ts, 1);
end

% if ismatrix(pre_sgm) && size(pre_sgm,1) == nbiom && size(pre_sgm,2) == nsubtype 
if ~isempty(pre_sgm)
    sgm = pre_sgm;
else
    sigma_init = parse_param(options, 'sigma_init', 1);
    switch options.sigma_type
        case {'1xK','1x1'}
            sgm = ones(1, nsubtype) * sigma_init;
        case 'BxK'
            sgm = ones(nbiom,nsubtype) * sigma_init;
    end
end

loglik_all = NaN(max_iter,1);
subtype_all = NaN(nsamp, max_iter);


% lambdas = parse_param(options, 'lambda_starts', 20*ones(1,nsubtype));


for iter = 1:max_iter

    %tic

    %     for j = 1:nbiom

    for k = 1:nsubtype
%         options.lambda_start = lambdas(k);

        dat_sub = dat(subtype == k,:);

        if sum(subtype == k) == 0 % cluster k is empty
            traj(:,:,k) = nan;
        elseif ~options.filter
%             [s,f,lambdas(k)] = ftr_base(dat_sub, num_int, 0);
            [s,f] = ftr_base(dat_sub, num_int, 0);
        else
            switch options.sigma_type
                case {'1xK','1x1'}
                    sgm1 = sgm(k);
                case 'BxK'
                    sgm1 = sgm(:,k);
            end

            [s,f] = ftr_base(dat_sub, num_int, sgm1, 0, biom_max, options);

        end
        traj(:,:,k) = f;
        
    end
    %     end

    %toc
    
    %% reparameterize f
    
    re_traj = traj;
    for k = 1:nsubtype
        traj_k = re_traj(:,:,k);
        re_traj(:,:,k) = reparam_traj(traj_k', size(traj,2), options)';
    end
    traj = re_traj;
    

    %% Calculate distance

    %tic

%     subtype_last = subtype;
    [stage,subtype,dist_samp_min] = cal_stage_subtype_kmeans(dat,PTID,traj,Ts);

    subtype_all(:,iter) = subtype;

    if options.filter || iter == max_iter
        % update sigma
        for k = 1:nsubtype
            dist_samp_min_k = dist_samp_min(subtype==k, :);
            
            switch options.sigma_type
                case '1xK'
                    sgm(:,k) = sqrt(mean(dist_samp_min_k(:)));
                case 'BxK'
                    sgm(:,k) = sqrt(mean(dist_samp_min_k,1));
            end
            
        end
        
        if strcmp(options.sigma_type, '1x1')
            sgm = sqrt(mean(dist_samp_min(:)));
            sgm = repmat(sgm, 1, nsubtype);
        end
    end
%     end
    
    if options.filter || iter == max_iter 
        proption = subtype_samp2proportion(subtype, PTID, nsubtype);
        proption = proption';
        if any(proption == 0) % one cluster becomes empty
            subtype(:) = NaN;
            break;
        end

        loglik = cal_loglik(dat,PTID,traj,proption,sgm);
        loglik_all(iter) = loglik;
    end

    %% stop if converged
%     if iter > 10 && all(all(subtype_all(:,iter-4:iter) == subtype_all(:, iter-5:iter-1)))
%         subtype_all(:,iter+1:end) = repmat(subtype_all(:,iter), 1, max_iter - iter);
%         loglik_all(iter+1:end) = repmat(loglik_all(iter), max_iter - iter, 1);
%         break;
%     end

    %toc

end

%% Reparametrization
for k = 1:nsubtype
    traj_k = traj(:,:,k);
%     traj_k = reparam_traj(traj_k')';
    re_traj(:,:,k) = traj_k;
end


extra = [];
extra.loglik_all = loglik_all;
extra.subtype_all = subtype_all;
% extra.lambdas = lambdas;
end