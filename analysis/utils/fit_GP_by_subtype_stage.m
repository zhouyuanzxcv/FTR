function re_traj = fit_GP_by_subtype_stage(data, stages, biomarker_names_pet)
K = length(unique(data.subtype));

nsubtype = K;

nbiom = length(biomarker_names_pet);
num_int = length(stages);

re_traj = zeros(nbiom,num_int,nsubtype);

meanfunc = [];
covfunc = @covSEiso; 
likfunc = @likGauss; 
hyp = struct('mean', [], 'cov', [0 0], 'lik', -1);

for k = 1:K
    lis = data.subtype == k;
    x = data.stage(lis);
    y = data{lis, biomarker_names_pet};
    xt = transpose(stages);
    hyp2 = minimize(hyp, @gp, -5, @infGaussLik, meanfunc, covfunc, likfunc,x, y);
    [yt ys] = gp(hyp2, @infGaussLik, meanfunc, covfunc, likfunc,x, y, xt);
    re_traj(:,:,k) = transpose(yt);
    %gprModel = fitrgp(data.stage(lis), data{lis,3:end});
    %re_traj(:,:,k) = predict(gprModel, transpose((0:num_int)/num_int));
    
end
end

