function biom1 = regress_out_stage_diagnosis(biom, stage, dx)
one = ones(length(stage),1);

% change dx to one-hot encoding with MCI as reference. If dx is nan, it is
% treated as MCI, and no correction with respect to the dx is made.
dx1 = double(dx == [0.5,0,1]);
dx1(isnan(dx), :) = NaN;
dx1(:,1) = [];

covs = [stage,dx1];
mean_covs = [mean(stage), 0, 0];

% regress biomarker values onto each stage and diagnosis
[b,~,~,~,stats] = regress(biom,[stage,one]);
p_cov(1) = stats(3);

[b,~,~,~,stats] = regress(biom, [dx1,one]);
p_cov(2:3) = stats(3);

% regress out covariates
ro_tf = p_cov < 0.05 * ones(size(p_cov));

[b,~,~,~,stats] = regress(biom,[covs(:,ro_tf),one]);

% if any(ro_tf == 0)
%     disp('');
% end

% if none of the covariates is significant
if sum(ro_tf) == 0
    biom1 = biom;
else
    biom1 = biom - (covs(:,ro_tf) - mean_covs(:,ro_tf)) * b(1:end-1);
end


end
