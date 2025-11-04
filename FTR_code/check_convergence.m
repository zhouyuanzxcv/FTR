function [outputArg1,outputArg2] = check_convergence(loglik_all,subtype_all)
%CHECK_CONVERGENCE Summary of this function goes here
%   Detailed explanation goes here

mismatch = double(subtype_all(:,2:end,:) ~= subtype_all(:,1:end-1,:));
mismatch_perc = sum(mismatch, 1) / size(mismatch, 1);
mismatch_perc = squeeze(mismatch_perc);

figure;
subplot(1,2,1);
plot(loglik_all);

subplot(1,2,2);
plot(mismatch_perc);
end

