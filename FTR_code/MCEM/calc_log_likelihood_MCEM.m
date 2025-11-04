function [loglik, sbtyp_nor, stage_nor, sbtyp_nor_it] = ...
    calc_log_likelihood_MCEM(dat, PTID, traj, proption, sgm)

[loglik, log_pdf, log_int, log_int_grouped] = cal_loglik(dat, PTID, traj,proption,sgm);

nsubj = size(log_int_grouped, 1);

loglik_sbtyp = exp(log_int_grouped-max(log_int_grouped,[],2)).*repmat(proption',nsubj,1);
sbtyp_nor = loglik_sbtyp./sum(loglik_sbtyp,2);

n_data_pts = size(log_int, 1);
loglik_sbtyp_it = exp(log_int-max(log_int,[],2)).*repmat(proption',n_data_pts,1);
sbtyp_nor_it = loglik_sbtyp_it./sum(loglik_sbtyp_it,2);

loglik_stage = exp(log_pdf-max(log_pdf,[],2));
stage_nor = loglik_stage./sum(loglik_stage,2);
end
