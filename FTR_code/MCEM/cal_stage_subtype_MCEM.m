function [stage,subtype,extra] = cal_stage_subtype_MCEM(dat,PTID,re_traj,proption,sgm,Ts)

[loglik_re, sbtyp_nor_re, stage_nor_re, sbtyp_nor_it] = calc_log_likelihood_MCEM( ...
        dat, PTID, re_traj, proption, sgm);
    
[~, subtype] = max(sbtyp_nor_re, [], 2);
subtype = repelem(subtype, Ts, 1);

[num_samp, num_int, nsubtype] = size(stage_nor_re);

if 0
    stage_prob = get_stage_prob_from_subtype(stage_nor_re, subtype');
    
else
    sbtyp_nor_it1 = reshape(sbtyp_nor_it, [num_samp, 1, nsubtype]);
    stage_prob = sum(stage_nor_re .* repmat(sbtyp_nor_it1, [1, num_int, 1]), 3);
end

[~, stage] = max(stage_prob, [], 2);

% for output, we need to return the stage in [0,1] instead of the position
stage = (stage-1)/(num_int-1);

extra = [];
extra.subtype_prob = repelem(sbtyp_nor_re, Ts, 1);
extra.stage_prob = stage_prob;
extra.subtype_prob_per_visit = sbtyp_nor_it;

end


