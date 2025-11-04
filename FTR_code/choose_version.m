function [traj,re_traj,subtype,stage,sigma,loglik,proption,extra,loglik_all,subtype_all] = ...
    choose_version(dat,PTID,nsubtype,pre_subtype,pre_sigma,options)

if strcmp(options.methods,'kmeans')
    [traj,re_traj,subtype,stage,sigma,loglik,proption,extra]...
        = kmeans_subtype(dat,PTID,nsubtype,pre_subtype,pre_sigma,options);
elseif strcmp(options.methods,'MCEM')
    [traj,re_traj,subtype,stage,sigma,loglik,proption,extra]...
        = MCEM_subtype(dat,PTID,nsubtype,pre_subtype,pre_sigma,options);
end

loglik_all = extra.loglik_all;
subtype_all = extra.subtype_all;

end