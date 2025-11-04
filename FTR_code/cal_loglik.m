function [loglik, log_pdf, log_int, log_int_grouped] = cal_loglik(dat, PTID, traj,proption,sgm)
[nsamp,nbiom] = size(dat);
[~,num_int,nsubtype] = size(traj);

if size(sgm,1) == 1 
    sgm = repmat(sgm, nbiom, 1);
end


step_mat = ones(num_int,1)/num_int;
log_pdf = zeros(nsamp,num_int,nsubtype);
log_int = zeros(nsamp,nsubtype);

for k = 1:nsubtype
    Sigma_k = diag(sgm(:,k).^2);
    for w = 1:num_int
        log_pdf(:,w,k) = logmvn(dat, traj(:,w,k)', Sigma_k);
    end
    log_int(:,k) = calc_log_mixture(log_pdf(:,:,k),step_mat');
end

if 0
    log_int_grouped = cell(1,nsubtype);
    unique_PTIDs = unique(PTID);
    for k = 1:nsubtype
        tmp = accumarray(PTID, log_int(:,k), [], [], [], 1);
        tmp1 = full(tmp(unique_PTIDs));
    %     tmp1 = nonzeros(tmp);
        log_int_grouped{k} = tmp1;
    end
    log_int_grouped = cell2mat(log_int_grouped);
else
    % This implementation is simpler and yields the same result as above
    [unique_PTIDs,ia,ic] = unique(PTID);
    [yy, xx] = ndgrid(ic, 1:nsubtype);
    log_int_grouped = accumarray([yy(:),xx(:)], log_int(:));
end


%[~,es_subtyp] = max(log_int,[],2);
loglik = sum(calc_log_mixture(log_int_grouped,proption'),1);

end
