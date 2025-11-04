function P = determine_trajectory_order(joindata, mmse, nsubtype)

% use MMSE progression rate to determine the order
dependentVariable = {'MMSE'};

covs = {'AGE_baseline', 'PTGENDER', 'PTEDUCAT', 'diagnosis_bl'};
independentVariables = [{'RID','years', 'subtype'}, covs];
joindata = joindata(:,[{'RID','subtype'},covs]);
joindata = unique(joindata);
joindata1 = outerjoin(joindata, mmse,'Type','Left','Keys','RID','MergeKeys',true);
joindata = joindata1;

% dx_bl1 = get_baseline_diagnosis(joindata);
dx_bl = joindata.diagnosis_bl;

data_sel_inds = dx_bl == 0.5 | dx_bl == 1;

data = joindata(data_sel_inds, [dependentVariable, independentVariables]);

data = rmmissing(data);

try
    [slopes_all, beta_all, intercept_all, beta_p_all] = lmem(data, dependentVariable, nsubtype);
catch me
    disp(['WARNING! LMEM cannot determine the order in determine_trajectory_order. ',...
        'Use AGE_baseline to determine the order']);
    msgText = getReport(me);
    disp(msgText);

    beta_all = determine_order_by_age_bl(joindata, nsubtype);
end


% the decrease of MMSE is sorted such that the largest decrease lies at the
% last position
[beta_all1, inds] = sort(beta_all, 'descend');

P = eye(nsubtype);
P = P(inds, :);

end

function beta_all = determine_order_by_age_bl(joindata, nsubtype)
[RIDs, ia] = unique(joindata.RID);
joindata = joindata(ia, :);
age_means = [];
for k = 1:nsubtype
    age_means(k) = mean(joindata.AGE_baseline(joindata.subtype == k));
end
beta_all = age_means;
end

