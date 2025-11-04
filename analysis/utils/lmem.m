function [slopes_all, beta_all, intercept_all, beta_p_all, beta_t_all] = lmem(...
    data, dependentVariable, nsubtype, reference_group)

if nargin < 4
    reference_group = 1;
end


subtype = data.subtype;

% make the reference group be the first element in the order
data.subtype = categorical(data.subtype, circshift((1:nsubtype), -reference_group+1));

% make 0.5 (MCI) be the first (reference) element in the order
data.diagnosis_bl = categorical(data.diagnosis_bl, circshift(unique(data.diagnosis_bl), 2));

% data.diagnosis_bl = string(data.diagnosis_bl);

% Combine the predictors

% demeaning years has effects on the results because years are used also in
% the random effects. Suppose there is only 1 subtype. Not demeaning years
% may make the model put more efforts on minimizing the random intercepts.
% Demeaning years may make the model put more efforts on minimizing the
% random slopes.
mean_years = mean(data.years);
data.years = data.years - mean_years;
% mean_years = 0;

% whether demeaning baseline age, gender, education has not effects on the
% results
data.AGE_baseline = data.AGE_baseline - mean(data.AGE_baseline);
data.PTGENDER = data.PTGENDER - mean(data.PTGENDER);
data.PTEDUCAT = data.PTEDUCAT - mean(data.PTEDUCAT);


% Specify the formula for the linear mixed-effects model
% formula = strcat(dependentVariable{1},' ~ years*subtype + AGE + ', ...
%     'PTGENDER + PTEDUCAT + diagnosis + stage + (1 | RID) + (years-1 | RID)');

if length(unique(data.diagnosis_bl)) == 1 % diagnosis_bl takes only 1 value
    formula = strcat(dependentVariable{1},' ~ years*subtype + AGE_baseline + ', ...
        'PTGENDER + PTEDUCAT + (1 | RID) + (years-1 | RID)');
else
    formula = strcat(dependentVariable{1},' ~ years*subtype + AGE_baseline + diagnosis_bl + ', ...
        'PTGENDER + PTEDUCAT + (1 | RID) + (years-1 | RID)');
end

% Fit the linear mixed-effects model
mdl = fitlme(data, formula,'DummyVarCoding','reference');
[beta, beta_names, stats] = fixedEffects(mdl);
[re, re_names] = randomEffects(mdl);

beta_ref = beta(strcmp(beta_names{:,1},'years'));
intercept_ref = beta(strcmp(beta_names{:,1}, '(Intercept)'));

slopes_all = {};
beta_all = [];
intercept_all = [];
beta_p_all = [];
beta_t_all = [];

for k = 1:nsubtype
    RID_k = data.RID(subtype == k);
    
    % get groupwise slope (pvalue) and intercept
    if k == reference_group
        beta_k = beta_ref;
        intercept_k = intercept_ref;
        beta_p = NaN;
        beta_t = NaN;
    else
        beta_ind = strcmp(beta_names{:,1},['years:subtype_',num2str(k)]);
        beta_k = beta(beta_ind);
        beta_p = mdl.Coefficients.pValue(beta_ind);
        beta_t = mdl.Coefficients.tStat(beta_ind);

        intercept_k = beta(strcmp(beta_names{:,1}, ['subtype_',num2str(k)]));


        beta_k = beta_k + beta_ref;
        intercept_k = intercept_k + intercept_ref;
    end
    
    beta_all(1,k) = beta_k;
    intercept_all(1,k) = intercept_k;
    beta_p_all(1,k) = beta_p;
    beta_t_all(1,k) = beta_t;

    % get individual slopes
    
%     [unique_RID, ia, ic] = unique(RID_k);
%     slopes = zeros(size(unique_RID));
    
    re_names_RID = cellfun(@(x) str2num(x), re_names{:,'Level'});
    valid_RID = ismember(re_names_RID, RID_k);
    valid_years = strcmp(re_names{:,'Name'}, 'years');
    slopes = re(valid_RID & valid_years);
    slopes = slopes + beta_k;
    
%     for i = 1:length(unique_RID)
%         each_RID = unique_RID(i);
%         
%         ind1 = each_RID == cellfun(@(x) str2num(x), re_names{:,'Level'});
%         ind2 = strcmp(re_names{:,'Name'}, 'years');
%         slopes(i) = re(ind1 & ind2);
%         slopes(i) = slopes(i) + beta_k;
%     end
    slopes_all{1,k} = slopes;
end

% since the intercept is calculated using the years demeaned data, the
% intercept should be adjusted to be the intercept at 0
for k = 1:nsubtype
    intercept_all(k) = intercept_all(k) - beta_all(k) * mean_years;
end

end