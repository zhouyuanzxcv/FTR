function input1 = normalize_to_zscore(input, idx_CN, biom_column, mri_column)
input1 = input;
nbiom = length(biom_column);

% normalize to z-score for all biomarkers
sigma = zeros(nbiom,1);
for i = 1:length(biom_column)
    bm = biom_column{i};
    pd = fitdist(input{idx_CN, bm},'Normal');
    input1{:, bm} = (input{:, bm} - pd.mu)/pd.sigma;
    sigma(i) = pd.sigma;
end

% add a minus sign for MRI
for bm = mri_column
    input1{:, bm} = -input1{:, bm};
end

end
