function disp_cells = show_subtype_diff(biomarker_names_pet, ps_all, ...
    ts_all, nsubtype, use_fdr_threshold)
disp_cells = cell(length(biomarker_names_pet), nsubtype + 1);
disp_cells(:,1) = biomarker_names_pet;

use_fdr_correction = use_fdr_threshold;

for k = 1:length(ps_all)
    ps = ps_all{k};
    ts = ts_all{k};

    if use_fdr_correction
        ind_reject = find(multi_test_correction(ps, 0.05, 'BH'));
    else
        ind_reject = (1:length(ps))';
    end

    significant_ts_k = ts(ind_reject);
%     disp_cells{ind_reject, k+1} = significant_ts_k; % error for cell array

    for i = 1:length(ind_reject)
        disp_cells{ind_reject(i), k+1} = significant_ts_k(i);
    end
end

end


