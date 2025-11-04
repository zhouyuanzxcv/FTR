function disp_cells = compare_mean_pattern(data, nsubtype, ...
    biomarker_names_pet, adjust_covs)


if nargin < 4
    adjust_covs = {'stage'};
end

data1 = data(:,[biomarker_names_pet, {'subtype'}, adjust_covs]);
data2 = rmmissing(data1);



disp('----- compare mean pattern (2-1,3-2,3-1) -----')
combs = [1,2;2,3;1,3];
for k = 1:size(combs,1)
    % fprintf('----- subtype %d ----- \n', k);
    
    k1 = combs(k,1);
    k2 = combs(k,2);
    data_k = data2(data2.subtype == k1 | data2.subtype == k2, :);
    data_k.subtype = double(data_k.subtype == k2);
    
    % make the reference group be the first element in the order
    data_k.subtype = categorical(data_k.subtype, [0,1]);

    % make 0.5 (MCI) be the first (reference) element in the order
%     data_k.diagnosis = categorical(data_k.diagnosis, circshift(unique(data_k.diagnosis), 2));

    
    ps = [];
    ts = [];
    for b = 1:length(biomarker_names_pet)
        adjust_covs_str = strjoin(adjust_covs, ' + ');
        formula = [biomarker_names_pet{b},' ~ subtype + ',adjust_covs_str];
        lm = fitlm(data_k, formula);

        ps(b) = lm.Coefficients{'subtype_1','pValue'};
        ts(b) = lm.Coefficients{'subtype_1','tStat'};
    end
%     ind_reject = multi_test_correction(ps, 0.05, 'BH');
% %     ind_reject = ps < 0.05;
%     significant_biomarkers_k = biomarker_names_pet(ind_reject)';
%     significant_ts_k = ts(ind_reject);
%     significant_biomarkers_k(significant_ts_k > 0)
%     significant_biomarkers_k(significant_ts_k < 0)

    ps_all{k} = ps;
    ts_all{k} = ts;
end

use_fdr_threshold = 1;
disp_cells = show_subtype_diff(biomarker_names_pet, ps_all, ts_all, nsubtype, use_fdr_threshold);

end

% function compare_mean_pattern_by_regressing_out_dx_stage(data, nsubtype, biomarker_names_pet)
% biom_all_ro = {};
% regress_out = 0;
% 
% for k = 1:nsubtype
%     data_k = data(data.subtype == k, :);
%     biom_k = [];
%     for b = 1:length(biomarker_names_pet)
%         biom = data_k{:,biomarker_names_pet(b)};
%         stage = data_k.stage;
%         dx = data_k.diagnosis;
%         if regress_out
%             biom1 = regress_out_stage_diagnosis(biom, stage, dx);
%         else
%             biom1 = biom;
%         end
%         biom_k(:,b) = biom1;
%     end
%     
%     biom_avg_k = mean(biom_k, 1);
%     
%     biom_all_ro{k} = biom_k;
% end
% 
% one_vs_others = 0;
% if one_vs_others
%     for k = 1:nsubtype
%         fprintf('----- subtype %d ----- \n', k);
%         subtype_k = biom_all_ro{k};
%         subtype_others = biom_all_ro(1:nsubtype ~= k);
%         subtype_others = cell2mat(subtype_others');
%         
%         ps = [];
%         ts = [];
%         for b = 1:size(subtype_k,2)
%             y1 = subtype_k(:,b);
%             y2 = subtype_others(:,b);
%             [h,p,ci,stats] = ttest2(y1,y2,'Vartype','unequal');
%             tscore = stats.tstat;
%             ps(b) = p;
%             ts(b) = tscore;
%         end
%         ind_reject = multi_test_correction(ps, 0.05, 'BH');
%     %     biomarker_names_pet(ind_reject)'
%     %     ts(ind_reject)
%     
%         significant_biomarkers_k = biomarker_names_pet(ind_reject)';
%         significant_ts_k = ts(ind_reject);
%         significant_biomarkers_k(significant_ts_k > 0)
%         significant_biomarkers_k(significant_ts_k < 0)
%     end
% else
%     disp('---- Compare mean patterns (2-1, 3-2, 3-1) ----');
%     ps_all = {};
%     ts_all = {};
%     comp_ks = [1,2;2,3;1,3];
%     for i = 1:size(comp_ks,1)
%         k1 = comp_ks(i,1);
%         k2 = comp_ks(i,2);
% %         k_1 = mod(k,nsubtype)+1;
% %         fprintf('----- subtype %d - subtype %d -----\n', k_1, k);
%         subtype_k = biom_all_ro{k1};
%         subtype_others = biom_all_ro(k2);
%         subtype_others = cell2mat(subtype_others');
%         
%         ps = [];
%         ts = [];
%         for b = 1:size(subtype_k,2)
%             y1 = subtype_k(:,b);
%             y2 = subtype_others(:,b);
%             [h,p,ci,stats] = ttest2(y2,y1,'Vartype','unequal');
%             tscore = stats.tstat;
%             ps(b) = p;
%             ts(b) = tscore;
%         end
% 
%         ps_all{i} = ps;
%         ts_all{i} = ts;
% %         ind_reject = multi_test_correction(ps, 0.05, 'BH');
%     %     biomarker_names_pet(ind_reject)'
%     %     ts(ind_reject)
%     
% %         significant_biomarkers_k = biomarker_names_pet(ind_reject)';
% %         significant_ts_k = ts(ind_reject);
% %         significant_biomarkers_k(significant_ts_k > 0)
% %         significant_biomarkers_k(significant_ts_k < 0)
%     end
%     show_subtype_diff(biomarker_names_pet, ps_all, ts_all, nsubtype);
% end
% 
% end
