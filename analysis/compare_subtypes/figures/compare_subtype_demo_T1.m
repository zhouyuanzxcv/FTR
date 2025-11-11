function [outputArg1,outputArg2] = compare_subtype_demo(data, method)
%COMPARE_SUBTYPE_DEMO Summary of this function goes here
%   Detailed explanation goes here
if nargin < 1
    data_names = {'ADNI_FSX_HS','OASIS3_ROD1_HS','NACC_HS'};
end

if nargin < 2
    method = {'FTR_MCEM'};
end

nsubtype = 3;
% data_sel = 3;

t_all = {};
for include_control = [0, 2]
    results = load_data_result(data_names, method, nsubtype, include_control, 1);
    
    t_all_data = table();
    for data_idx = [1,2,3]
        fprintf('----- include_control = %d, data_idx = %d\n', include_control, data_idx);
        data = results{data_idx,1}.joindata;
        mmse = results{data_idx,1}.mmse;

        % calculate slopes
        joindata = data(:,{'RID','subtype'});
        joindata = unique(joindata);
        joindata1 = outerjoin(joindata, mmse,'Type','Left','Keys','RID','MergeKeys',true);
        joindata = joindata1;
        [RID_unique, ia_mmse, ic_mmse] = unique(joindata.RID);
        mmse_bl = joindata(ia_mmse, :);
        betas = [];
        Ts_longitudinal_MMSE = zeros(size(mmse_bl.RID));
        for i = 1:length(mmse_bl.RID)
            mmse_i = joindata(joindata.RID == mmse_bl.RID(i), {'years', 'MMSE'});
            mmse_i = rmmissing(mmse_i);
            if size(mmse_i, 1) == 1
                betas(i) = NaN;
            else
                beta = regress(mmse_i.MMSE, [ones(size(mmse_i,1),1), mmse_i.years]);
                betas(i) = beta(2);
                Ts_longitudinal_MMSE(i) = size(mmse_i, 1);
            end
        end
        disp('Number of individuals having longitudinal MMSE')
        length(find(~isnan(betas)))
        MMSE_slope = betas';
        mmse_bl = addvars(mmse_bl, MMSE_slope);

        disp('Total time points with MMSE')
        sum(Ts_longitudinal_MMSE)



        num_patients_total = length(unique(data.RID));

        % calculate statistics of baseline age, sex, education, APOE4
        demo_all = [];

        baseline_ages = {};

        column_names = {};

        for k = 1:nsubtype
            demo = [];
            row_names = {};
            column_names{end+1} = sprintf('Subtype %d (%s)', k, data_names{data_idx});

            data_k = data(data.subtype == k, :);
            [unique_RID, ia, ic] = unique(data_k.RID);

            % count number of patients
            num_patients = size(ia,1);
            demo(end+1) = num_patients;
            row_names{end+1} = 'Number of patients';

            demo(end+1) = num_patients / num_patients_total;
            row_names{end+1} = '% of patients';

            % get age statistics
            baseline_age = data_k.AGE_baseline(ia);
            demo(end+1) = mean(baseline_age);
            row_names{end+1} = 'Mean baseline age';

            demo(end+1) = std(baseline_age);
            row_names{end+1} = 'Std baseline age';
            baseline_ages{k} = baseline_age;

            % get gender statistics
            num_male = length(find(data_k.PTGENDER(ia) == 1));
            demo(end+1) = num_male;
            row_names{end+1} = 'Male #';
            demo(end+1) = num_male / size(ia, 1);
            row_names{end+1} = 'Male Percentage';

            % get education statistics
            educ_k = data_k.PTEDUCAT(ia);
            demo(end+1) = mean(educ_k);
            row_names{end+1} = 'Education mean';

            demo(end+1) = std(educ_k);
            row_names{end+1} = 'Education std';

            % get APOE statistics
            apoe4 = data_k.APOE4(ia);
            apoe4_stats = [length(find(apoe4 == 0)), length(find(apoe4 == 1)), length(find(apoe4 == 2))];
            demo(end+1) = apoe4_stats(1);
            row_names{end+1} = '# of APOE4 == 0';
            demo(end+1) = apoe4_stats(2);
            row_names{end+1} = '# of APOE4 == 1';
            demo(end+1) = apoe4_stats(3);
            row_names{end+1} = '# of APOE4 == 2';

            demo(end+1:end+3) = apoe4_stats / sum(apoe4_stats);
            row_names{end+1} = '% of APOE4 == 0';
            row_names{end+1} = '% of APOE4 == 1';
            row_names{end+1} = '% of APOE4 == 2';

            % baseline diagnosis
            bl_dx = data_k.diagnosis_bl(ia);
            bl_dx = [length(find(bl_dx == 0)), length(find(bl_dx == 0.5)), length(find(bl_dx == 1))];
            demo(end+1) = bl_dx(1);
            demo(end+1) = bl_dx(2);
            demo(end+1) = bl_dx(3);
            row_names{end+1} = '# of diagnosis bl == 0';
            row_names{end+1} = '# of diagnosis bl == 0.5';
            row_names{end+1} = '# of diagnosis bl == 1';

            demo(end+1:end+3) = bl_dx / sum(bl_dx);
            row_names{end+1} = '% of diagnosis bl == 0';
            row_names{end+1} = '% of diagnosis bl == 0.5';
            row_names{end+1} = '% of diagnosis bl == 1';

            % diagnosis
            dx_counts = [length(find(data_k.diagnosis == 0)), ...
                length(find(data_k.diagnosis == 0.5)), ...
                length(find(data_k.diagnosis == 1))];
            demo(end+1:end+3) = dx_counts;
            row_names{end+1} = '# of dx == 0';
            row_names{end+1} = '# of dx == 0.5';
            row_names{end+1} = '# of dx == 1';

            demo(end+1:end+3) = dx_counts / sum(dx_counts);
            row_names{end+1} = '% of dx == 0';
            row_names{end+1} = '% of dx == 0.5';
            row_names{end+1} = '% of dx == 1';

            for cognitive_name = {'MMSE_bl', 'CDRSB_bl'}
                cognitive_name = cognitive_name{:};
                cognitive_k = mmse_bl(mmse_bl.subtype == k, :);
                mmse_bl_k = cognitive_k.(cognitive_name);
                mmse_bl_k = rmmissing(mmse_bl_k);
                demo(end+1) = mean(mmse_bl_k);
                demo(end+1) = std(mmse_bl_k);
                row_names{end+1} = [cognitive_name, ' (mean)'];
                row_names{end+1} = [cognitive_name, ' (std)'];
            end

            %     mmse_slope_k = cognitive_k.MMSE_slope;
            %     mmse_slope_k = rmmissing(mmse_slope_k);
            %     demo(end+1) = median(mmse_slope_k);
            %     demo(end+1) = prctile(mmse_slope_k, 25);
            %     demo(end+1) = prctile(mmse_slope_k, 75);
            %     row_names{end+1} = 'MMSE slope (median)';
            %     row_names{end+1} = 'MMSE slope (Q1)';
            %     row_names{end+1} = 'MMSE slope (Q3)';

            demo_all(:,k) = demo;
        end

        % disp('# of subjects, age mean, age std, male %, educ mean, educ std, # APOE4 = 0, 1, 2');
        demo_table = table(demo_all(:,1),demo_all(:,2),demo_all(:,3), ...
            'rownames', row_names, 'VariableNames', column_names);
        %disp(demo_table);

        t_all_data = cat(2, t_all_data, demo_table);

        [h, p12_bl_age] = ttest2(baseline_ages{1}, baseline_ages{2}, 'vartype', 'unequal');
        [h, p23_bl_age] = ttest2(baseline_ages{2}, baseline_ages{3}, 'vartype', 'unequal');
        [h, p13_bl_age] = ttest2(baseline_ages{1}, baseline_ages{3}, 'vartype', 'unequal');
        disp('p-values for differences between baseline ages (1 vs. 2, 2 vs. 3, 1 vs. 3)');
        [p12_bl_age, p23_bl_age, p13_bl_age]


        compare = [1,2;1,3;2,3];
        ps = [];
        tstats = [];
        for compare_ind = 1:size(compare, 1)
            k1 = compare(compare_ind, 1);
            k2 = compare(compare_ind, 2);
            cognitive_k1 = mmse_bl(mmse_bl.subtype == k1, :);
            cognitive_k2 = mmse_bl(mmse_bl.subtype == k2, :);
            mmse_bl_k1 = cognitive_k1.MMSE_bl;
            mmse_bl_k1 = rmmissing(mmse_bl_k1);
            mmse_bl_k2 = cognitive_k2.MMSE_bl;
            mmse_bl_k2 = rmmissing(mmse_bl_k2);

            [h,p,ci,stats] = ttest2(mmse_bl_k1, mmse_bl_k2, 'Vartype','unequal');
            ps(compare_ind) = p;
            tstats(compare_ind) = stats.tstat;
        end
    end

    t_all{end+1} = t_all_data;
end

disp('Discovery set');
disp(t_all{1});
disp('Validation set');
disp(t_all{2});

end

