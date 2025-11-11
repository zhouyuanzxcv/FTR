function [elapsed_ts1,sim_mats1] = compare_method_time(inputArg1,inputArg2)
%COMPARE_METHOD_TIME Summary of this function goes here
%   Detailed explanation goes here
close all


num_of_region = {'13','26','41','82'};
data = {'ADNI_FSX_LM', 'ADNI_FSX_LS','ADNI_FSX_HM', 'ADNI_FSX_HS'; ...
    'OASIS3_ROD1_LM', 'OASIS3_ROD1_LS','OASIS3_ROD1_HM', 'OASIS3_ROD1_HS'; ...
    'NACC_LM', 'NACC_LS', 'NACC_HM', 'NACC_HS'};

method = {'FTR_MCEM','sustain'};


[elapsed_ts, sim_mats] = get_time_stability(data(:), method);

elapsed_ts1 = reshape(elapsed_ts', [size(data,1)*length(method), size(data,2)]);
sim_mats1 = reshape(sim_mats', [size(data,1)*length(method), size(data,2)]);


ts_stats = create_stats_table(data, method, num_of_region, 'runtime', elapsed_ts1, 'mean-std');
disp(ts_stats);

stability_stats = create_stats_table(data, method, num_of_region, 'stability', sim_mats1, 'median-IQR');
disp(stability_stats);

end

function stats_table = create_stats_table(data_names, method_names, num_of_region, ...
    data_type, elapsed_ts1, mode)
stats_all = cell(1, size(data_names,2));

column_names = num_of_region;

row_names = {};
for i = 1:size(elapsed_ts1, 1)
    for j = 1:size(elapsed_ts1, 2)
        row_name = sprintf('%s_%s_%s', data_type, data_names{round(i/2),1}(1:end-3), ...
                    method_names{mod(i+1,2)+1});

        if strcmp(mode, 'mean-std')
            stats_all{j}(end+1,:) = mean(elapsed_ts1{i,j});        
            stats_all{j}(end+1,:) = std(elapsed_ts1{i,j});
            
            if j == 1
                row_name1 = [row_name, '_mean'];
                row_names{end+1} = row_name1;
                row_name1 = [row_name, '_std'];
                row_names{end+1} = row_name1;
            end
        elseif strcmp(mode, 'median-IQR')
            stats_all{j}(end+1,:) = median(elapsed_ts1{i,j});
            stats_all{j}(end+1,:) = prctile(elapsed_ts1{i,j}, 25);
            stats_all{j}(end+1,:) = prctile(elapsed_ts1{i,j}, 75);

            if j == 1
                row_name1 = [row_name, '_median'];
                row_names{end+1} = row_name1;
                row_name1 = [row_name, '_Q1'];
                row_names{end+1} = row_name1;
                row_name1 = [row_name, '_Q3'];
                row_names{end+1} = row_name1;
            end
        end
    end
end

stats_table = table(stats_all{:}, 'rownames', row_names, 'VariableNames', column_names);
end