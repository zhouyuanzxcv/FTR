function [elapsed_ts, sim_mats] = get_time_stability(dataset_name, method_name)
options = [];

nsubtype = 3;
include_control = 0;

subtypes = {};
elapsed_ts = {};

num_runs = 5;
for i = 1:num_runs
    options.postfix = ['_run',num2str(i)];
    results = load_data_result(dataset_name, method_name, nsubtype, include_control, options);
    
    for data_sel = 1:length(dataset_name)
        for method_sel = 1:length(method_name)
            if ~isempty(results{data_sel, method_sel})
                joindata = results{data_sel, method_sel}.joindata;
                subtypes{data_sel,method_sel}(:,i) = joindata.subtype;
                elapsed_ts{data_sel,method_sel}(:,i) = results{data_sel, method_sel}.elapsed_t;
            end
        end
    end    
end

sim_mats = {};
for data_sel = 1:length(dataset_name)
    for method_sel = 1:length(method_name)
        if ~isempty(subtypes{data_sel, method_sel})
            sim_mat = calc_ARI(subtypes{data_sel, method_sel}, num_runs);
            sim2 = sim_mat(triu(true(num_runs),1));
            sim_mats{data_sel,method_sel} = sim2;
        end
    end
end
end

function sim_mat = calc_ARI(subtypes, num_runs)
for run1 = 1:num_runs
    for run2 = 1:num_runs
        sim_mat(run1, run2) = rand_index(subtypes(:,run1), subtypes(:,run2));
    end
end

% disp('ARI between partitions');
% sim_mat
end

