function [outputArg1,outputArg2] = check_stability(dataset_name, method_name, options)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
tic

if nargin < 1
    dataset_name = 'OASIS3_ROD1_HS';
end

if nargin < 3
    options = [];
end

if nargin < 2
    method_name = 'FTR_MCEM';
end

subtypes = [];

num_runs = 5;
for i = 1:num_runs
    [subtype_stage, mdl, options, train_data, test_data] = ...
        run_algo(dataset_name, method_name, options);
    subtype_i = subtype_stage.subtype(options.train_inds);
%     RID_i = train_data.RID;
%     [~,ia] = unique(RID_i);
%     subtypes(:,i) = subtype_i(ia);
    subtypes(:,i) = subtype_i;

    % save result folder
    save_dir = ['./output/',options.output_file_name];
    save_dir1 = [save_dir, '/',num2str(options.nsubtype),'_subtypes'];

    save_dir2 = [save_dir1, '_run', num2str(i)];

    copyfile(save_dir1, save_dir2, 'f');

    calc_ARI(subtypes, i);
end

toc



end

function calc_ARI(subtypes, num_runs)
for run1 = 1:num_runs
    for run2 = 1:num_runs
        sim_mat(run1, run2) = rand_index(subtypes(:,run1), subtypes(:,run2));
    end
end

disp('ARI between partitions');
sim_mat
end