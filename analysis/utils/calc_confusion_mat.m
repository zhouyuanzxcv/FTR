function [conf_matrix, consistent_subject_IDs, inconsistent_subject_IDs] ...
    = calc_confusion_mat(train_PTID, train_subtype, ...
    test_PTID, test_subtype, num_subtypes)
%CALC_CONFUSION_MAT Summary of this function goes here
%   Detailed explanation goes here
conf_matrix = zeros(num_subtypes);

conf_PTIDs = cell(num_subtypes);

% Fill in the confusion matrix based on intersection of PTIDs between train and test subtypes
for i = 1:num_subtypes
    for j = 1:num_subtypes
        % Find PTIDs of train_data with subtype i and test_data with subtype j
        train_PTID_subtype_i = unique(train_PTID(train_subtype ==  i));
        test_PTID_subtype_j = unique(test_PTID(test_subtype ==  j));

        % Count the number of common PTIDs between train and test subtypes
        common_PTIDs = intersect(train_PTID_subtype_i, test_PTID_subtype_j);
        conf_matrix(i, j) = numel(common_PTIDs);
        conf_PTIDs{i, j} = common_PTIDs;
    end
end

consistent_subject_IDs = [];
inconsistent_subject_IDs = [];

for i = 1:num_subtypes
    for j = 1:num_subtypes
        tmp = conf_PTIDs{i, j};
        if ~isempty(tmp)
            if i == j
                consistent_subject_IDs = cat(1, consistent_subject_IDs, tmp);
            else
                inconsistent_subject_IDs = cat(1, inconsistent_subject_IDs, tmp);
            end
        end
    end
end

end

