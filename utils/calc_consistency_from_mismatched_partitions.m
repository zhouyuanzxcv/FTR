function consistency = calc_consistency_from_mismatched_partitions(clustering1, clustering2)

% length of clustering 1 should be less than the length of clustering 2
if length(clustering1) > length(clustering2)
    tmp = clustering2;
    clustering2 = clustering1;
    clustering1 = tmp;
end

% Get unique labels from each clustering
unique_labels1 = unique(clustering1);
unique_labels2 = unique(clustering2);

% Ensure both clusterings have the same number of clusters
if length(unique_labels1) ~= length(unique_labels2)
    error('Clusterings must have the same number of unique labels');
end

% Generate all possible permutations of the second clustering's labels
num_clusters = length(unique_labels1);
perms = generate_permutations(unique_labels2);

% Initialize variables to track best permutation and alignment
best_ari = -inf;
best_perm = [];
best_indices1 = [];
best_indices2 = [];

lcs_length_best = -inf;
permuted_clustering2_best = [];

for i = 1:size(perms, 1)
    % Create a mapping from original labels to permuted labels
    label_map = containers.Map(unique_labels2, perms(i, :));
    
    % Apply the permutation to the second clustering
    permuted_clustering2 = zeros(size(clustering2));
    for j = 1:length(clustering2)
        permuted_clustering2(j) = label_map(clustering2(j));
    end
    
    % Find the longest common subsequence (LCS) between the two clusterings
    [lcs_length, lcs_indices1, lcs_indices2] = findLCS(clustering1, permuted_clustering2);
    
    % Only proceed if there's a non-empty LCS
    if lcs_length > lcs_length_best
        % Compute ARI for this permutation and alignment
%         current_ari = rand_index(clustering1(lcs_indices1), permuted_clustering2(lcs_indices2));
        
        % Update best permutation and alignment
%         if current_ari > best_ari
        lcs_length_best = lcs_length;
%             best_ari = current_ari;
        best_perm = perms(i, :);
        best_indices1 = lcs_indices1;
        best_indices2 = lcs_indices2;

        permuted_clustering2_best = permuted_clustering2;
%         end
    end
end


consistency = lcs_length_best / length(clustering1);

% match = zeros(length(clustering1), length(clustering2));
% for i = 1:length(best_indices1)
%     match(best_indices1(i), best_indices2(i)) = 1;
% end
% 
% c1 = clustering1(best_indices1);
% c2 = clustering2(best_indices2);

end


function [lcs_length, indices1, indices2] = findLCS(seq1, seq2)
m = length(seq1);
n = length(seq2);

% Initialize LCS matrix
L = zeros(m+1, n+1);

% Fill LCS matrix
for i = 1:m
    for j = 1:n
        if seq1(i) == seq2(j)
            L(i+1, j+1) = L(i, j) + 1;
        else
            L(i+1, j+1) = max(L(i+1, j), L(i, j+1));
        end
    end
end

% Backtrack to find the indices of the LCS
lcs_length = L(m+1, n+1);
indices1 = zeros(1, lcs_length);
indices2 = zeros(1, lcs_length);
i = m;
j = n;
k = lcs_length;

while i > 0 && j > 0
    if seq1(i) == seq2(j)
        indices1(k) = i;
        indices2(k) = j;
        i = i - 1;
        j = j - 1;
        k = k - 1;
    elseif L(i, j+1) > L(i+1, j)
        i = i - 1;
    else
        j = j - 1;
    end
end
end

function perms = generate_permutations(arr)
% Generate all permutations of an array using recursion
n = length(arr);
perms = zeros(factorial(n), n);

if n == 1
    perms(1, 1) = arr;
else
    index = 1;
    for i = 1:n
        % Fix the i-th element and permute the rest
        rest = [arr(1:i-1), arr(i+1:end)];
        restPerms = generate_permutations(rest);

        % Add the fixed element to the front of each permutation of the rest
        for j = 1:size(restPerms, 1)
            perms(index, :) = [arr(i), restPerms(j, :)];
            index = index + 1;
        end
    end
end
end
