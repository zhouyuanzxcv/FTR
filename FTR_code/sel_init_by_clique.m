function [best_partition, best_idx]  = sel_init_by_clique(partition_matrix, ari_threshold)
% Inputs:
%   partition_matrix: matrix where each column is a partition vector
%   ari_threshold: similarity threshold for binarizing the ARI matrix
% Output:
%   best_partition: the selected representative partition vector

%% Step 1: Compute ARI similarity matrix
ari_matrix = create_partition_similarity(partition_matrix);

%% Step 2: Threshold to create binary similarity matrix
binary_matrix = ari_matrix >= ari_threshold;

%% Step 3: Find all maximal cliques (size ≥ 2)
G = graph(binary_matrix);
bins = conncomp(G);

% Count members in each component and filter
unique_bins = unique(bins);
% n_cliques = 0;
clique_indices = cell(length(unique_bins), 1);
clique_sizes = zeros(length(unique_bins), 1);

for i = 1:length(unique_bins)
    current_clique = find(bins == unique_bins(i));
    clique_sizes(i) = length(current_clique);
    clique_indices{i} = current_clique;
end

% clique_size_thresh = round(max(clique_sizes)/2);
[sorted_size, sort_ind] = sort(clique_sizes, 'descend');

% keep_inds = clique_sizes > clique_size_thresh;
keep_inds = sort_ind(1:5);

clique_indices = clique_indices(keep_inds);
clique_sizes = clique_sizes(keep_inds);

n_cliques = length(clique_indices);

%% Step 4: Find the center of each clique (based on in-clique similarity)
clique_centers = zeros(n_cliques, 1);
for i = 1:n_cliques
    current_clique = clique_indices{i};
    if isempty(current_clique)
        continue;
    end

    % Compute centrality scores (sum of ARIs to others IN THE CLIQUE)
%     centrality_scores = zeros(length(current_clique), 1);
%     for j = 1:length(current_clique)
%         idx = current_clique(j);
%         centrality_scores(j) = sum(ari_matrix(idx, current_clique));
%     end
    centrality_scores = sum(ari_matrix(current_clique, current_clique), 2);

    % Find the most central partition WITHIN THE CLIQUE
    [~, max_idx] = max(centrality_scores);
    clique_centers(i) = current_clique(max_idx);
end

%% Step 5: Select the best partition (based on global similarity)
% Compute global centrality for each clique center
% global_centrality = zeros(n_cliques, 1);
% for i = 1:n_cliques
%     if clique_centers(i) == 0
%         global_centrality(i) = -Inf;
%         continue;
%     end
%     % Sum of ARIs between this center and ALL OTHER PARTITIONS
%     global_centrality(i) = sum(ari_matrix(clique_centers(i), :));
% end

ari_matrix1 = sort(ari_matrix(clique_centers, :), 2, 'descend');
global_centrality = sum(ari_matrix1(:, 1:round(end/2)), 2);

% Choose the center with highest global centrality
[~, best_clique_idx] = max(global_centrality);

best_idx = clique_centers(best_clique_idx);
best_partition = partition_matrix(:, best_idx);
end