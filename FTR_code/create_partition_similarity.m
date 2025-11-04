function ari_matrix = create_partition_similarity(partition_matrix)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
n_partitions = size(partition_matrix, 2);
ari_matrix = zeros(n_partitions);

for i = 1:n_partitions
    for j = i:n_partitions
        ari_matrix(i,j) = rand_index(partition_matrix(:,i), partition_matrix(:,j));
        ari_matrix(j,i) = ari_matrix(i,j); % Symmetric
    end
end
end

