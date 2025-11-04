function show_trajectory_in_PCA(data, subtypes_gt, subtypes_est, re_traj)
%SHOW_TRAJECTORY_IN_PCA Summary of this function goes here
%   Detailed explanation goes here
% if nargin < 3
%     re_traj = [];
% end
% 
% if nargin < 4
%     fig_h = figure;
% end

if isempty(subtypes_gt)
    num_subtype = length(unique(subtypes_est));
else
    num_subtype = length(unique(subtypes_gt));
end

marker_types = {'^','square','o'};

traj_colors = cat(3, spring, summer, winter);

c = squeeze(traj_colors(128,:,:))';

[coeffs,data_pca,latent,tsquared,explained,mu] = pca(data, 'NumComponents', 3);

hold on;

traj_hs = [];
data_hs = [];
dot_sz = 8;
for k = 1:num_subtype

    if isempty(subtypes_gt)

        data_k = data(subtypes_est == k, :);
        data_pca = (data_k - repmat(mu, [size(data_k,1), 1])) * coeffs;
        data_hs(k) = scatter3(data_pca(:,1), data_pca(:,2), data_pca(:,3),dot_sz,c(k,:), ...
            'filled','o','MarkerFaceAlpha',0.5);
        
    else
        for k_gt = 1:num_subtype
            data_k = data(subtypes_est == k & subtypes_gt == k_gt,:);
            data_pca = (data_k - repmat(mu, [size(data_k,1), 1])) * coeffs;
            data_hs(k) = scatter3(data_pca(:,1), data_pca(:,2), data_pca(:,3),dot_sz,c(k,:), ...
                'filled',marker_types{k_gt},'MarkerFaceAlpha',0.5);
        end
    end
end

for k = 1:num_subtype
    if ~isempty(re_traj)
        traj_pca = (re_traj(:,:,k)' - repmat(mu, [size(re_traj(:,:,k),2), 1])) * coeffs;

        for j = 1:size(traj_pca, 1)-1
            traj_hs(k) = plot3(traj_pca(j:j+1,1), traj_pca(j:j+1,2), traj_pca(j:j+1,3), ...
                '-','Color',traj_colors(round(j/size(traj_pca,1)*256),:,k),'LineWidth',5);
        end
    end
end

subtype_names = cellfun(@(x) ['Subtype ',int2str(x),' (est.)'], num2cell((1:num_subtype)),'UniformOutput',false);
legend(data_hs, subtype_names);

end

