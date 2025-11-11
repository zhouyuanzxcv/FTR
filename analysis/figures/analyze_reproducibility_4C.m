function [outputArg1,outputArg2] = analyze_reproducibility_4C(inputArg1,inputArg2)
%ANALYZE_REPRODUCIBILITY Summary of this function goes here
%   Detailed explanation goes here
close all


datas = {{'OASIS3_ROD1_HS','ADNI_FSX_HS'}, {'NACC_HS','ADNI_FSX_HS'}, ...
    {'NACC_HS','OASIS3_ROD1_HS'}};
subplot_idx = [1,4,7];

labels = {{'OASIS','ADNI'},{'NACC','ADNI'},{'NACC','OASIS'}};


f = figure;
% f.Position = get_figure_position();

result_mats = [];

tiledlayout(3,3, "TileSpacing", "loose")
for dataIdx = 1:length(datas)

    data = datas{dataIdx};

    method = {'FTR_MCEM'};

    nsubtypes = [3;3];

    results = load_data_result(data, method, nsubtypes, 0, 1);

    biomarker_names = results{1,1}.biomarker_names;

    % Check that the two datasets have the same biomarkers
    same_bioms = false;
    if length(biomarker_names) == length(results{2,1}.biomarker_names)
        tmp = [strrep(biomarker_names', '-', '_'), results{2,1}.biomarker_names'];
        tmp = cellfun(@lower, tmp, 'UniformOutput', false);
        same_bioms = all(strcmp(tmp(:,1), tmp(:,2)));
    end
    if ~same_bioms
        disp('WARNING: The 2 datasets do not have the same set of biomarkers. Use common biomarkers');
        biomarker_names = intersect(results{1,1}.biomarker_names, results{2,1}.biomarker_names);
        inds_sel1 = ismember(results{1,1}.biomarker_names, biomarker_names);
        inds_sel2 = ismember(results{2,1}.biomarker_names, biomarker_names);
        inds_sel = {inds_sel1; inds_sel2};
    else
        inds_sel1 = true(size(biomarker_names));
        inds_sel2 = true(size(biomarker_names));
        inds_sel = {inds_sel1; inds_sel2};
    end



    [nbiom,num_int,~] = size(results{1,1}.mdl.re_traj);
    stage = (0:num_int-1)/(num_int-1);

    num_datasets = length(data);


    %% calculate trajectory similarity measured by PCC
    cef_mat = zeros(num_datasets,num_datasets,nsubtypes(1),nsubtypes(2));
    dist_mat = zeros(num_datasets,num_datasets);

    for i = 1:num_datasets
        for j = 1:num_datasets
            data1_traj = results{i,1}.mdl.re_traj;
            data2_traj = results{j,1}.mdl.re_traj;

            joindata1 = results{i,1}.joindata;
            joindata2 = results{j,1}.joindata;

            biomarker_names1 = results{i,1}.biomarker_names;
            biomarker_names2 = results{j,1}.biomarker_names;

            % select common biomarkers
            data1_traj = data1_traj(inds_sel{i},:,:);
            data2_traj = data2_traj(inds_sel{j},:,:);

            biomarker_names1 = biomarker_names;
            biomarker_names2 = biomarker_names;

            dist_all = abs(data1_traj - data2_traj);
            dist_mat(i,j) = mean(dist_all(:));

            for k1 = 1:nsubtypes(i)
                for k2 = 1:nsubtypes(j)
                    traj_k1 = squeeze(data1_traj(:,:,k1));
                    traj_k2 = squeeze(data2_traj(:,:,k2));

                    %% calculate MAD

                    mads(i,j,k1,k2) = mean(abs(traj_k1(:) - traj_k2(:)));

                    % calculate trajectory PCC
                    cef = corrcoef(traj_k1(:),traj_k2(:));
                    %             cef = corrcoef(traj_mid(k,:,i),traj_mid(k,:,j));

                    cef_mat(i,j,k1,k2) = cef(1,2);

                    %% calculate order correlation
                    thresh = 1;
                    [~,reach_stages1] = calc_order(traj_k1, stage, thresh);
                    [~,reach_stages2] = calc_order(traj_k2, stage, thresh);

                    kcc = corr(reach_stages1, reach_stages2, 'Type','Kendall');
                    kcc_mat(i,j,k1,k2) = kcc;

                    %% calculate mean atrophy correlation
                    mean_atrophy1 = mean(joindata1{joindata1.subtype == k1, biomarker_names1}, 1);
                    mean_atrophy2 = mean(joindata2{joindata2.subtype == k2, biomarker_names2}, 1);

                    cef_mean = corrcoef(mean_atrophy1, mean_atrophy2);

                    cef_mean_mat(i,j,k1,k2) = cef_mean(1,2);

                end
            end
        end
    end

    result_mat = squeeze(cef_mat(1,2,:,:));
    result_mats(:,:,dataIdx) = result_mat;
end
    
for dataIdx = 1:length(datas)
    % subplot(1,2,dataIdx)
    nexttile(subplot_idx(dataIdx));
    result_mat = result_mats(:,:,dataIdx);

    imagesc(result_mat);
    colormap('parula');
    clim([min(result_mat(:)), max(result_mat(:))]);

    axis equal
    xlim([0.5 3.5])
    ylim([0.5 3.5])

    set(gca, 'XTick', 1:3, 'XTickLabel', {'S1','S2','S3'}, 'YTick', 1:3, 'YTickLabel', {'S1','S2','S3'});
    textStrings = num2str(result_mat(:), '%.2f');
    textStrings = strtrim(cellstr(textStrings));
    [x, y] = meshgrid(1:3, 1:3);
    text(x(:), y(:), textStrings(:), 'HorizontalAlignment', 'center');
    ylabel(labels{dataIdx}{1});
    xlabel(labels{dataIdx}{2});
end


export_fig './figures/all/4C.jpg' -r500 -transparent
end

