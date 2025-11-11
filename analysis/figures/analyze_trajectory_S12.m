function [outputArg1,outputArg2] = analyze_trajectory_S12(inputArg1,inputArg2)
%ANALYZE_TRAJECTORY Summary of this function goes here
%   Detailed explanation goes here

close all


data_names = {'ADNI_FSX_HS','OASIS3_ROD1_HS','NACC_HS'};
method = {'FTR_MCEM'};

nsubtype = 3;


results = load_data_result(data_names, method, nsubtype, 0);
for sel_idx = 1:length(data_names)
    joindata = results{sel_idx,1}.joindata;
    biomarker_names = [results{sel_idx,1}.biomarker_names];
    % data = joindata{:, biomarker_names};
    traj = results{sel_idx,1}.mdl.re_traj;
    thres = 1;
    ord = cal_order(traj,thres,biomarker_names, data_names{sel_idx});
end

end

function ord = cal_order(traj,thres,bio_name, subfig_name)

[nbiom,num_int,nsubtype] = size(traj);
ord = zeros(nsubtype,nbiom);
s = (1:num_int)/num_int;


[colors, linestyles] = get_color_linestyle(bio_name);
for k = 1:nsubtype
    
    comp = sum(traj(:,:,k) < thres,2);
    [~,ord(k,:)] = sort(comp);
    f = figure();
    f.Position = get_figure_position('figS12');
    % subplot(3, 3, 3 * (idx_row - 1) + k)
    legendText = get_legend_texts(bio_name);

    for j = 1:nbiom
        plot(s,traj(ord(k,j),:,k),linestyles{j},'Color',colors(j,:),'LineWidth',1);
        hold on;
    end
    xlabel("Stage")
    ylabel("z-score (sigma)")

    plot([0,1],[thres,thres],':','Color','k','LineWidth',1);

    set(gca,"Position",[0.1300 0.1 0.7750 0.50])
    lgd_disp_num = 20;
    ord_k = ord(k, :);
    lgd = legend(legendText(ord_k(1:lgd_disp_num)),'Location', 'northoutside', 'Interpreter', 'none', 'NumColumns', 2);
    lgd.Position = [0.1300 0.7 0.7750 0.15];
    hold off;
    export_fig(['./figures/all/S12_',subfig_name,'_',int2str(k),'.jpg'],'-r500','-transparent')
    close all
end

end

function [colors, linestyles] = get_color_linestyle(bio_name)
linestyles = repmat({'-','--',':','-.'}, 1, ceil(length(bio_name) / 4));
linestyles = linestyles(1:length(bio_name));

colors = turbo(length(bio_name));
colors = flipud(colors);
end

function legendText = get_legend_texts(bio_name)
for j = 1:length(bio_name)
    if isempty(bio_name)
        legendText{j} = sprintf(int2str(j-1) , j);
    else        
        legendText{j} = sprintf(char(bio_name(j)), j);
    end
end

end
