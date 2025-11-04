function ord = cal_order(traj,thres,bio_name,save_dir)
if nargin < 4
    save_dir = [];
end

[nbiom,num_int,nsubtype] = size(traj);
ord = zeros(nsubtype,nbiom);
s = linspace(0, 1, num_int);


[colors, linestyles] = get_color_linestyle(bio_name);
for k = 1:nsubtype
    comp = sum(traj(:,:,k) < thres,2);
    [~,ord(k,:)] = sort(comp);

    h = figure;

    legendText = get_legend_texts(bio_name);

    for j = 1:nbiom
        plot(s,traj(ord(k,j),:,k),linestyles{j},'Color',colors(j,:),'LineWidth',1);
        hold on;
    end
    xlabel("Stage")
    ylabel("z-score")
    
    plot([0,1],[thres,thres],':','Color','k','LineWidth',1);
    %lgd = legend(legendText,'Location', 'northwest');
    lgd = legend(legendText(ord(k,:)),'Location', 'eastoutside', 'Interpreter', 'none');
    hold off;
    if ~isempty(save_dir)
        saveas(h,[save_dir,'/trajectory_subtype',int2str(k),'.png']);
    end
end

subtype = cell(nsubtype,1);
for i = 1:nsubtype
    subtype(i) = cellstr(['Subtype',int2str(i)]);
    row_table = cell2table(subtype);
end
ord_table = array2table(ord);
ord_name_table = array2table(bio_name(ord));

if ~isempty(save_dir)
    writetable([row_table,ord_table], [save_dir,'/biomarker_order.csv'])
    writetable([row_table,ord_name_table], [save_dir,'/biomarker_name_order.csv'])
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
        % YZ: The second argument j in sprintf may not be needed
        legendText{j} = sprintf(char(bio_name(j)), j);
    end
end

end