function ord = show_trajectory(traj,thres,bio_name)

[nbiom,num_int,nsubtype] = size(traj);

if nargin < 3
    bio_name = cell(1, nbiom);
    for j = 1:nbiom
        bio_name{j} = ['biomarker ',int2str(j)];
    end
end

if nargin < 2
    thres = 1;
end

ord = zeros(nsubtype,nbiom);
s = linspace(0, 1, num_int);


[colors, linestyles] = get_color_linestyle(bio_name);
for k = 1:nsubtype
    compare1 = sum(traj(:,:,k) < thres,2);
    compare2 = mean(traj(:,:,k), 2);
    [~,ord(k,:)] = sortrows([compare1,compare2], {'ascend','descend'});

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
end

end

function [colors, linestyles] = get_color_linestyle(bio_name)
linestyles = repmat({'-','--',':','-.'}, 1, ceil(length(bio_name) / 4));
linestyles = linestyles(1:length(bio_name));
colors = colormap(turbo(length(bio_name)));
colors = flipud(colors);
end

function legendText = get_legend_texts(bio_name)
for j = 1:length(bio_name)
    if isempty(bio_name)
        legendText{j} = sprintf(int2str(j-1));
    else
        % YZ: The second argument j in sprintf may not be needed
        legendText{j} = sprintf(char(bio_name(j)));
    end
end

end