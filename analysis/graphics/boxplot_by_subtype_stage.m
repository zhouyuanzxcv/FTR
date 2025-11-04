function hs = boxplot_by_subtype_stage(slopes_all,stage_names,yname)
%BOXPLOT_BY_SUBTYPE_STAGE Summary of this function goes here
%   Detailed explanation goes here
% Input:
%   slopes_all - num_data x num_subtype cell array, where each cell
%                contains a vector of numbers
%   stage_names - 1 x num_data cell array, each cell contains a string
%                 indicating the name
%   yname - a string or a cell array to be input into ylabel

nsubtype = size(slopes_all,2);

subtype_offsets = [-floor(nsubtype/2):floor(nsubtype/2)];
subtype_name = cellfun(@(x) ['Subtype ', num2str(x)], num2cell(1:nsubtype), 'UniformOutput', 0);

hold on;

for j = 1:length(stage_names)
    for k = 1:nsubtype
        y = slopes_all{j,k};
        
        stage_width = (nsubtype+1);
        x = (j-1)*stage_width + 1 + subtype_offsets(k);
        x = repmat(x, size(y,1), 1);
        
        swarmchart(x, y, 5, get_subtype_color(k), 'filled', 'MarkerFaceAlpha', 0.4);
        hs(k) = boxchart(x, y, 'boxfacecolor', ...
            get_subtype_color(k)/2, 'whiskerlinecolor', get_subtype_color(k)/2, ...
            'markerstyle', 'none');
        
        
    end
end

% legend(hs, subtype_name);

xticks((0:nsubtype-1)*stage_width + 1);
xticklabels(stage_names);
ylabel(yname);
end

