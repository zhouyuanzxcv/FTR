function [h_box, h_swarm] = boxplot_w_dots(x, y, mk_sz, color, box_width)
%BOXPLOT_W_DOTS Summary of this function goes here
%   Detailed explanation goes here
if length(y) ~= length(x)
    x = repmat(x, size(y));
end

h_swarm = swarmchart(x, y, mk_sz, color, 'filled','XJitterWidth', 0.9 * (box_width / 0.5));
h_box = boxchart(x, y, 'boxfacecolor', ...
    color, 'whiskerlinecolor', color, ...
    'markerstyle', 'none', 'BoxWidth', box_width);
end

