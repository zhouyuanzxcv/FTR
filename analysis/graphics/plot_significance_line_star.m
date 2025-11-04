function plot_significance_line_star(line_x,line_y, p1, spacing)
%PLOT_SIGNIFICANCE_LINE_STAR Summary of this function goes here 
%   Detailed explanation goes here
line(line_x, line_y,'linewidth',0.5,'color',[0,0,0]);
text(mean(line_x), mean(line_y)+spacing, convert_ps2stars(p1), ...
    'fontsize',14,'HorizontalAlignment','center');
ylim1 = ylim;
length = ylim1(2) - ylim1(1);
line([line_x(1),line_x(1)], [line_y(1),line_y(1)-length/(8*5)],'color',[0,0,0]);
line([line_x(2),line_x(2)], [line_y(2),line_y(2)-length/(8*5)],'color',[0,0,0]);

end

