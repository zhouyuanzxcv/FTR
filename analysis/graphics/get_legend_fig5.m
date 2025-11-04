function get_legend_fig5()

close all
colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], ...
    [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560], ...
    [0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330], [0.6350 0.0780 0.1840]};
shapes = {'o','square','diamond','v','^','>','<'};
method_labels = {'FTR (MCEM)','SuStaIn','Hierarchical Clustering','CTV-HV-Ratio','DEBM','EBM','SVM'};

% Create figure
f = figure;
f.Position = get_figure_position(3);
hold on

% Create empty plots for legend
handles = [];
for methodIdx = 1:length(shapes)
    h = plot(nan, nan, shapes{methodIdx}, 'Color', colors{methodIdx}, ...
             'MarkerFaceColor', colors{methodIdx}, 'MarkerEdgeColor', colors{methodIdx});
    handles = [handles, h];
end
hold off

% Hide axes
axis off;

% Create legend
legend(method_labels, 'Location', 'westoutside', 'Interpreter', 'none', ...
       'Orientation', 'horizontal', 'NumColumns', 7);

% Export the figure with just the legend
export_fig './figures/all/5_legend.jpg' -r500 -transparent

end
