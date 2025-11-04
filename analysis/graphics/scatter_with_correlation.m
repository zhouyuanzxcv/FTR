function [outputArg1,outputArg2] = scatter_with_correlation(x,y,title, jitter_y)
%SCATTER_WITH_ Summary of this function goes here
%   Detailed explanation goes here
if nargin < 4
    jitter_y = 0;
end

dot_color = [0,0,0];
alphaValue = 1;
line_color = [0,0,0];
markerSize = 5;

if jitter_y
    y_max = max(abs(y));
    y1 = y + (rand(size(y))-0.5)/100*y_max;
else
    y1 = y;
end
scatter_handle = scatter(x, y1, markerSize, 'MarkerFaceColor', dot_color, 'MarkerEdgeColor', 'none');
scatter_handle.MarkerFaceAlpha = alphaValue; % Set opacity for the marker

xlabel('First argument');
ylabel('Second argument');

hold on;

% Fit linear regression
mdl = fitlm(x, y);

% Plot the linear regression line
x_range = linspace(min(x), max(x), 100);
y_fit = predict(mdl, x_range');
plot(x_range, y_fit, 'Color', line_color, 'LineWidth', 1.5);

% Remove rows with NaN in MMSE and stage columns
valid_indices = ~isnan(x) & ~isnan(y);
cleaned_x = x(valid_indices, :);
cleaned_y = y(valid_indices, :);

% Calculate correlation coefficient (r-value) and p-value
[r, p] = corr(cleaned_x, cleaned_y);

% Add subtype names, r-value, and p-value to the legend
legend_info = sprintf('%s (r = %.3f, p = %.3e)', title, r, p);
legend({legend_info});


end

