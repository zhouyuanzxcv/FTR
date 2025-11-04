function show_longitudinal_stability_in_heatmap(result_mat)
result_mat1 = result_mat * 255/256;
cm = parula(256);
cm(end+1,:) = [0.7 0.7 0.7];
% Create the heatmap
h = imagesc(result_mat1);
colormap(cm);
clim([0, 1]);  % Set color limits (adjust if needed)

% Set NaN values to display as gray
nan_mask = isnan(result_mat1);  % Mask of NaN locations
hold on;
h_nan = imagesc(nan_mask);     % Overlay NaN regions
set(h_nan, 'AlphaData', nan_mask); % Make only NaN regions visible

colormap(gca, cm); % Add gray to colormap
hold off;

% Axis settings
set(gca, 'YDir', 'reverse');
axis equal;
xlim([0.5 4.5]);
ylim([0.5 8.5]);



% Add text labels (skip NaNs)
textStrings = num2str(result_mat(:), '%.2f');  % Convert to strings
textStrings = strtrim(cellstr(textStrings));
textStrings(strcmp(textStrings, 'NaN')) = {''}; % Remove 'NaN' strings
[x, y] = meshgrid(1:4, 1:8);
valid = ~isnan(result_mat(:));  % Indices of non-NaN values
text(x(valid), y(valid), textStrings(valid), 'HorizontalAlignment', 'center');

end
