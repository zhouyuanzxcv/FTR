function plot_stage_vs_MMSE(data, nsubtype)

stage = data.stage;
mmse = data.MMSE;

mmse1 = mmse + (rand(size(data.MMSE))-0.5);

% scatter_with_correlation(stage, mmse, 'All stages');

% Fit linear regression
x = stage;
y = mmse;

mdl = fitlm(x, y);

% Plot the linear regression line
x_range = linspace(min(x), max(x), 100);
y_fit = predict(mdl, x_range');
h_all = plot(x_range, y_fit, 'Color', [0, 0, 0], 'LineWidth', 1.5);

hold on;

% Remove rows with NaN in MMSE and stage columns
valid_indices = ~isnan(x) & ~isnan(y);
cleaned_x = x(valid_indices, :);
cleaned_y = y(valid_indices, :);

% Calculate correlation coefficient (r-value) and p-value
[r_all, p] = corr(cleaned_x, cleaned_y);

% % Add subtype names, r-value, and p-value to the legend
% legend_info = sprintf('%s (r = %.3f, p = %.3e)', ['Subtype ',num2str(k)], r, p);
% legend({legend_info});

xlabel('Stage');
ylabel('MMSE');


hs = [];

% for subject = 1:size(stage, 1)
%     k = data.subtype(subject);
%     x = stage(subject);
%     y = mmse(subject);
%     scatter(x + rand(size(x))/100, y, 20, get_subtype_color(k), ...
%         'filled', 'MarkerFaceAlpha', 0.35);
% end


control_idx = data.group == 0;
case_idx = data.group == 1;

colors = get_subtype_color(data.subtype);

colors(control_idx,:) = repmat([0.4, 0.4, 0.4], sum(control_idx), 1);

scatter(stage+rand(size(stage))/100, mmse1, 20, colors, 'filled', 'MarkerFaceAlpha', 0.35);
% 
% for k = 1:nsubtype
%     % Fit linear regression
%     x = stage(data.subtype==k);
%     y = mmse(data.subtype==k);
% 
%     scatter(x + rand(size(x))/100, y, 20, get_subtype_color(k), ...
%         'filled', 'MarkerFaceAlpha', 0.35);
% 
% end

rs = [];
for k = 1:nsubtype
% Fit linear regression
    x = stage(data.subtype==k & case_idx);
    y = mmse(data.subtype==k & case_idx);

    mdl = fitlm(x, y);
    
    % Plot the linear regression line
    x_range = linspace(min(x), max(x), 100);
    y_fit = predict(mdl, x_range');
    h = plot(x_range, y_fit, 'Color', get_subtype_color(k), 'LineWidth', 2);
    hs = [hs, h];
    
    % Remove rows with NaN in MMSE and stage columns
    valid_indices = ~isnan(x) & ~isnan(y);
    cleaned_x = x(valid_indices, :);
    cleaned_y = y(valid_indices, :);
    
    % Calculate correlation coefficient (r-value) and p-value
    [r, p] = corr(cleaned_x, cleaned_y);
    rs(end+1) = r;
    
    % Add subtype names, r-value, and p-value to the legend
    % legend_info = sprintf('%s (r = %.3f, p = %.3e)', ['Subtype ',num2str(k)], r, p);
    % legend({legend_info});
end
xlim([0 1])
ylim([0, 30.5])

legend_str = {sprintf('All $r=%.2f$\n', r_all), sprintf('S1 $r=%.2f$\n',rs(1)), ...
    sprintf('S2 $r=%.2f$\n',rs(2)), sprintf('S3 $r=%.2f$\n',rs(3))};
legend([h_all,hs],legend_str, 'Interpreter','latex','location','best');
end
