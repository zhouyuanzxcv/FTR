function cutoff_value = get_csf_cutoff(data)
%GET_CSF_CUTOFF Summary of this function goes here
%   Detailed explanation goes here
close all;

data(isnan(data.ABETA),:) = [];
x = data.ABETA;

gmm = fitgmdist(x, 2, 'Replicates', 100);
mu = gmm.mu;
sigma = sqrt(gmm.Sigma);
pi = gmm.ComponentProportion;

binWidth = 30; % Adjust bin width as needed
minValue = min(x);
maxValue = max(x);
edges = (minValue:binWidth:maxValue)';

% centers = 0.5 * (edges(2:end) + edges(1:end-1));
centers = (minValue:maxValue)';

comp1 = normpdf(centers, mu(1), sigma(:,:,1));
comp1 = comp1 * pi(1);

comp2 = normpdf(centers, mu(2), sigma(:,:,2));
comp2 = comp2 * pi(2);

figure;
% tiledlayout(1, 2);

% nexttile;
histogram(x, edges, 'Normalization', 'pdf', 'FaceColor', [0.5,0.5,0.5], ...
    'EdgeColor', 'none');
xlabel('CSF A\beta_{1-42} (pg/ml)');
ylabel('Probability density');
title({'NACC CSF A\beta_{1-42} distribution',''},'fontweight','normal');
% histogram(x, edges, 'Normalization', 'pdf', 'FaceColor', 'b', ...
%     'EdgeColor', 'none', 'DisplayName', 'All');
hold on;
h_c1 = plot(centers, comp1, 'r--', 'linewidth', 1);
plot(centers, comp2, 'r--', 'linewidth', 1);
h_gmm = plot(centers, comp1 + comp2, 'r-', 'linewidth', 1);

cutoff_value = 747;

h_line = line([cutoff_value, cutoff_value], [0, 0.0015], ...
    'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.5 ...
    );
legend([h_c1,h_gmm,h_line], {'Fitted GMM component','Fitted GMM',['Cutoff = ',num2str(cutoff_value)]});


if ismember('diagnosis', data.Properties.VariableNames)
    % Extract ABETA values for diagnosis == 0 and diagnosis == 1
    ABETA_CN = data(data.diagnosis == 0, 'ABETA').ABETA;
    ABETA_MCI = data(data.diagnosis == 0.5, 'ABETA').ABETA;
    ABETA_AD = data(data.diagnosis == 1, 'ABETA').ABETA;

    % Define histogram properties
    % binWidth = 50; % Adjust bin width as needed
    % minValue = min([min(ABETA_CN), min(ABETA_AD)]);
    % maxValue = max([max(ABETA_CN), max(ABETA_AD)]);
    % edges = minValue:binWidth:maxValue;

    % Plot histograms with different colors
    nexttile;
    histogram(ABETA_CN, edges, 'FaceColor', 'b', 'EdgeColor', 'none', 'DisplayName', 'CN');
    hold on;
    histogram(ABETA_MCI, edges, 'FaceColor', 'g', 'EdgeColor', 'none', 'DisplayName', 'MCI');
    histogram(ABETA_AD, edges, 'FaceColor', 'r', 'EdgeColor', 'none', 'DisplayName', 'AD');
    line_cutoff_y = max([histcounts(ABETA_CN, edges), ...
        histcounts(ABETA_AD, edges), histcounts(ABETA_MCI, edges)]);


    % scale = (length(ABETA_CN) + length(ABETA_MCI) + length(ABETA_AD));
    % comp1_count = comp1 ./ sum(comp1) * scale;
    % comp2_count = comp2 ./ sum(comp2) * scale;
    % plot(centers, comp1_count, '-k', 'linewidth', 1);
    % plot(centers, comp2_count, '-k', 'linewidth', 1);
    % plot(centers, comp1_count + comp2_count, '-k', 'linewidth', 1);


    % Add vertical line for cutoff at ABETA = 880
    % cutoff_value = 880;
    line([cutoff_value, cutoff_value], [0, line_cutoff_y], ...
        'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.5, 'DisplayName', ['Cutoff = ',num2str(cutoff_value)]);

    xlabel('CSF A\beta_{1-42} (pg/ml)');
    ylabel('Count');
    title({'NACC CSF A\beta_{1-42} distribution',''},'fontweight','normal');
    % title(['Histogram of CSF ABETA Values for CN, MCI, AD with Cutoff at ABETA = ',num2str(cutoff_value)]);
    legend;
    hold off;
end

% export_fig('../../figures/ADNI_CSF_cutoff.jpg','-transparent','-r500','-nocrop');

end

