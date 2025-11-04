function plot_stage_distribution_by_diagnosis(data)
stage = data.stage;
stag_CN = stage(data.diagnosis == 0);
stag_MCI = stage(data.diagnosis == 0.5);
stag_AD = stage(data.diagnosis == 1);
binrng = 0:0.1:1;
counts(1,:) = histcounts(stag_CN, binrng);
counts(2,:) = histcounts(stag_MCI, binrng);
counts(3,:) = histcounts(stag_AD, binrng);

h = bar((binrng(2:end) + binrng(1:end-1))/2, counts);
set(h(1),'FaceColor',[0.2 0.2 0.8]);
set(h(2),'FaceColor',[0.2 0.8 0.2]);
set(h(3),'FaceColor',[0.8 0.2 0.2]);
legend('CN','MCI','AD')
xlabel({'Stage'});
ylabel({'Number of data points'});

end

function plot_stage_distribution_by_diagnosis1(data)
stage = data.stage;

% stag_MCI = stage(data.diagnosis == 0.5);
% stag_AD = stage(data.diagnosis == 1);
binrng = 0:0.1:1;
width = 0.18;

colors = [0.4, 0.4, 0.4];
colors = cat(1, colors, get_subtype_color([1,2,3]));

dxs = [0, 0.5, 1];

pos_offset = [-0.03, 0, 0.03];

edgecolors = [0.2 0.2 0.8; 0.2 0.8 0.2; 0.8 0.2 0.2];

for i = 1:length(dxs)

    stag_CN0 = stage(data.diagnosis == dxs(i) & data.group == 0);
    stag_CN1 = stage(data.diagnosis == dxs(i) & data.group == 1 & data.subtype == 1);
    stag_CN2 = stage(data.diagnosis == dxs(i) & data.group == 1 & data.subtype == 2);
    stag_CN3 = stage(data.diagnosis == dxs(i) & data.group == 1 & data.subtype == 3);
    counts_CN(1,:) = histcounts(stag_CN0, binrng);
    counts_CN(2,:) = histcounts(stag_CN1, binrng);
    counts_CN(3,:) = histcounts(stag_CN2, binrng);
    counts_CN(4,:) = histcounts(stag_CN3, binrng);

    h1 = bar((binrng(2:end) + binrng(1:end-1))/2 + pos_offset(i), counts_CN, width, ...
        'stacked','facecolor','flat');
    hold on;

    for k = 1:4
        h1(k).CData = colors(k,:);
    end

    set(h1, 'EdgeColor', edgecolors(i,:), 'linewidth', 2);
end

% counts2 = histcounts(stag_MCI, binrng);
% h2 = bar((binrng(2:end) + binrng(1:end-1))/2, counts2, width);
% 
% counts3 = histcounts(stag_AD, binrng);
% h3 = bar((binrng(2:end) + binrng(1:end-1))/2 + 0.03, counts3, width);



% set(h1,'EdgeColor',[0.2 0.2 0.8],'LineWidth',1);
% set(h2,'EdgeColor',[0.2 0.8 0.2],'LineWidth',1);
% set(h3,'EdgeColor',[0.8 0.2 0.2],'LineWidth',1);
xlim([0,1]);

legend('CN','MCI','AD')
xlabel({'Stage'});
ylabel({'Number of data points'});

end
