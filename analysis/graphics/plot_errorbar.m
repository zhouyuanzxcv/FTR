    % cd 'F:\Research\Project\1_Filtered_trajectory_recovery'
addpath("FTR_code\utils")
addpath("FTR_code")
dir = '\Result\Multiple Trajectory\';

figure;
noise_list = [0.01,0.1,0.2];

name = 'Sigmoid';
result_RI = readmatrix(strcat(dir,name,'_RI.xlsx'));

for i = 1:length(noise_list)

    options = [];
    handle = parse_param(options, 'figure', gcf);
    alpha = parse_param(options, 'alpha', 0.5);
    line_width = parse_param(options, 'line_width', 2);
    center_type = parse_param(options, 'center_type', 'mean');
    error_type = parse_param(options, 'error_type', 'std');
    plot_type = parse_param(options, 'plot_type', 'plot');

    m = 6*(i-1);
    subplot(3,3,3+i);

    x1 = 5:5:30;
    x1 = x1';
    options.color_area = [70 240 240]./255;    % Blue theme
    options.color_line = [0 130 200]./255;
    y1 = result_RI(m+1:m+6,16:20);
    plot_uncertainty(x1, y1, options)
    hold on;

    x2 = 5:5:30;
    x2 = x2';
    options.color_area = [250 190 212]./255;    % Red theme
    options.color_line = [230 25 75]./255;
    y2 = result_RI(m+1:m+6,9:13);
    plot_uncertainty(x2, y2, options)
    hold on;    

    x3 = 5:5:15;
    x3 = x3';
    options.color_area = [220 190 255]./255;    % Purple theme
    options.color_line = [145 30 180]./255;
    y3 = result_RI(m+1:m+3,2:6);
    plot_uncertainty(x3, y3, options)
    hold on;

    plot(x3,mean(y3,2),':|','Color',[145 30 180]./255,'LineWidth',1,'MarkerSize',3)
    hold on;
    plot(x1,mean(y1,2),':x','Color',[0 130 200]./255,'LineWidth',1,'MarkerSize',3)
    hold on;
    plot(x2,mean(y2,2),'Color',[230 25 75]./255,'LineWidth',1,'MarkerSize',3)
    hold on;

    xlim([5 30])
    xlabel('Number of biomarkers');
    if i == 1
        ylabel({'Sigmoid function','Rand index'});
    else
        ylabel('Rand index');
    end
    if i == 3
        ylim([0.6 1])
    end
    title({['Noise = ',num2str(noise_list(i))],''})
end

name = 'SuStaIn';
result_RI = readmatrix(strcat(dir,name,'_RI.xlsx'));

for i = 1:length(noise_list)

    options = [];
    handle = parse_param(options, 'figure', gcf);
    alpha = parse_param(options, 'alpha', 0.5);
    line_width = parse_param(options, 'line_width', 2);
    center_type = parse_param(options, 'center_type', 'mean');
    error_type = parse_param(options, 'error_type', 'std');
    plot_type = parse_param(options, 'plot_type', 'plot');

    m = 6*(i-1);
    subplot(3,3,6+i);

    x1 = 5:5:30;
    x1 = x1';
    options.color_area = [70 240 240]./255;    % Blue theme
    options.color_line = [0 130 200]./255;
    y1 = result_RI(m+1:m+6,16:20);
    plot_uncertainty(x1, y1, options)
    hold on;

    x2 = 5:5:30;
    x2 = x2';
    options.color_area = [250 190 212]./255;    % Red theme
    options.color_line = [230 25 75]./255;
    y2 = result_RI(m+1:m+6,9:13);
    plot_uncertainty(x2, y2, options)
    hold on;    

    x3 = 5:5:15;
    x3 = x3';
    options.color_area = [220 190 255]./255;    % Purple theme
    options.color_line = [145 30 180]./255;
    y3 = result_RI(m+1:m+3,2:6);
    plot_uncertainty(x3, y3, options)
    hold on;


    plot(x3,mean(y3,2),':|','Color',[145 30 180]./255,'LineWidth',1,'MarkerSize',3)
    hold on;
    plot(x1,mean(y1,2),':x','Color',[0 130 200]./255,'LineWidth',1,'MarkerSize',3)
    hold on;
    plot(x2,mean(y2,2),'Color',[230 25 75]./255,'LineWidth',1,'MarkerSize',3)
    hold on;

    xlim([5 30])
    xlabel('Number of biomarkers');
    if i == 1
        ylabel({'Event permutation','Rand index'});
    else
        ylabel('Rand index');
    end
    if i == 3
        ylim([0.6 1])
    end
    title({['Noise = ',num2str(noise_list(i))],''})
end

legend({'','','','SuStaIn','FTR  w.o.f.','FTR',},'Location','bestoutside','Orientation','horizontal')

%export_fig('F:\Research\Project\1_Filtered_trajectory_recovery\multiple_errorbar.jpg','-transparent','-r500')