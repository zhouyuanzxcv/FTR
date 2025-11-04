function sustain_traj1 = get_sustain_traj(joindata, biomarker_names, Z_vals)
subtype = joindata.subtype;
stage = joindata.stage;
nsubtype = length(unique(subtype));

% sustain has stages in the range (0, 1, ..., length(Z_vals)*B + 1)
stages = (0: length(Z_vals)*length(biomarker_names)+1);
num_stage = length(stages);

if 1 % moving average
    stage_unique = unique(stage);

    sustain_traj = nan(length(biomarker_names), num_stage, nsubtype);
    for k = 1:nsubtype
        for j = 1:length(stage_unique)
            s = stage_unique(j);
            bioms = mean(joindata{subtype == k & stage == s, biomarker_names}, 1);
            idx_stage = stage_unique(j) == stages;
            sustain_traj(:,idx_stage,k) = bioms;
        end
    end

    sustain_traj1 = sustain_traj;
    for k = 1:nsubtype
        sustain_traj1(:,:,k) = move_average(sustain_traj(:,:,k), 3:2:9);
    end
else % gaussian process
    sustain_traj1 = fit_GP_by_subtype_stage(joindata, stages, biomarker_names);
end

end


function y_smooth = move_average(y, windowSizes)
% Sample data
% x = linspace(0, 10, 100);
% y = sin(x) + 0.4*randn(size(x));

% Possible window sizes (odd numbers)
% windowSizes = 3:2:7;
mse = zeros(size(windowSizes));

% Leave-one-out cross validation
for i = 1:length(windowSizes)
    ws = windowSizes(i);
    pred = zeros(size(y));
    
    for j = 1:size(y,2)
        % Create indices for window (handling edges)
        idx = max(1, j-floor(ws/2)):min(size(y,2), j+floor(ws/2));
        idx = idx(idx~=j); % Exclude current point
        
        % Simple average prediction
        pred(:,j) = mean(y(:,idx), 2);
    end
    
    mse(i) = mean((y - pred).^2, 'all','omitnan');
end

% Find window with minimum MSE
[~, bestIdx] = min(mse);
optimalWindow = windowSizes(bestIdx);

fprintf('Optimal window size %d\n', optimalWindow);

% Apply optimal moving average
y_smooth = movmean(y, optimalWindow, 2);

% % Plot results
% figure;
% subplot(2,1,1);
% plot(windowSizes, mse, 'bo-');
% xlabel('Window Size');
% ylabel('Mean Squared Error');
% title('Cross-Validation MSE vs. Window Size');
% 
% subplot(2,1,2);
% plot(x, y, 'b.', 'MarkerSize', 10);
% hold on;
% plot(x, y_smooth, 'r-', 'LineWidth', 2);
% title(['Optimal Moving Average (Window = ' num2str(optimalWindow) ')']);
% legend('Original', 'Smoothed');
end