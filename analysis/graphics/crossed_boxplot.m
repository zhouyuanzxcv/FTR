function crossed_boxplot(stabCell, timeCell, datasetNames, methodNames)
% CROSSED_BOXPLOT  Crossed 2D boxplots: vertical stability box + runtime dot±std
%
% Usage:
%   crossed_boxplot(stabCell, timeCell)
%   crossed_boxplot(stabCell, timeCell, datasetNames, methodNames)
%
% Inputs:
%   stabCell (6x2) : each cell either a numeric vector (stability) or []
%   timeCell (6x2) : each cell either a numeric vector (runtime, positive) or []
%   datasetNames    : (optional) cellstr of 6 dataset names
%   methodNames     : (optional) cellstr of length 2 method names
%
% Notes:
% - runtime is log10-transformed internally so box widths are consistent.
% - function is split into small helper functions; change any helper without touching others.

% -------------- Defaults & parameters --------------
if nargin < 3 || isempty(datasetNames)
    datasetNames = arrayfun(@(k) sprintf('Dataset %d',k), 1:size(stabCell,1), 'UniformOutput', false);
end
if nargin < 4 || isempty(methodNames)
    methodNames = { 'Method 1', 'Method 2' };
end

colors = lines(2);
markers = {'o','s'};            % marker style per method
boxWidth = 0.03;               % width in log10(runtime) units
labelFontSize = 8;

% -------------- Validate and gather stats for all items --------------
stats = collect_stats(stabCell, timeCell); % returns struct array for each valid pair
if isempty(stats)
    warning('No valid data to plot.'); return;
end

% check runtimes positive
if any([stats.runVals] <= 0)
    error('All runtime values must be positive for log transform.');
end

% compute global limits in log-space (x) and stability-space (y)
allRunLog = [stats.runLogVals];
allStab = [stats.stabVals];
allRunLog = allRunLog(:);
allStab = allStab(:);
xlimVals = compute_xlim(allRunLog);
ylimVals = compute_ylim(allStab);

% -------------- Precompute all obstacle rectangles (box + runtime bar) --------------
% We'll compute bpRect and barRect for every item BEFORE drawing anything.
obstacles = []; % Mx4 rectangles [x y w h] in data units (x = log10(runtime))
for k = 1:numel(stats)
    s = stats(k);
    % vertical boxplot rectangle
    bpRect = [ s.runMean - boxWidth/2, s.stabMin, boxWidth, (s.stabMax - s.stabMin) ];
    stats(k).bpRect = bpRect;
    obstacles = [obstacles; bpRect];

    % runtime bar rectangle (thin height)
    barH = 0.02 * range(ylimVals); % small vertical thickness for obstacle
    barRect = [ s.runMean - s.runStd, s.stabMed - barH/2, (2*s.runStd), barH ];
    stats(k).barRect = barRect;
    obstacles = [obstacles; barRect];
end

% -------------- Draw everything (boxplots first, then runtime markers) --------------
ax = gca(); hold(ax,'on'); grid(ax,'on');
xlim(ax, xlimVals); ylim(ax, ylimVals);
xlabel(ax, 'Running time (min)'); ylabel(ax, 'Algorithmic stability (pairwise ARI)');
%title(ax, 'Stability vs Runtime (crossed boxplots)');

% Draw boxplots (vertical) for every item
for k = 1:numel(stats)
    s = stats(k);
    draw_vertical_box(ax, s.stabVals, s.runMean, colors(s.methodIdx,:), boxWidth);
end
% Ensure whiskers solid for readability (some MATLAB versions tag whiskers)
try
    w = findobj(ax,'Tag','Whisker'); set(w,'LineStyle','-');
catch
    % ignore if tags differ by version
end

% Draw runtime markers (dot + horizontal bar)
for k = 1:numel(stats)
    s = stats(k);
    draw_runtime_marker(ax, s.runMean, s.stabMed, s.runStd, colors(s.methodIdx,:), markers{s.methodIdx});
end

% -------------- Place labels (now that all obstacles exist) --------------
% Place labels in order of *descending box height* (more constrained first)
[~, order] = sort(arrayfun(@(s) s.stabMax - s.stabMin, stats), 'descend');
placedRects = obstacles; % start with box & bar obstacles (so labels avoid them)
labelHandles = gobjects(numel(stats),1);
for idx = order
    s = stats(idx);
    labelStr = sprintf('%s\n%s', methodNames{s.methodIdx}, datasetNames{s.datasetIdx});
    [lx, ly, lrect] = place_label_near_box(ax, labelStr, s.bpRect, placedRects, labelFontSize);
    % draw label
    labelHandles(idx) = text(ax, lx, ly, labelStr, ...
        'FontSize', labelFontSize, 'Color', colors(s.methodIdx,:), ...
        'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
        'BackgroundColor','none');
    % add label rect to placedRects so next labels avoid it
    placedRects = [placedRects; lrect];
end

% -------------- Legend (dummy handles) --------------
hold(ax,'on');
hleg = gobjects(1,2);
for m = 1:2
    hleg(m) = plot(ax, NaN, NaN, markers{m}, 'MarkerFaceColor', colors(m,:), ...
        'MarkerEdgeColor','k', 'LineStyle','none', 'DisplayName', methodNames{m});
end
%legend(ax, hleg, methodNames, 'Location','best');

% -------------- xticks as 10^p (LaTeX) --------------
apply_log_xticks(ax, allRunLog);

hold(ax,'off');
end

%% -------------------- Helper functions -------------------- %%

function stats = collect_stats(stabCell, timeCell)
% Returns struct array with fields:
% datasetIdx, methodIdx, stabVals, runVals, runLogVals, runMean, runStd, stabMin, stabMed, stabMax
stats = struct([]);
[nR, nC] = size(stabCell);
k = 0;
for i = 1:nR
    for j = 1:nC
        svals = stabCell{i,j};
        rvals = timeCell{i,j};
        if isempty(svals) || isempty(rvals), continue; end
        k = k + 1;
        stats(k).datasetIdx = i;
        stats(k).methodIdx  = j;
        stats(k).stabVals   = svals(:);
        stats(k).runVals    = rvals(:);
        stats(k).runLogVals = log10(rvals(:));
        stats(k).runMean    = mean(stats(k).runLogVals);
        stats(k).runStd     = std(stats(k).runLogVals);
        stats(k).stabMin    = min(svals(:));
        stats(k).stabMed    = median(svals(:));
        stats(k).stabMax    = max(svals(:));
    end
end
end

function xr = compute_xlim(allRunLog)
mn = min(allRunLog); mx = max(allRunLog);
if mn == mx
    xr = [mn-0.5, mx+0.5];
else
    pad = 0.12 * (mx - mn);
    xr = [mn - pad, mx + pad];
end
end

function yr = compute_ylim(allStab)
mn = min(allStab); mx = max(allStab);
if mn == mx
    yr = [mn-0.5, mx+0.5];
else
    pad = 0.1 * (mx - mn);
    yr = [mn - pad, mx + pad];
end
end

function draw_vertical_box(ax, stabVals, xPos, color, widthVal)
% wrapper around MATLAB's boxplot to draw a vertical box at xPos
holdState = ishold(ax); hold(ax,'on');

% boxplot(ax, stabVals, 'positions', xPos, 'Widths', widthVal, ...
%     'Colors', color, 'Symbol', '', 'Whisker', 1.5);
% % ensure thin lines (optional)
% hp = findobj(ax,'Tag','Box');
% set(hp,'LineWidth',0.9);

% Draw vertical boxchart
bc = boxchart(ax, repmat(xPos, size(stabVals)), stabVals, ...
    'BoxWidth', widthVal, 'MarkerStyle','none');
bc.BoxFaceColor = color;
bc.LineWidth = 1;  % make boxes and whiskers solid

if ~holdState, hold(ax,'off'); end
end

function draw_runtime_marker(ax, xMean, yMed, xStd, color, marker)
% draw horizontal bar ±std at yMed and a center marker at (xMean,yMed)
line(ax, [xMean - xStd, xMean + xStd], [yMed yMed], 'Color', color, 'LineWidth', 1.5);
plot(ax, xMean, yMed, marker, 'MarkerFaceColor', color, 'MarkerEdgeColor', 'k', 'MarkerSize', 7);
end

function apply_log_xticks(ax, allRunLog)
minPow = floor(min(allRunLog));
maxPow = ceil(max(allRunLog));
xt = minPow:maxPow;
xticks(ax, xt);
labs = arrayfun(@(p) sprintf('$10^{%d}$', p), xt, 'UniformOutput', false);
xticklabels(ax, labs);
set(ax, 'TickLabelInterpreter', 'latex');
end

%% Label placement that is aware of all obstacles
function [cx, cy, rect] = place_label_near_box(ax, labelStr, bpRect, obstacles, fontSize)
% Try anchored positions around bpRect. Use label Extent in data units.
% Returns chosen center (cx,cy) and rect [x y w h].

% measure label size (Extent) in data units
tmp = text(ax, 0, 0, labelStr, 'Units', 'data', 'Visible', 'off', 'FontSize', fontSize);
ext = get(tmp, 'Extent'); delete(tmp);
lw = ext(3); lh = ext(4);

% box center
bx = bpRect(1); by = bpRect(2); bw = bpRect(3); bh = bpRect(4);
boxCenter = [bx + bw/2, by + bh/2];

% axis extents
axExt = axis(ax); % [xmin xmax ymin ymax]
baseMarginX = 0.01 * (axExt(2)-axExt(1));
baseMarginY = 0.01 * (axExt(4)-axExt(3));

% Candidate centers anchored just outside the box edges/corners (progressively increasing margin)
marginScales = [1, 1.25, 1.6, 2.2, 3];
best = struct('score', -Inf, 'cx', [], 'cy', [], 'rect', []);
for ms = marginScales
    mx = baseMarginX * ms;
    my = baseMarginY * ms;
    cands = candidate_centers_from_box(bpRect, [lw, lh], [mx, my]);
    for r = 1:size(cands,1)
        cxTry = cands(r,1); cyTry = cands(r,2);
        rectTry = [cxTry - lw/2, cyTry - lh/2, lw, lh];
        % ensure inside axes
        if rectTry(1) < axExt(1) || rectTry(1)+rectTry(3) > axExt(2) || ...
                rectTry(2) < axExt(3) || rectTry(2)+rectTry(4) > axExt(4)
            continue;
        end
        % check overlap with obstacles
        [overlaps, minClear] = clearance_to_obstacles(rectTry, obstacles);
        if overlaps
            continue; % reject overlapping candidates at this margin
        end
        % score = clearance - 0.08*distance-to-box-center (prefer closer)
        dist = hypot(cxTry - boxCenter(1), cyTry - boxCenter(2));
        score = minClear - 0.08 * dist;
        if score > best.score
            best.score = score;
            best.cx = cxTry; best.cy = cyTry; best.rect = rectTry;
        end
    end
    if ~isempty(best.cx), break; end
end

% fallback: if no non-overlapping candidate found, choose candidate with minimal overlap
if isempty(best.cx)
    % iterate all candidates across larger margins and pick minimal-overlap one
    allCands = [];
    for ms = marginScales
        mx = baseMarginX * ms; my = baseMarginY * ms;
        allCands = [allCands; candidate_centers_from_box(bpRect, [lw lh], [mx my])]; %#ok<AGROW>
    end
    bestOverlap = inf;
    for r = 1:size(allCands,1)
        cxTry = allCands(r,1); cyTry = allCands(r,2);
        rectTry = [cxTry - lw/2, cyTry - lh/2, lw, lh];
        % clamp to axes
        rectTry = clamp_rect_to_axes(rectTry, axExt);
        [overlaps, minClear] = clearance_to_obstacles(rectTry, obstacles);
        if overlaps
            % compute total overlapping area as penalty (approx)
            overlapArea = 0;
            for o=1:size(obstacles,1)
                if rect_intersects(rectTry, obstacles(o,:))
                    overlapArea = overlapArea + rect_intersection_area(rectTry, obstacles(o,:));
                end
            end
        else
            overlapArea = 0;
        end
        if overlapArea < bestOverlap
            bestOverlap = overlapArea;
            best.cx = cxTry; best.cy = cyTry; best.rect = rectTry;
        end
    end
end

% ensure we always return something
if isempty(best.cx)
    % place north of box with large margin
    [best.cx, best.cy] = north_of_box(bpRect, [lw, lh], [3*baseMarginX, 3*baseMarginY]);
    best.rect = [best.cx - lw/2, best.cy - lh/2, lw, lh];
    best.rect = clamp_rect_to_axes(best.rect, axis(ax));
end

cx = best.cx; cy = best.cy; rect = best.rect;
end

%% ---------- small geometric helpers ----------
function cands = candidate_centers_from_box(bpRect, labelSize, margins)
% bpRect = [x y w h], labelSize = [lw lh], margins = [mx my]
x = bpRect(1); y = bpRect(2); w = bpRect(3); h = bpRect(4);
lw = labelSize(1); lh = labelSize(2);
mx = margins(1); my = margins(2);

xc = x + w/2; yc = y + h/2;
xl = x - lw/2 - mx;            xr = x + w + lw/2 + mx;
yb = y - lh/2 - my;            yt = y + h + lh/2 + my;

cands = [
    xc, yt;  % N
    xc, yb;  % S
    xr, yc;  % E
    xl, yc;  % W
    xr, yt;  % NE
    xl, yt;  % NW
    xr, yb;  % SE
    xl, yb   % SW
    ];
end

function [overlaps, minClear] = clearance_to_obstacles(L, obstacles)
% L = [x y w h]; obstacles = Nx4
if isempty(obstacles)
    overlaps = false; minClear = inf; return;
end
overlaps = false; dmin = inf;
for i = 1:size(obstacles,1)
    O = obstacles(i,:);
    if rect_intersects(L, O)
        overlaps = true; minClear = 0; return;
    else
        d = rect_distance(L, O);
        if d < dmin, dmin = d; end
    end
end
minClear = dmin;
end

function yes = rect_intersects(A, B)
yes = ~(A(1)+A(3) <= B(1) || B(1)+B(3) <= A(1) || ...
    A(2)+A(4) <= B(2) || B(2)+B(4) <= A(2));
end

function d = rect_distance(A, B)
dx = max([B(1) - (A(1)+A(3)), A(1) - (B(1)+B(3)), 0]);
dy = max([B(2) - (A(2)+A(4)), A(2) - (B(2)+B(4)), 0]);
d = hypot(dx, dy);
end

function R = clamp_rect_to_axes(R, ext)
% R = [x y w h], ext = [xmin xmax ymin ymax] from axis()
if R(1) < ext(1), R(1) = ext(1); end
if R(2) < ext(3), R(2) = ext(3); end
if R(1)+R(3) > ext(2), R(1) = ext(2) - R(3); end
if R(2)+R(4) > ext(4), R(2) = ext(4) - R(4); end
end

function [nx, ny] = north_of_box(bpRect, labelSize, margins)
x = bpRect(1); y = bpRect(2); w = bpRect(3); h = bpRect(4);
lw = labelSize(1); lh = labelSize(2);
mx = margins(1); my = margins(2);
nx = x + w/2;
ny = y + h + lh/2 + my;
end

function A = rect_intersection_area(R1, R2)
x1 = max(R1(1), R2(1));
x2 = min(R1(1)+R1(3), R2(1)+R2(3));
y1 = max(R1(2), R2(2));
y2 = min(R1(2)+R1(4), R2(2)+R2(4));
if x2>x1 && y2>y1
    A = (x2-x1)*(y2-y1);
else
    A = 0;
end
end
