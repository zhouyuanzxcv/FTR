function [s,f] = ftr_approx(X, rec_pts, sigma, xmin, xmax)
[F_y,F_x] = ecdf(X);
F_x(1) = F_x(1) - 1e-8;

if sigma == 0
    % take the inverse
    s = linspace(0, 1, rec_pts);
    f = interp1(F_y, F_x, s);
    f(f < xmin) = xmin;
    return;
end
% F_x = F_x(2:end);
% F_y = F_y(2:end);

[xs, s1] = deconvolution(F_x, F_y, rec_pts, sigma);
s = linspace(0, 1, rec_pts);
f = interp1(s1, xs, s);

f(f < xmin) = xmin;

stage_thresh = 0.98;
xmax = max(f(s <= stage_thresh));
f(f > xmax) = xmax;
end

function [xs, s1] = deconvolution(F_x, F_y, rec_pts, sigma)
% set the number of points for deconvolution 

cdf_x = linspace(min(F_x), max(F_x), rec_pts);
delta_x = cdf_x(2) - cdf_x(1);

cdf_y = interp1(F_x, F_y, cdf_x);

if 0
    nabla2_cdf = calc_nabla2(cdf_y, delta_x, 10);
elseif 0
    r = round(rec_pts/20);

    cdf_y1 = [zeros(1,r), cdf_y, ones(1,r)];
    inds_y = repmat((1:2*r+1)',1,rec_pts) + repmat((0:rec_pts-1),2*r+1,1);
    y_inds = reshape(cdf_y1(inds_y(:)), size(inds_y));
    x = [(1:2*r+1)' * delta_x, ones(2*r+1,1)];
    betas = inv(x'*x)*x'*y_inds;
    beta1 = betas(1,:);
else
    gauss_sigma = 5;
    pad_sz = 3*gauss_sigma;
    cdf_y1 = [zeros(1,pad_sz), cdf_y, ones(1,pad_sz)];
    cdf_y1 = imgaussfilt(cdf_y1, gauss_sigma, 'Padding', 'replicate');

    % cdf_y_left = cdf_y1([1,1:end-1]);
    % cdf_y_right = cdf_y1([2,2:end]);
    % nabla2_cdf = (cdf_y_left + cdf_y_right - 2*cdf_y1) / (delta_x^2);

    nabla2_cdf = calc_nabla2(cdf_y1, delta_x, 10);

    nabla2_cdf = nabla2_cdf(pad_sz+1:end-pad_sz);
end

cdf_y1 = cdf_y - sigma^2/2*nabla2_cdf;

% figure, plot(cdf_y);
% hold on; plot(cdf_y1);

xs = cdf_x;
s1 = cdf_y1;

% if xmax > xmin
    % xs = linspace(xmin, xmax, n);
% else
%     % if all the points are negative, i.e. max(X) <= 0, the trajectory must
%     % be flat
%     xs = [xmin, xmin + 1e-4];
%     s1 = [0, 1];
%     return;
% end

end

function nabla2_cdf = calc_nabla2(cdf_y, delta_x, h)
rec_pts = length(cdf_y);

left_inds = (1:rec_pts) - h;
inds_pos = (left_inds > 0);
cdf_y_left = [zeros(1, length(find(~inds_pos))), cdf_y(left_inds(inds_pos))];

right_inds = (1:rec_pts) + h;
inds_valid = right_inds <= rec_pts;
cdf_y_right = [cdf_y(right_inds(inds_valid)), ones(1, length(find(~inds_valid)))];

nabla2_cdf = (cdf_y_left + cdf_y_right - 2*cdf_y) / ((h*delta_x)^2);
end