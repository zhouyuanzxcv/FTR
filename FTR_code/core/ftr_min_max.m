function [s,f] = ftr_min_max(X, rec_pts, sigma, xmin, xmax, options)
if nargin < 6
    options = [];
end

[F_y,F_x] = ecdf(X);

if sigma == 0
    % take the inverse
    s = linspace(0, 1, rec_pts);
    f = interp1(F_y, F_x, s);
    f(f < xmin) = 0;
    return;
end
F_x = F_x(2:end);
F_y = F_y(2:end);

[xs, s1] = deconvolution_with_known_minmax(F_x, F_y, xmin, xmax, sigma, options);
s = linspace(0, 1, rec_pts);
f = interp1(s1, xs, s);

% f(f < xmin) = xmin;

stage_thresh = parse_param(options, 'stage_thresh', 0.97);
xmax = max(f(s <= stage_thresh));
f(f > xmax) = xmax;

end

function [xs, s1] = deconvolution_with_known_minmax(F_x, F_y, xmin, xmax, sigma, options)
% set the number of points for deconvolution 
n = 101;

if xmax > xmin
    xs = linspace(xmin, xmax, n);
else
    % if all the points are negative, i.e. max(X) <= 0, the trajectory must
    % be flat
    xs = [xmin, xmin + 1e-4];
    s1 = [0, 1];
    return;
end

lambda = parse_param(options, 'lambda', 20);
L1 = calc_Laplacian_1d(ones(1, n-1));
L1 = L1(2:end-1,:);

delta_x = xs(2) - xs(1);
[PSF_y, ~] = create_gaussian_psf(delta_x, sigma);
l = floor(length(PSF_y) / 2);

% construct Fy for deconvolution
Fx_left = delta_x * (-l+1:0) + xs(1);
Fx_center = xs(2:end-1);
Fx_right = delta_x * (0:l-1) + xs(end);
Fx = [Fx_left, Fx_center, Fx_right];
Fy = interp1(F_x, F_y, Fx, 'linear', NaN);
m = length(Fy);
% set the outside points to 0 or 1 based on whether the points are on the
% left side or the right side
inter_pt = find(Fy>0.5);
Fy(isnan(Fy) & (1:length(Fy)) < inter_pt(1)) = 0;
Fy(isnan(Fy) & (1:length(Fy)) > inter_pt(1)) = 1;

% create the kernel and b
K = zeros(m, m+2*l);
for i = 1:m
    K(i, i-1 + (1:2*l+1)) = PSF_y;
end
s1 = calc_s(K, L1, Fy, m, l, lambda);
xs = Fx(l+1:end-l);
s1 = [0, s1', 1];
xs = [Fx(l), xs, Fx(end-l+1)];
end

function [s, A, b, L, lambda] = calc_s(K, L1, Fy, n, l, lambda)
K2 = K(:,2*l+1:end-2*l);
K3 = K(:,end-2*l+1:end);
b = Fy' - K3*ones(2*l,1);
A = K2;

% create L
L = L1(1:n-2*l-2,1:n-2*l);

%% fmincon/quadprog optimization
m = n-2*l;
D = create_difference(m);
c = -0.0001*ones(m-1,1);
options = optimoptions('quadprog','Display','off');
s = quadprog(A'*A + lambda*(L'*L),-A'*b,D,c,[],[],0.0001*ones(m,1),0.9999*ones(m,1),[],options);
end

function D = create_difference(n)
D = zeros(n-1,n);
for i = 1:n-1
    D(i,i) = 1;
    D(i,i+1) = -1;
end
end