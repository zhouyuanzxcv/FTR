function f1 = reparam_traj(f, N, options)
%REPARAM_TRAJ Reparameterize the trajectory according to the Euclidean
%distance.
% Input:
%   f - N1 x B matrix where each column is the trajectory of a biomarker.
%   N - number of points in the reparameterized trajectory
%
% Output:
%   f1 - N x B reparameterized trajectory

% x = x(:); y = y(:);

[N1,B] = size(f);  

if nargin < 2
    N = N1;
end

if nargin < 3
    options = [];
end

reparam_mode = parse_param(options, 'reparam_mode', 'L2');

% make it a circular list since we are dealing with closed contour
% x = [x;x(1)];
% y = [y;y(1)];

% dx = x([2:N+1])- x(1:N);
% dy = y([2:N+1])- y(1:N);
% d = sqrt(dx.*dx+dy.*dy);  % compute the distance from previous node for point 2:N+1

f0 = [zeros(1,B); f];

if strcmp(reparam_mode, 'L2')
    df = f0(2:end,:) - f0(1:end-1,:);
    d = sqrt(sum(df.^2, 2));
elseif strcmp(reparam_mode, 'L1')
    df = f0(2:end,:) - f0(1:end-1,:);
    d = sum(abs(df), 2);
end

d = [0;d];   % point 1 to the origin is 0

% now compute the arc length of all the points to point 1
% we use matrix multiply to achieve summing 
% M = length(d);
% d = (d'*uppertri(M,M))';

d = cumsum(d);

% now ready to reparametrize the closed curve in terms of arc length
maxd = d(end);

% if (maxd/RES<3)
%    error('RES too big compare to the length of original curve');
% end

di = linspace(0, maxd, N);

[~, idx, ~] = unique(d);
d = d(idx);
f0 = f0(idx,:);

f1 = zeros(N,B);
for i = 1:B
    f1(:,i) = interp1(d, f0(:,i), di);
end
    
% xi = interp1(d,x,di');
% yi = interp1(d,y,di');

% N = length(xi);

% if (maxd - di(length(di)) <RES/2)  % deal with end boundary condition
%    xi = xi(1:N-1);
%    yi = yi(1:N-1);
% end

end

