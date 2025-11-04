function [s,f] = ftr_base(X, rec_pts, sigma, xmin, xmax, options)
% Input:
%   sigma - standard deviation of the noise

if nargin < 6
    options = [];
end

nbiom = size(X,2);

if numel(sigma) == 1
    sigma = repmat(sigma, nbiom, 1);
end

if numel(xmin) == 1
    xmin = repmat(xmin, 1, nbiom);
end

params = [];

params.options = options;
params.X = X;
params.rec_pts = rec_pts;
params.sigma = sigma;
params.xmin = xmin;
params.xmax = xmax;

% if all(sigma == 0)
if 1 % use fixed lambda
    f = ftr_base_impl(params, []);
else
    lambdas = [0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000];
    start_lambda = parse_param(options, 'lambda_start', 10);
    eval_fun = @eval_obj_fun;
    [lambda, extra_best, obj_vals, best_ind] = search_for_lambda_with_step(...
        params, lambdas, start_lambda, eval_fun);
    f = extra_best.f;
end


% s is the same for all the biomarkers
s = linspace(0,1,rec_pts);

end


%% The hyperparameters are used here

function f = ftr_base_impl(params, lambda)
options = params.options;
X = params.X;
rec_pts = params.rec_pts;
sigma = params.sigma;
xmin = params.xmin;
xmax = params.xmax;

nbiom = size(X,2);

% solve the problem 
method = options.deconv;

f = zeros(nbiom, rec_pts);

if ~isempty(lambda)
    options.lambda = lambda;
end

for j = 1:nbiom
    switch method
        case 'search_max'
            [s,f(j,:)] = ftr_search_max(X(:,j), rec_pts, sigma(j), options);
        case 'min_max'
            [s,f(j,:)] = ftr_min_max(X(:,j), rec_pts, sigma(j), xmin(j), xmax(j), options);
        case 'approx'
            [s,f(j,:)] = ftr_approx(X(:,j), rec_pts, sigma(j), xmin(j), xmax(j));
    end
end

end

function [obj_val, extra] = eval_obj_fun(params, lambda)

X = params.X;

sigma = params.sigma;

% solve the problem 
f = ftr_base_impl(params, lambda);

extra = [];
extra.f = f;

% evaluate the solution
PTID = (1:size(X,1))';
proption = 1;
[loglik, log_pdf, log_int, log_int_grouped] = cal_loglik(X, PTID, f, proption,sigma);
obj_val = -loglik;

end
