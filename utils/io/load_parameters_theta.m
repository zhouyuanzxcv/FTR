function mdl = load_parameters_theta(save_dir)
%LOAD_PARAMETERS_THETA Summary of this function goes here
%   Detailed explanation goes here

sigma = readmatrix(strcat(save_dir,'/pred_sigma.csv'));
nsubtype = size(sigma, 2);

for p = 1:nsubtype
    re_traj(:,:,p) = readmatrix([save_dir,'/trajectory',int2str(p),'.csv'])';
end

mdl.sigma = sigma;
mdl.re_traj = re_traj;

proption_path = strcat(save_dir,'/proportion.csv');

if exist(proption_path, 'file')
    proption = readmatrix(proption_path);
    mdl.proption = proption;
end

end

