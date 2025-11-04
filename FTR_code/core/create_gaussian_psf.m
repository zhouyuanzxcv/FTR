function [PSF_y, PSF_x] = create_gaussian_psf(delta_x, sigma)
%CREATE_GAUSSIAN_PSF Summary of this function goes here
%   Detailed explanation goes here

% Changed it from 2 sigma to 3 sigma. Sometimes, for 2 sigma, the left
% and the right part of trajectory gets unexpected jumps
size = ceil(3*sigma / delta_x);

PSF_x = (-size:size)*delta_x;
%PSF_y = exp(-(PSF_x-mu).^2 / (2*sigma^2));
PSF_y = exp(-PSF_x.^2 / (2*sigma^2));
PSF_y = PSF_y / sum(PSF_y);

end

