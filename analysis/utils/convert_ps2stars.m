function stars = convert_ps2stars(ps)
%CONVERT_PS2STARS Summary of this function goes here
%   Detailed explanation goes here

stars = cell(size(ps));

for i = 1:length(ps)
    stars{i} = p2star(ps(i));
end


end

function s = p2star(p)
if p < 0.001
    s = '***';
elseif p < 0.01
    s = '**';
elseif p < 0.05
    s = '*';
else
    s = '';
end

end