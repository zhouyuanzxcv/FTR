function color = get_subtype_color(k)

traj_colors = cat(3, spring, summer, winter);

subtype_colors = squeeze(traj_colors(128,:,:))';

% subtype_colors = {[28 128 130] / 255; [229 87 9] / 255; ...
%     [80 29 138] / 255; [217 74 193] / 255};

color = subtype_colors(k,:);




end