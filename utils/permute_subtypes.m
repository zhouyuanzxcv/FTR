function [mdl, joindata] = permute_subtypes(mdl, joindata, P)
%PERMUTE_SUBTYPES Summary of this function goes here
%   Detailed explanation goes here

% permute the trajectories

nsubtype = size(P, 1);

traj_names = {'traj','re_traj'};
for i = 1:length(traj_names)
    if isfield(mdl, traj_names{i})
        traj_idx = mdl.(traj_names{i});
        [nbiom, num_int, nsubtype] = size(traj_idx);
        traj_i_perm = P * reshape(traj_idx, [num_int*nbiom, nsubtype])';
        traj_idx1 = reshape(traj_i_perm', [nbiom,num_int, size(traj_i_perm,1)]);
        mdl.(traj_names{i}) = traj_idx1;
    end
end

% permute the sigma and proportion
if isfield(mdl, 'sigma')
    mdl.sigma = mdl.sigma * P';
end

if isfield(mdl, 'proption')
    mdl.proption = P * mdl.proption;
end


% permute the subtype labels

subtype1 = joindata.subtype==1:nsubtype;
[~,ind] = max((subtype1 * P'), [], 2);
joindata.subtype = ind;

% permute the subtype probabilities

fn_name = @(x) ['subtype', num2str(x)];

joindata = permute_subtype_prob(joindata, nsubtype, fn_name, P);

fn_name = @(x) ['subtype', num2str(x), '_visit'];

joindata = permute_subtype_prob(joindata, nsubtype, fn_name, P);


end

function joindata = permute_subtype_prob(joindata, nsubtype, fn_name, P)
subtype_names = cellfun(fn_name, num2cell(1:nsubtype), 'UniformOutput', 0);

if all(ismember(subtype_names, joindata.Properties.VariableNames))
    permuted_subtype_prob = joindata{:,subtype_names} * P';
    if size(permuted_subtype_prob, 2) > nsubtype
        subtype_names1 = cellfun(fn_name, ...
            num2cell(1:size(permuted_subtype_prob, 2)), 'UniformOutput', 0);
    else
        subtype_names1 = subtype_names;
    end
    joindata{:,subtype_names1} = permuted_subtype_prob;
end
end