function [subtype_avg, sgm_avg, subtype_perms, sgm_perms] = average_subtypes(subtype_sel, sgm_sel)
subtype_template = subtype_sel(:,1);

nsubtype = max(subtype_sel(:));
template = double(subtype_template == 1:nsubtype);

subtype_perms = subtype_sel;
sgm_perms = sgm_sel;

subtype_summed = template;
sgm_summed = sgm_sel(:,:,1);
for k = 2:size(subtype_sel,2)
    subtype_k = double(subtype_sel(:,k) == 1:nsubtype);
    [P,subtype_k1] = permute_endmembers(template', subtype_k');
    sgm_k = P * sgm_sel(:,:,k)';
    
    [~, subtype_perm_k] = max(subtype_k1', [], 2);
    subtype_perms(:,k) = subtype_perm_k;
    sgm_perms(:,:,k) = sgm_k';
    
    subtype_summed = subtype_summed + subtype_k1';
    sgm_summed = sgm_summed + sgm_k';
end

subtype_averaged = subtype_summed / size(subtype_sel, 2);
[~, subtype_avg] = max(subtype_averaged, [], 2);

sgm_avg = sgm_summed / size(sgm_sel, 3);
end

