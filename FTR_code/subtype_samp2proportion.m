function proption = subtype_samp2proportion(subtype, PTID, nsubtype)

if 0
    tmp = accumarray(PTID, subtype, [], @mean, [], 1);
    tmp1 = nonzeros(tmp);
    
    freq = tabulate(tmp1);
    proption = freq(:,3)/100;
else
    subtype = double(subtype == 1:nsubtype);
    
    [unique_PTIDs,ia,ic] = unique(PTID);
    [yy, xx] = ndgrid(ic, 1:nsubtype);
    subtype_grouped = accumarray([yy(:),xx(:)], subtype(:), [], @mean);
    proption = mean(subtype_grouped, 1);
end

end