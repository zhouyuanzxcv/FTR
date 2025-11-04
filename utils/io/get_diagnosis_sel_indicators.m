function rem_inds1 = get_diagnosis_sel_indicators(PTID,labels)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
nsamp = length(PTID);
PTID_abn = unique(PTID(labels==1));
rem_inds1 = zeros(nsamp,1);
for i = 1:length(PTID_abn)
    rem_inds1(PTID==PTID_abn(i)) = 1;
end

end

