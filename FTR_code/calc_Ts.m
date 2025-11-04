function Ts = calc_Ts(PTID)
%CALC_TS Summary of this function goes here
%   Detailed explanation goes here
unique_PTIDs = unique(PTID);
Ts = zeros(1, length(unique_PTIDs));
for i = 1:length(Ts)
    Ts(i) = length(find(unique_PTIDs(i) == PTID));
end

% [unique_RID, ia] = unique(PTID);
% Ts = [ia(2:end) - ia(1:end-1); length(PTID) - ia(end) + 1];

end

