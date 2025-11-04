function diagnosis_bl = convert_NACCUDSD_to_dx(diagnosis_bl)
% treat normal cognition and cognitive impared (not MCI) both as CN
cn_idx = diagnosis_bl == 1 | diagnosis_bl == 2;
% MCI
mci_idx = diagnosis_bl == 3;
% AD
ad_idx = diagnosis_bl == 4;

diagnosis_bl(cn_idx) = 0;
diagnosis_bl(mci_idx) = 0.5;
diagnosis_bl(ad_idx) = 1;

end