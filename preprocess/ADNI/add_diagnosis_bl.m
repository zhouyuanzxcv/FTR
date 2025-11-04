function data = add_diagnosis_bl(data)
%ADD_DIAGNOSIS_BL Summary of this function goes here
%   Detailed explanation goes here
[years1, diagnosis_bl] = get_years_from_examdate(data); 

data = addvars(data,diagnosis_bl,'After','years');

end

