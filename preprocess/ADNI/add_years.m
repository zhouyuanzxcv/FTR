function data = add_years(data)

[years1, diagnosis_bl] = get_years_from_examdate(data); 

if 0
% check dates from examdate is the same as the viscode
years2 = get_years_from_viscode(data);

inds_conv = ~isnan(years1) & ~isnan(years2);
figure, histogram(abs(years1(inds_conv) - years2(inds_conv)));
end

% for fsx MRI, years1 has 203 nans, years2 has 3094 nans
years = years1;

data = addvars(data,years,'After','RID');

data(isnan(data.years), :) = [];

end
