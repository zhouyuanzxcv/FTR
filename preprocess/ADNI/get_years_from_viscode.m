function [years] = get_years_from_viscode(data)
%GET_YEARS_FROM_VISCODE Summary of this function goes here
%   Detailed explanation goes here

VISCODE = data.VISCODE;

if ismember("VISCODE2", data.Properties.VariableNames)
    VISCODE2 = data.VISCODE2;
else
    VISCODE2 = VISCODE;
end

% Initialize an array to store the converted data
years = zeros(size(VISCODE));

% Loop through the data and apply conversions
for i = 1:numel(VISCODE)
    value = VISCODE{i};
    value2 = VISCODE2{i};

    if startsWith(value, 'm')
        years(i) = str2double(extractAfter(value, 'm'))/12;
    elseif startsWith(value2, 'm')
        years(i) = str2double(extractAfter(value2, 'm'))/12;
    else
        years(i) = nan;
    end
end

end

