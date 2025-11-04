function combined_data = combine_regions(data, map, biomarker_ver, modality)
%COMBINE_REGIONS Summary of this function goes here
%   Detailed explanation goes here
% Create a new table to store the combined data
biomarkers = unique(map.(biomarker_ver));

combined_data = table();

% Iterate through the unique biomarkers
for i = 1:numel(biomarkers)
    comb_biomarker = biomarkers{i};
    
    % Find the rows in map where 'BIOMARKER' matches the current biomarker
%     matching_rows1 = strcmp(map.BIOMARKER, biomarker);
    
%     if any(matching_rows1)
        % Get the corresponding 'FLDNAME' 
%         biomerker = map.BIOMARKER(matching_rows1);
%         comb_biomerker = map.(biomarker_ver)(matching_rows1);
        
    % Identify the corresponding columns in data
    matching_rows2 = strcmp(map.(biomarker_ver), comb_biomarker);
    columns_to_combine = map.FLDNAME(matching_rows2);

    % Sum the columns in data and store the result in the new column
    column_inds = ismember(columns_to_combine, data.Properties.VariableNames);
    valid_columns = any(column_inds);

    if valid_columns
        columns_to_combine = columns_to_combine(column_inds);
        if strcmp(modality, 'mri')
            % Sum up the specified columns and store the result in combined_data
            combined_data.(comb_biomarker) = sum(data{:, columns_to_combine}, 2);
        elseif strcmp(modality, 'pet')
            % \TODO: For column Abeta-positive, columns_weights are
            % directly the original SUMMARYSUVR_WHOLECEREBNORM_1.11CUTOFF.
            % Hence, if the original value is 0, weights become NaN and
            % the entries become NaN
            columns_weights = strrep(columns_to_combine, '_SUVR', '_VOLUME');
            weights = data{:, columns_weights};
            values = data{:, columns_to_combine};
            weights = weights ./ repmat(sum(weights, 2), 1, size(weights,2));
            combined_data.(comb_biomarker) = sum(values .* weights, 2);
        end
    end
   
end

end

