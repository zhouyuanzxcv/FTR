function data = add_baseline_info(data)
%GET_YEARS_FROM_EXAMDATE Summary of this function goes here
%   Detailed explanation goes here

data_adni_merge = readtable('./data/nacc_clinical.csv',VariableNamingRule='preserve');


clinical_sorted = sortrows(data_adni_merge, {'NACCID','NACCVNUM'});
[~, ia] = unique(clinical_sorted.NACCID);
bl_dates = clinical_sorted(ia, ["NACCID","NACCVNUM","VISITYR","VISITMO","VISITDAY",...
    "NACCUDSD","BIRTHYR","BIRTHMO","SEX","EDUC","NACCNE4S"]);

% Verify that all the NACCVNUMs are 1
assert(all(bl_dates.NACCVNUM == 1));

num_data = size(data,1);
converted_dates = NaN(num_data, 1);
diagnosis_bl = NaN(num_data, 1);
age_bl = NaN(num_data, 1);
gender = NaN(num_data, 1);
educ = NaN(num_data, 1);
ne4s = NaN(num_data, 1);

for i = 1:size(data, 1)
    ind = strcmp(data{i, "naccid"}, bl_dates{:, "NACCID"});
    bl_date_i = bl_dates{ind, {'VISITYR', 'VISITMO', 'VISITDAY'}};
    bl_dx = bl_dates{ind, "NACCUDSD"};
    bl_birth = bl_dates{ind, {'BIRTHYR','BIRTHMO'}};
    bl_gender = bl_dates{ind, 'SEX'};
    bl_educ = bl_dates{ind, 'EDUC'};
    bl_ne4s = bl_dates{ind, 'NACCNE4S'};

    if any(ind) % there is a record in the demographic table
        
        educ(i) = bl_educ;
        gender(i) = double(bl_gender == 1);
        diagnosis_bl(i) = bl_dx;
        ne4s(i) = bl_ne4s;
    
        % calculate years from baseline
        bl_date_i = datetime(bl_date_i(1,:));
    
        curr_date_i = data{i, "scandt"};
        years_diff = yearfrac(bl_date_i, curr_date_i, 0);
        
        converted_dates(i) = years_diff;

        % calculate age at baseline
        birthday = datetime([bl_birth, 15]);
        age_bl(i) = yearfrac(birthday, bl_date_i, 0);
    end
end

years = converted_dates;

data = addvars(data, years, 'After', 'naccid');

diagnosis_bl = convert_NACCUDSD_to_dx(diagnosis_bl);
data = addvars(data, diagnosis_bl);

AGE_baseline = age_bl;
data = addvars(data, AGE_baseline);

PTGENDER = gender;
data = addvars(data, PTGENDER);

PTEDUCAT = educ;
data = addvars(data, PTEDUCAT);

APOE4 = ne4s;
APOE4(APOE4==9) = NaN;
data = addvars(data, APOE4);

end

