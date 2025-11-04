function [train_data,test_data,all_data,biomarker_name,options] = load_dataset(dataset_name, train_inds, test_inds, options)
%LOAD_DATASET Summary of this function goes here
%   Detailed explanation goes here

%% load csv file

input_file_name = options.input_file_name;

if ~strcmp(dataset_name, 'custom')
    input_data = readtable(['./input/',input_file_name],'VariableNamingRule','preserve');
else
    input_data = readmatrix(['./input/',input_file_name]);
    [nsamp, nbioms] = size(input_data);
    RID = (1:nsamp)';
    years = zeros(nsamp, 1);
    diagnosis = ones(nsamp, 1);
    biomarker_names = cellfun(@(x) ['Biomarker',num2str(x)], num2cell(1:nbioms), ...
        'UniformOutput', false);
    input_data1 = num2cell(input_data, 1);
    input_data = table(RID, years, diagnosis, input_data1{:}, 'VariableNames', ...
        [{'RID','years','diagnosis'},biomarker_names]);
end

input_data1 = sortrows(input_data, {'RID','years'});

if ~all(all(input_data{:,{'RID','years'}} == input_data1{:,{'RID','years'}}))
    warning('Input data is not sorted according to RID and years. Sort it in load_dataset.');
end
input_data = input_data1;

PTID= input_data.RID;

nsamp = size(input_data,1);

labels = input_data.diagnosis;

%% select the case group (rows)

rem_inds = ones(nsamp,1);

fprintf('Original number of individuals: %d. \n', length(unique(PTID)));

if options.group_sel % not false
    if options.group_sel == 2
%         error('Setting group_sel = 2 to select both the case and control groups has been deprecated');
        disp('Select the validation group for analysis')

        idx_validation = input_data.group == 2;
        idx_control = input_data.group == 0;

        show_data_statistics(input_data, idx_validation, idx_control);

        rem_inds1 = idx_validation;

        rem_inds = rem_inds & rem_inds1;

        fprintf('After group selection, %d individuals left. \n', length(unique(PTID(rem_inds))));
    elseif options.group_sel == 10
        disp('Select both the case group and control group for analysis');

        idx_case = input_data.group == 1;
        idx_control = input_data.group == 0;

        show_data_statistics(input_data, idx_case, idx_control);

        rem_inds1 = idx_case | idx_control;

        rem_inds = rem_inds & rem_inds1;

        fprintf('After group selection, %d individuals left. \n', length(unique(PTID(rem_inds))));

    else
        disp('Select only the case group for analysis');

        idx_case = input_data.group == 1;
        idx_control = input_data.group == 0;

        show_data_statistics(input_data, idx_case, idx_control);

        rem_inds1 = idx_case;

        rem_inds = rem_inds & rem_inds1;

        fprintf('After group selection, %d individuals left. \n', length(unique(PTID(rem_inds))));
    end

end

if options.diagnosis_sel
    % Keep time series that have at least 1 AD time point
    disp('Further select subjects based on diagnosis')
    rem_inds1 = get_diagnosis_sel_indicators(PTID,labels);

    rem_inds = rem_inds & rem_inds1;

    show_data_statistics(input_data, rem_inds, idx_control);

    fprintf('After diagnosis selection, %d individuals left\n', ...
        length(unique(PTID(rem_inds))));

end

%% select the biomarkers (columns)

biomarker_type = parse_param(options, 'biomarker_type', []);

% default biomarkers, i.e. regional brain volumes
if isempty(biomarker_type)
    bcr = parse_param(options, 'biomarker_column_range', [4,0]);
    if bcr(2) <= 0
        biomarker_name = input_data.Properties.VariableNames(bcr(1):end+bcr(2));
        bioms = input_data{:,bcr(1):end+bcr(2)};
    else
        biomarker_name = input_data.Properties.VariableNames(bcr(1):bcr(2));
        bioms = input_data{:,bcr(1):bcr(2)};
    end
elseif strcmp(biomarker_type, 'ctv_hv')
    biomarker_name = {'CTV','HV'};
    bioms = input_data{:, biomarker_name};
end

% check if there is nan in the data
if any(any(isnan(bioms)))
    error('There is NaN in the data in load_dataset');
end


%% split the data

data_split = parse_param(options, 'data_split', '');

% split into training and testing
if isempty(train_inds) || isempty(test_inds)
    rem_list = find(rem_inds == 1);
    switch data_split
        case 'cross_patients'
            PTID_uniq = unique(PTID(rem_list));
            N = length(PTID_uniq);
            [trainInd,testInd] = dividerand(N,0.8,0.2);
            train_inds = search_id(PTID_uniq(trainInd),PTID);
            test_inds = search_id(PTID_uniq(testInd),PTID);
        case 'cross_points'
            N = sum(rem_inds);
            [trainInd,testInd] = dividerand(N,0.8,0.2);
            train_inds = arr2vec(rem_list(trainInd),nsamp);
            test_inds = arr2vec(rem_list(testInd),nsamp);           
        case 'last_point'
            PTID_uniq = unique(PTID(rem_list));
            N = length(PTID_uniq);
            test_inds = zeros(nsamp,1);
            for i = 1:N
                inds = find(PTID == PTID_uniq(i) & rem_inds == 1);
                % If a subject has only 1 point, this point is assigned
                % to the training set
                if length(inds) > 1
                    test_inds(inds(end)) = 1;
                end
            end
            train_inds = rem_inds - test_inds;
        case 'baseline'
            PTID_uniq = unique(PTID(rem_list));
            N = length(PTID_uniq);
            train_inds = zeros(nsamp,1);
            % If a subject has only 1 point, this point is assigned to the
            % training set
            for i = 1:N
                inds = find(PTID == PTID_uniq(i) & rem_inds == 1);
                train_inds(inds(1)) = 1;
            end
            test_inds = rem_inds - train_inds;
        case 'first_AD'
            PTID_uniq = unique(PTID(rem_list));
            N = length(PTID_uniq);
            train_inds = zeros(nsamp,1);
            for i = 1:N
                inds = find(PTID == PTID_uniq(i) & rem_inds == 1);
                % find the point when the patient is diagnosed as AD for
                % the first time
                k = find(labels(inds) == 1, 1);
                if ~isempty(k)
                    train_inds(inds(k)) = 1;
                end
            end
            test_inds = rem_inds - train_inds;
        otherwise
            train_inds = rem_inds;
            test_inds = zeros(nsamp,1);
    end
end

options.train_inds = train_inds;
options.test_inds = test_inds;

train_data = read_inds(input_data,bioms,train_inds);
test_data = read_inds(input_data,bioms,test_inds);
all_data = read_inds(input_data,bioms,ones(nsamp,1));

end

function data = read_inds(data_all,bioms,inds)
data = [];
data.RID = data_all.RID(inds==1);
data.years = data_all.years(inds==1);
data.labels = data_all.diagnosis(inds==1);
data.vols = bioms(inds==1,:);
end

function vec = arr2vec(arr,len)
    vec = zeros(len,1);
    for i = arr
        vec(i) = 1;
    end
end

function vec = search_id(PTID_sel,PTID)
    vec = zeros(length(PTID),1);
    for i = 1:length(PTID_sel)
        vec(PTID == PTID_sel(i)) = 1;
    end
end

function show_data_statistics(input_data, idx_case, idx_control)
PTID= input_data.RID;

fprintf('Case points/subjects: %d/%d; control points/subjects: %d/%d \n', ...
        length(find(idx_case)), length(unique(PTID(idx_case))), ...
        length(find(idx_control)), length(unique(PTID(idx_control))));

ages_AD = input_data.AGE(idx_case);
ages_CN = input_data.AGE(idx_control);

if ismember('PTGENDER', input_data.Properties.VariableNames)
    genders = input_data.PTGENDER;
elseif ismember('gender', input_data.Properties.VariableNames)
    genders = input_data.gender;
end

genders_AD = genders(idx_case);
genders_CN = genders(idx_control);

% each row has 5 numbers: min, Q1, median, Q3, max. First row: number of time points per
% subject. Second row: duration (year) per subject. Third row: interval in year
% between adjacent time points.
stats_control = get_group_statistics(input_data(idx_control,:));
stats_case = get_group_statistics(input_data(idx_case,:));

disp(['Case (row 1)/control (row 2) age mean, age std, gender(male%),',...
    ' (Q1, median, Q3) for Ts, duration, interval']);
[mean(ages_AD), std(ages_AD), length(find(genders_AD == 1))/length(genders_AD), ...
    stats_case(1,2:4), stats_case(2,2:4), stats_case(3,2:4); ...
    mean(ages_CN), std(ages_CN), length(find(genders_CN == 1))/length(genders_CN), ...
    stats_control(1,2:4), stats_control(2,2:4), stats_control(3,2:4)]
end

function stats = get_group_statistics(data)
[C,ia,ic] = unique(data.RID);

Ts = accumarray(ic, ones(size(data.RID)));
% num_pts_stats = [mean(Ts), std(Ts)];
num_pts_stats = [min(Ts), prctile(Ts, 25), median(Ts), prctile(Ts, 75), max(Ts)];

duration = data.years([ia(2:end)-1;end]) - data.years([ia(1:end)]);
% duration_stats = [mean(duration), std(duration)];
duration_stats = [min(duration), prctile(duration, 25), median(duration), ...
    prctile(duration, 75), max(duration)];

interval = data.years(2:end) - data.years(1:end-1);
interval(ia(2:end)-1) = [];
% interval_stats = [mean(interval), std(interval)];
interval_stats = [min(interval), prctile(interval, 25), median(interval), ...
    prctile(interval, 75), max(interval)];

stats = [num_pts_stats; duration_stats; interval_stats];
end