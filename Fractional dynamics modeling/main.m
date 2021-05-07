% define a common frequency (readings per second)
COMMON_FREQUENCY = 10;

% define the sample length
SAMPLE_LENGTH_IN_MINUTES = 10;

% how many extra random samples from the patient
EXTRA_RANDOM_SAMPLES = 200;

% specify the EDF file to be read
% all EDF files MUST follow this convention: if the patient is having the
% condition the file starts with 1_ to 4_ for COPD stages, otherwise if the patient is a control/
% healthy patient, then the file starts with 0_

FILE_PATH = './1_test_00.edf';
[~,fname,~] = fileparts(FILE_PATH);
fparts = strsplit(fname,'_');
y = str2num(fparts{1});
if y < 0 || y > 4
    disp('Wrong filename format!');
    return
end


% specify the signals we want to extract
signals = { 'Abdomen', 'Activity', 'Elevation', 'Flow' ,...
            'Nasal Pressure', 'Pleth', 'Pulse', 'Resp Rate', ...
            'RIP Sum', 'SpO2', 'SpO2 B-B', 'Thorax'};

% read the signals
[~, data] = inception_edfread(FILE_PATH, ...
    'targetSignals', signals, 'commonSampleRate', COMMON_FREQUENCY);

disp('Deleting outlier data...');
% define min-max for signals
outliers = [
    -1, 1;  NaN, 6;  NaN, NaN;  NaN, NaN; ...
    NaN, NaN;  NaN, NaN;  20, 230;   5, 40; ...
    -1, 1;  65, 100;   65, 100;  -1, 1];

% delete outlier data
indexes_to_delete = [];
for row=1:size(outliers, 1)
    min_cutoff = outliers(row, 1);
    max_cutoff = outliers(row, 2);
    
    if ~isnan(min_cutoff)
        indexes_to_delete = [indexes_to_delete find(data(row, :) < min_cutoff)];
    end
    
    if ~isnan(max_cutoff)
        indexes_to_delete = [indexes_to_delete find(data(row, :) > max_cutoff)];
    end
end

indexes_to_delete = unique(indexes_to_delete);
data(:, indexes_to_delete) = [];


% rescale between 0 and 1
disp('Rescaling data...');
rowmin = min(data, [], 2);
rowmax = max(data, [], 2);
L_BOUND = 0;
U_BOUND = 1;
data_scaled = rescale(data, L_BOUND, U_BOUND, 'InputMin',rowmin,'InputMax',rowmax);

disp('Collecting samples...');

sample_length = 60 * SAMPLE_LENGTH_IN_MINUTES * COMMON_FREQUENCY;

sample_points = int32(linspace(1, size(data_scaled, 2), floor(size(data_scaled, 2)/sample_length)+1));

% adjust for non-overlapping iteration
if sample_points(1) == 1
    sample_points(1) = 0;
end

diff = setdiff(0:size(data_scaled, 2)-sample_length, sample_points);
extra_random_samples = datasample(diff, EXTRA_RANDOM_SAMPLES, 'Replace', false);

% how many runs
runs = length(sample_points)-1;

% this is the matrix that holds the results
X = zeros(runs + EXTRA_RANDOM_SAMPLES, length(signals)^2+1);

fprintf('Running with %d samples in total, from which %d are spanning across all records and %d are random...\n',...
    runs + EXTRA_RANDOM_SAMPLES, runs, EXTRA_RANDOM_SAMPLES);

% running the 'span' samples
for i = 1:runs
    
    offset_start = sample_points(i)+1;
    offset_end = sample_points(i+1);
    
    fprintf('Processing span sample %d/%d with offset: %d - %d\n', i, runs, offset_start, offset_end);
    % subsample
    subsample = transpose(data_scaled(:,offset_start:offset_end));
    
    % run Gaurav's algorithm
    [Aout, B, order, u, relErr] = ...
        modelEst('sensInd', 1:length(signals), 'numInp', floor(length(signals)/2),...
            'data', subsample, 'silentFlag', 0);
        
    % make a row-vector out of the result and add it to the main dataset
    X(i,1:length(signals)^2) = reshape(Aout', 1, []);
end

% running the random samples
for i = 1:EXTRA_RANDOM_SAMPLES
    
    offset_start = extra_random_samples(i);
    offset_end = offset_start + sample_length;
    
    fprintf('Processing random sample %d/%d with offset: %d - %d\n', i, EXTRA_RANDOM_SAMPLES,...
        offset_start, offset_end);
    % subsample
    subsample = transpose(data_scaled(:,offset_start:offset_end));
    
    % run Gaurav's algorithm
    [Aout, B, order, u, relErr] = ...
        modelEst('sensInd', 1:length(signals), 'numInp', floor(length(signals)/2),...
            'data', subsample, 'silentFlag', 0);
        
    % make a row-vector out of the result and add it to the main dataset
    X(i + runs,1:length(signals)^2) = reshape(Aout', 1, []);
end

% add the y-label
X(:, length(signals)^2+1) = y;

% save the dataset as CSV
csvwrite(strcat(fname, '.csv'), X);

