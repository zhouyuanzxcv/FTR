function [outputArg1,outputArg2] = run_OASIS_pipeline(inputArg1,inputArg2)
%RUN_OASIS_PIPELINE Summary of this function goes here
%   Detailed explanation goes here

[status, pythonPath] = system('which python3');
pythonPath = strtrim(pythonPath);  % remove newline

if status ~= 0 || isempty(pythonPath)
    error('Python not found on PATH. Please install or add to PATH.');
end

% install required python packages
cmd = sprintf('"%s" -m pip install --user pandas', pythonPath);
system(cmd);

cmd = sprintf('"%s" -m pip install --user scipy', pythonPath);
system(cmd);

% run the preprocess python code
cmd = sprintf('"%s" preprocess_oasis.py', pythonPath);
system(cmd);

preprocess_data;

copy_files_to_destination_OASIS();

end

