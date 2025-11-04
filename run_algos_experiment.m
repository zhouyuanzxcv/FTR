function [outputArg1,outputArg2] = run_algos_experiment(postfix)
%RUN_ALGOS_EXPERIMENT Summary of this function goes here
%   Detailed explanation goes here
data_names = {'ADNI_FSX_HS','ADNI_FSX_HM','ADNI_FSX_LS','ADNI_FSX_LM',...
    'ADNI_FSX_HS_DS','ADNI_FSX_HM_DS','ADNI_FSX_LS_DS','ADNI_FSX_LM_DS',...
    'OASIS3_ROD1_HS','OASIS3_ROD1_HM','OASIS3_ROD1_LS','OASIS3_ROD1_LM',...
    'NACC_HS','NACC_HM','NACC_LS','NACC_LM',...
    };

data_names = {'ADNI_FSX_HS','ADNI_FSX_HS_DS','OASIS3_ROD1_HS','NACC_HS'};
% method_names = {'FTR_MCEM','sustain'};
method_names = {'FTR_MCEM'};
% method_names = {'sustain'};

re_subtype_staging = 1;

skip_combinations = { ...
%     'ADNI_FSX_LM','sustain'; ...
%     'ADNI_FSX_LS','sustain'; ...
    'ADNI_FSX_HS','sustain'; ...
    'ADNI_FSX_HM','sustain'; ...
    'ADNI_FSX_HS_DS','sustain'; ...
    'ADNI_FSX_HM_DS','sustain'; ...
    'OASIS3_ROD1_HS','sustain'; ...
    'OASIS3_ROD1_HM','sustain'; ...
    'NACC_HS','sustain'; ...
    'NACC_HM','sustain' ...
    };

% nsubtypes = [1,3];
nsubtypes = [3];

if nargin < 1
    postfix = '';
end

for i = 1:length(data_names)
    for j = 1:length(method_names)
        for k = nsubtypes
            skip = false;
            for l = 1:size(skip_combinations,1)
                if strcmp(data_names{i}, skip_combinations{l,1}) && ...
                        strcmp(method_names{j}, skip_combinations{l,2}) 
                    skip = true;
                    break;
                end
            end

            if skip
                continue;
            end

            disp(['=========================== Run ', ...
                method_names{j},' on ',data_names{i},' with ',num2str(k), ...
                ' subtypes =========================== ']);

            options = [];
            options.re_subtype_staging = re_subtype_staging;
            options.nsubtype = k;
            options.postfix = postfix;
            run_algo(data_names{i}, method_names{j}, options);
        end
    end
end

end

