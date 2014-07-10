function s_cat_connectomes
% Load a series of conenctomes precomputed and fitted, eliminate the
% 0-weight fibers and concatenate them.
%
% Copyright (c) Franco Pestilli Stanford University 2014

% Get the base directory for the data
datapath = '/marcovaldo/frk/2t1/predator/';
subjects = {...
    'KK_96dirs_b2000_1p5iso', ...    'JW_96dirs_b2000_1p5iso', ...
    'HT_96dirs_b2000_1p5iso', ...
    'KW_96dirs_b2000_1p5iso', ...
    'FP_96dirs_b2000_1p5iso', ...
    'MP_96dirs_b2000_1p5iso', ...
    };

if notDefined('saveDir'), savedir = fullfile('/marcovaldo/frk/Dropbox','pestilli_etal_revision',mfilename);end

for isbj = 1:length(subjects)
    % File to load
    connectomesPath   = fullfile(datapath,subjects{isbj},'connectomes');
    feFileToLoad      = dir(fullfile(connectomesPath,sprintf('*recomputed*.mat')));
    keyboard
    
    fname = feFileToLoad(probIndex).name(1:end-4);
    feFileToLoad = fullfile(connectomesPath,fname);
    fprintf('[%s] Loading: \n%s\n ======================================== \n\n',mfilename,feFileToLoad)
    load(feFileToLoad);
    fprintf('[%s] Extracting info: \n%s\n ======================================== \n\n',mfilename,feFileToLoad)
    keyboard
    
end