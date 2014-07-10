function s_pestilli_etal_fig8_major_fascicles_save_fg()
%
% load an plots the 20 major fascicles generated with AFQ saves them out as
% independent Fiber Groups
%
% Copyright Franco Pestilli 2014 Stanford University

% Get the base directory for the data
datapath = '/marcovaldo/frk/2t1/predator/';
subjects = {...
    'FP_96dirs_b2000_1p5iso', ...'JW_96dirs_b2000_1p5iso', ...
    'KK_96dirs_b2000_1p5iso', ...
    'MP_96dirs_b2000_1p5iso', ...
    'HT_96dirs_b2000_1p5iso', ... 
    'KW_96dirs_b2000_1p5iso', ...
    };

for isbj = 1:length(subjects)
    fasciclesPath  = fullfile(datapath,subjects{isbj},'afq');
    fasciclesFiles = dir(fullfile(fasciclesPath,'fascicles_*recomputed.mat'));
    savedir        = fullfile(fasciclesPath,'major_fascicles');
    mkdir(savedir);
    fprintf('[%s] Subject: %s\n',mfilename,subjects{isbj})
    load(fullfile(fasciclesPath,fasciclesFiles(1).name))
    for iFas = 1:length(fascicles)
        classification.names{iFas}(classification.names{iFas}==' ')='_';
        tract_name = classification.names{iFas};
        fprintf('[%s] Saving Fascicle: %s\n',mfilename,tract_name)
        fgWrite(fascicles(iFas),fullfile(savedir,tract_name),'mat');
    end
    clear fascicles classification
end

end % Main function

