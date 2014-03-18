function s_fe_make_wm_mask
%
% This script makes the white-matter mask used to track the connectomes in
% Pestilli et al., LIFE paper.
%
% Retinotopy folder: /biac4/wandell/data/anatomy/
%
% Written by Franco Pestilli (c) Stanford University, Vista Team 2013

% Get the base directory for the data
datapath = '/home/frk/2t1/predator/';
subjects = {...   
    'KK_96dirs_b2000_1p5iso', ...
    'JW_96dirs_b2000_1p5iso', ... 
    'HT_96dirs_b2000_1p5iso', ...
    'KW_96dirs_b2000_1p5iso', ...
    'FP_96dirs_b2000_1p5iso', ...
    'MP_96dirs_b2000_1p5iso', ...
    };

if notDefined('saveDir'),      savedir      = fullfile(datapath,'anatomy');end
if notDefined('trackingType'), trackingType = 'lmax10';end
if notDefined('hemisphere'),   hemisphere   = {'left','right'};end
if notDefined('plotAnatomy'),  plotAnatomy  = 0;end
md_percentile = 95;

for isbj = 1:length(subjects)
   % Load the class file
   classFile  = fullfile(datapath,subjects{isbj},'anatomy', 't1_class.nii.gz'); 
   wmMaskFile = fullfile(datapath,subjects{isbj},'anatomy','wm_mask.nii.gz'); 
   mdFile     = fullfile(datapath,subjects{isbj},'dtiInit','dt6.mat');
   wmMask     = feMakeWMmask(classFile,mdFile,wmMaskFile,md_percentile);
end

end % Main function
