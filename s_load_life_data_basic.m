%% Example load data for Cesar

% (1) Download the data from:
http://purl.stanford.edu/cs392kv3054

% (2) UNZIP/UNTAR the file.
% (3) Add the unzipped/untarred Data folder to your matlab search path. To do so in the MatLab prompt type:
addpath(genpath('/my/path/to/the/life_data_demo/folder/'))

% (4) Load the fe structure
feFileToLoad = fullfile('/path/to/fe-structure.mat')
load(feFileToLoad); 

% (5) Load the diffusion data
dwiFile       = fullfile(lifeDemoDataPath('diffusion'),'life_demo_scan1_subject1_b2000_150dirs_stanford.nii.gz');
dwi_nii = niftiRead(dwiFile);
% This file contains a field called "data" (dwi_nii.data), data stored the diffusion data as [x,y,z,dirs+nbvals], 
% where x,y,z are the spatial coordinates in the brain. the last dimension is the number of direction

% So the indices to the diffusion directions are:
idxDiff = 1:96;

% The indices to the S0 are 
idxSzero = 97:106;

% See (8) below for more information.

% (6) Store the affine trasformation for the fibers to index into the data.
xform.acpc2img = inv(dwi_nii.qto_xyz);

% (7) Extract the fibers wtih coordinates addresing the data.
fgFileName    = fullfile(lifeDemoDataPath('tractography'), ...
                'life_demo_mrtrix_csd_lmax10_probabilistic.mat');
fgAcpc = fgRead(fgFileName); % Fibers are by default stored in AC-PC orientation in our software.
% We will now transform them back into indices of the data 
fgImg  = dtiXformFiberCoords(fgAcpc, xform.acpc2img,'img'); % After this each x,y,z coordinate in fg.fiber{i} is an index in dwi_nii[x,y,x,:] NB it is not rounded

% After this you do not really need the fgAcpc any longer.
clear fgAcpc

% (8) Load the bvecs and bvals
% BVECS are the x,y,z coodinates of the unit vectors indicating the directions along a sphere where the diffusion measurements were acquired
bvecs = dlmread(lifeDemoDataPath('diffusion'),'life_demo_scan1_subject1_b2000_150dirs_stanford.bvecs')

% BVALS is a single number set by the scanner
bvals = dlmread(lifeDemoDataPath('diffusion'),'life_demo_scan1_subject1_b2000_150dirs_stanford.bvals')
% BVALS indices inside diw_nii.data(x,y,z,97:end); It is always 2000 for this dat aset but it coudl be different later. 1000, 2000 etc.



