% Directory structure
%
%                 S0X
%        __________|__________
%       |                     |
%      raw                 anatomy
%       |                     |
% SOX_dwi.nii.gz         S0X_t1.nii.gz
% SOX_dwi.bvec
% S0X_dwi.bval
%


dir = '/home/samir/CNI_Lab/Imaging/Tractography'

addpath(genpath('/opt/vistasoft/'));
addpath(genpath('/opt/life/'));
addpath(genpath(fullfile(dir, 'life_scripts')));

resolution = 1.5;

dwi_file = 'S0X_dwi';
	dataDir = fullfile(dir, 'Data/S0X');

dwRawFileName = fullfile(dataDir, 'raw', sprintf('%s.nii.gz', dwi_file));
t1FileName = fullfile(dataDir, 'anatomy', 'S0X_t1.nii.gz');

dwp = dtiInitParams;
dwp.eddyCorrect = false;
dwp.phaseEncodeDir = 2;
dwp.rotateBvecsWithCanXform = 1;
dwp.dwOutMm = [resolution resolution resolution];
dwp.dt6BaseName = '';
dwp.bvecsFile = fullfile(dataDir, 'raw', sprintf('%s.bvec', dwi_file));
dwp.bvalsFile = fullfile(dataDir, 'raw', sprintf('%s.bval', dwi_file));

[dtFile, outBaseDir] = dtiInit(dwRawFileName, t1FileName, dwp);

% MATLAB OUTPUT
% >> preproc

% dir =

% /home/samir/CNI_Lab/Imaging/Tractography

% Loading raw data...
% Skipping eddy-current correction. Rigid-body motion correction.
% Data will be saved to: /home/samir/CNI_Lab/Imaging/Tractography/Data/S0X 
% Dims = [107 141 101 106] 
% Data Dir = /home/samir/CNI_Lab/Imaging/Tractography/Data/S0X/raw 
% Output Dir = /home/samir/CNI_Lab/Imaging/Tractography/Data/S0X 
% t1FileName = /home/samir/CNI_Lab/Imaging/Tractography/Data/S0X/anatomy/S0X_t1.nii.gz;
% freq_dim not set correctly in NIFTI header.
% phase_dim not set correctly in NIFTI header.
% Undefined function 'spm_get_defaults' for input arguments of type 'char'.

% Error in dtiRawComputeMeanB0 (line 14)
% estParams        = spm_get_defaults('coreg.estimate');

% Error in dtiInit (line 155)
% if computeB0, dtiRawComputeMeanB0(dwRaw, bvals, dwDir.mnB0Name); end

% Error in preproc (line 37)
% [dtFile, outBaseDir] = dtiInit(dwRawFileName, t1FileName, dwp);
 
