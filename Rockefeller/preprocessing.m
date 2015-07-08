% Preprocessing stream to integrate probabilistic tractography and LiFE for Ilaria and Winrich. 
%
%
% Processing steps:
%  * The script starts with raw diffusion data (NIFTI-1 files)
%  * Preprocess the data using dtiInit (motion compensation, eddy currents correction alignment)
%  * Generate a whole-brain candidate connectome using mrTrix
%  * Build and fit a LiFE model of the whole-brain candidate connectome.
%  * Reduce the connectome by keeping onlyfibers with positive weight.
%  
% Dependencies:
% - github.com/brain-life/lifebid
% - github.com/vistalab/vistasoft
% - github.com/francopestilli/life_scripts
%
% Copyright 2015 Franco Pestilli Indiana University pestillifranco@gmail.com

%% Set parameters for the reconstruction and preprocessign of the diffusion data
% We assume that the data are saved as NIFTI-1 files.

% We decie a numbe of fascicles to tract in the brain.
nFascicles = 120000;

% Build a new file name for the dt6 folder
res    = num2str(resolution); 
idx    = strfind(res,'.'); 
if ~isempty(idx), res(2) ='p';end

% File name of the input diffusion data
dwi_raw_file  = 'nhp01_paap_aligned_trilin';
dt6_dir_name  = sprintf('dt6_%s_%smm',dwi_raw_file(1:14),res);
dataDir       = fullfile('/write/correct/path/to/the/data/here');
dwRawFileName = fullfile(dataDir,dt6_dir_name, 'raw', sprintf('%s.nii.gz',dwi_raw_file));
if ~exist(dwRawFileName,'file'), error('No dwRaw file %s',dwRawFileName); end
t1FileName    = fullfile(dataDir,dt6_dir_name, 'anatomy','t1.nii.gz');
if ~exist(t1FileName,'file'), error('No T1 File %s',t1FileNames); end

% Initialization parameters
dwp = dtiInitParams;
dwp.eddyCorrect    = false;
dwp.phaseEncodeDir = 2;
dwParams.rotateBvecsWithCanXform = 1;
% We first set the resoltuion for the preprocessed files, this resolution
% (Ideally the resolutuon could be the closest isotropic resolution to the
% acquisition resolution.)
dwp.dwOutMm = [0.250 0.250 0.254]; % For the NHP data collected by Henning Voss acquisition resolution was [0.250, 0.250, 0.254] 
dwp.dt6BaseName = '';
dwp.bvecsFile   = fullfile(dataDir,dt6_dir_name, 'raw',sprintf('%s.bvecs',dwi_raw_file));
dwp.bvalsFile   = fullfile(dataDir,dt6_dir_name, 'raw',sprintf('%s.bvals',dwi_raw_file));

% Run the preprocessing
[dtFile, outBaseDir] = dtiInit(dwRawFileName, t1FileName, dwp);
