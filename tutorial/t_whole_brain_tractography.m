%% t_whole_brin_tractography.m
%
% This is the script used to preprocess diffusion data using the VISTASOFT tools and generate a whole brain set of fascicles using mrTrix.
%
% This is data colelcted at Stanford University and published with Petilli et al., Nature Methods 2014.
%
% The following are the steps necessary to rung the code"
% A. Download LiFE, life_scripts, VISTASOFT and MBA from github.com
% B. Start MatLab.
% C. Add LiFE, life_scripts, VISTASOFT and MBA (Matlab Brain Anatomy).
% D. Make sure that MRTRix 0.2.X is installed in yoru OS System on Ubuntu and Mint (apt-get install mrtrix)
%
% Dependencies:
%  - VISTASOFT: https://github.com/vistalab/vistasoft/
%  - LiFE: https://github.com/francopestilli/life
%  - life_scripts: https://github.com/francopestilli/life_scripts
%  - MBA: https://github.com/francopestilli/mba
%  - MRTRIX: http://www.brain.org.au/software/mrtrix/ 
%  - For practice: (Data to practice with this can be found here: https://stacks.stanford.edu/file/druid:cs392kv3054/life_demo_data.tar.gz)
% 
% The following lines are necessary to set up the path in MatLab to make all the functions visible. 
%    >> addpath(genpath('/my/path/to/the/VISTASOFT/folder/'))
%    >> addpath(genpath('/my/path/to/the/life_data_demo/folder/'))
%    >> addpath(genpath('/my/path/to/the/life/folder/'))
%    >> addpath(genpath('/my/path/to/the/life_scripts/folder/'))
%
% After this the code below should work and walk thorught the main steps. 
%
% Copyright 2015 Franco Pestilli Indiana University, pestillifranco@gmail.com

%% (0) Make sure the files are organized as necessary (file names are relative):
%    <subject directory>
%       _____|_____
%      |           |  
%  t1.nii.gz      raw 
%                  |
%               dwi.nii.gz
%               dwi.bval
%               dwi.bvec

%% (1) Initianlize a set of folders and files.
% We will compute the diffusion tensor for a single brain data set (NIFTI)
% 

% cube mm of the spatial resolution for the output diffusion data
resolution = 2; % This resolution should be isotropic and (in most cases) match the acquisition resolution

% File name of the input diffusion data. Match file names below to the ones used in setting up the folder strucutre above (0)
dt6_dir_name  = sprintf('dt6_%s_%imm',dwi_raw_file(1:14),res);
dwRawFileName = fullfile(dataDir,dt6_dir_name, 'raw', sprintf('%s.nii.gz',dwi_raw_file));
t1FileName    = fullfile(dataDir,dt6_dir_name, 'anatomy','t1.nii.gz');

% Initialization parameters. These parameters will be used for intializing the robust tensor fit.
dwp = dtiInitParams; % Help dtiInit.m
dwp.eddyCorrect    = false; % False if you collected a rephase sequence or if your data were previously AP/PA combined using FSL.
dwp.phaseEncodeDir = 2; % Make sure this is correct
dwp.rotateBvecsWithCanXform = 1;
dwp.dwOutMm = [resolution resolution resolution];
dwp.dt6BaseName = '';
dwp.bvecsFile   = fullfile(dataDir,dt6_dir_name, 'raw',sprintf('%s.bvecs',dwi_raw_file));
dwp.bvalsFile   = fullfile(dataDir,dt6_dir_name, 'raw',sprintf('%s.bvals',dwi_raw_file));

% Run the actuall preprocessing.
% This will motion correct, eddy current correct (if requested) and fit the tensor.
[dtFile, outBaseDir] = dtiInit(dwRawFileName, t1FileName, dwp);

%% (2) Run MRtrix tractography
tic, fibersFolder = 'mrtrix_fascicles';
if ~exist(fibersFolder,'dir'), mkdir(fibersFolder); end
dt6_file = fullfile(dataDir,dt6_dir_name,'dti96trilin','dt6.mat');
nFascicles = 120000;
tractographyType = {'prob'}; % THis can be probabilistic and deterministic
lmax = mrtrix_findlmax(48); % Make sure that the number of directiosn here is set to the acquired number of directions
[status, ~, fg] = feTrack(tractogrpahyType, dt6_file,fibersFolder,,nFascicles);
toc
