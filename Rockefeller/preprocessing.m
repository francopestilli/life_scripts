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

% Run FSL topup preprocessing of dti images.
% First we pull out all the B0 volumes:
% AP

% (1) Extract the B0 files from the NIFTI containing B0 and diffusion weighted
% imaging
niAP  = nifitLoad('path/to/niftiAP.nii.gz');
bvals = dlmread('path/to/bvals.bval'); 
B0_AP_indices = bvals==0; % THese are the indices to the B0 

% We repeate the same for the PA
niPA  = nifitLoad('path/to/niftiPA.nii.gz');
bvals = dlmread('path/to/bvals.bval'); 
B0_PA_indices = bvals==0; % These are the indices to the B0 

% (2) Combine the PA and AP B0 files into a signel nifti
niB0 = niPA; 
niB0.data = [];
niB0.data = cat(niPA.data(:,:,:,B0_PA_indices),niAP.data(:,:,:,B0_AP_indices));
% Modify the rest of the NIFTI file to mach the new file size: niB0.DOMORE
% Store the file name: fileName = '/path/to/file/name/B0_AP_and_PA.nii.gz'

% (3) Run TOPUP on the B0 images
% TOPUP is a FLS command that will remove distortions from the images due
% to the phase encode direction
cmd = sprintf('topup', ... 
              '--imain=/path/to/nifti/B0_AP_and_PA.nii.gz' , ...
              '--datain=/path/to/dti_acqparams.txt', ... 
              '--config=/path/to/topup/configuration/params/topup_b02b0.cnf', ...
              '--out=/path/to/output/file/B0_AP_and_PA_out.nii.gz', ...
              '--fout=/path/to/output/file/B0_AP_and_PA_TOPUP_fout', ...
              '--iout=/path/to/output/file/B0_AP_and_PA_TOPUP_iout');
mysystem(cmd)

% (4) Apply topup to the rest of the images, these are the diffusion weighted images
cmd = sprintf('applytopup --imain=/path/to/B0_AP,/path/to/B0_PA', ...
              '--topup=dti/%s_TOPUP', ...
              '--datain=/path/to/dti_acqparams.txt',...
              '--inindex=1,5', ...
              '--out=path/to/_hifib0');
mysystem(cmd)

% ** Franco make sure we need to do BET here: **
% (4.1) Run BET to create a brain_mask
cmd =sprintf('bet /path/to/B0_ dti/%s_hifib0_brain -R -m -f 0.2',subject,subject);
mysystem(cmd)

% (5) Use eddy to unwarp the images using the topup fields. This will take
%     a few hours...
cmd = sprintf('eddy --imain=/path/to//%s_Big4D', ...
              '--acqp=/biac4/wandell/users/sajina/code/dwi-vf-patients/dti_acqparams_SARA.txt',... 
              '--mask=/path/to/brain/mask/bet_brain_mask.nii.gz', ...
              '--index=/path/to/inde/for/AP/PA/Index.txt', ...
              '--bvals=/path/to/bvals', ...
              '--bvecs=/path/to/bvecs', ... 
              '--out=/path/to/eddy_unwarped_images', ...
              '--topup=/path/to/outputdirectory/TOPUP', ... 
              '--fwhm=5 --flm=quadratic');
mysystem(cmd)

% (5) Preprocess the data using mrDiffusion, and dtiInit
dwParams = dtiInitParams;
dwParams.bvecsFile = '/path/to/merged.bvec';
dwParams.bvalsFile = '/path/to/merged.bval';

% (5.1) Build a new file name for the dt6 folder
res    = num2str(resolution); 
idx    = strfind(res,'.'); 
if ~isempty(idx), res(2) ='p';end

% (5.2) Set file name of the input diffusion data
dwi_raw_file  = 'nhp01_paap_aligned_trilin';
dt6_dir_name  = sprintf('dt6_%s_%smm',dwi_raw_file(1:14),res);
dataDir       = fullfile('/write/correct/path/to/the/data/here');
dwRawFileName = fullfile(dataDir,dt6_dir_name, 'raw', sprintf('%s.nii.gz',dwi_raw_file));
if ~exist(dwRawFileName,'file'), error('No dwRaw file %s',dwRawFileName); end
t1FileName    = fullfile(dataDir,dt6_dir_name, 'anatomy','t1.nii.gz');
if ~exist(t1FileName,'file'), error('No T1 File %s',t1FileNames); end

% (5.3) Initialization parameters
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

% (6) Run the preprocessing, this will take a few hours...
[dtFile, outBaseDir] = dtiInit(dwRawFileName, t1FileName, dwp);
