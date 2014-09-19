%% A script showing how to perform a virtual lesion
%
% Processing steps:
%  * The script starts with raw diffusion data
%  * Passes it through the dtiInit preprocessing steps
%  * Tracks using mrTrix
%  * Uses AFQ to identify tracts of interest
%  * Sets up a Life structure
% 
%
% FP Vistasoft Lab, Stanford, 2014

%% Set the path and run dtiInit

% cube mm of the spatial resolution for the output diffusion data
resolution = 4; 

% Build the new file name for the dt6 folder
res    = num2str(resolution); 
idx    = strfind(res,'.'); 
if ~isempty(idx), res(2) ='p';end

% File name of the input diffusion data
%dwi_raw_file  = 'dwi_data_b2000_1p25mm';
dwi_raw_file  = 'run01_fliprot_aligned_trilin';
dt6_dir_name  = sprintf('dt6_%s_%smm',dwi_raw_file(1:14),res);
%dataDir       = '/home/frk/2t1/HCP/105115_data_variability/';
dataDir       = fullfile(mrvDataRootPath,'life');
dwRawFileName = fullfile(dataDir,dt6_dir_name, 'raw', sprintf('%s.nii.gz',dwi_raw_file));
if ~exist(dwRawFileName,'file'), error('No dwRaw file %s',dwRawFileName); end
t1FileName    = fullfile(dataDir,dt6_dir_name, 'anatomy','t1.nii.gz');
if ~exist(t1FileName,'file'), error('No T1 File %s',t1FileNames); end

% Initialization parameters
dwp = dtiInitParams;
dwp.eddyCorrect    = false;
dwp.phaseEncodeDir = 2;
%dwParams.rotateBvecsWithCanXform = 1;
  
% Set the spatial resolution of the output diffusion data
dwp.dwOutMm = [resolution resolution resolution];

dwp.dt6BaseName = '';
dwp.bvecsFile   = fullfile(dataDir,dt6_dir_name, 'raw',sprintf('%s.bvecs',dwi_raw_file));
dwp.bvalsFile   = fullfile(dataDir,dt6_dir_name, 'raw',sprintf('%s.bvals',dwi_raw_file));

% Run the preprocessing
[dtFile, outBaseDir] = dtiInit(dwRawFileName, t1FileName, dwp);

%% Run MRtrix tractography
tic, fibersFolder = 'mrtrix_fascicles';
if ~exist(fibersFolder,'dir'), mkdir(fibersFolder); end
dt6_file = fullfile(dataDir,dt6_dir_name,'dti96trilin','dt6.mat');
nFascicles = 120000;
[status, ~, fg] = feTrack({'prob'}, dt6_file,fibersFolder,mrtrix_findlmax(90),nFascicles);
toc
       
%% Run AFQ

% Segment the fibers using AFQ
tic, [fg_classified,~,classification]= AFQ_SegmentFiberGroups(dt6_file, fg);
toc

%% Split the fiber groups into individual groups
tic, fascicles = fg2Array(fg_classified);
toc

% Write the fascicles down to disk as independent files
tic, afqFolder = fullfile(dataDir,dt6_dir_name, 'AFQ');
if ~exist(afqFolder,'dir'), mkdir(afqFolder); end
for iif = 1:length(fascicles)
    fgWrite(fascicles(iif),    fullfile(afqFolder,[fascicles(iif).name,'_uncleaned']),'mat')
end

% Clean the fibers, we apply the same trhesholds to all fiber
% groups this is the default thrshold used by AFQ. This is done by
% not passing opts
tic, [fascicles1, classification] = feAfqRemoveFascicleOutliers(fascicles,classification);
toc

% Write the fascicles down to disk as independent files
tic, afqFolder = fullfile(dataDir,dt6_dir_name, 'AFQ');
if ~exist(afqFolder,'dir'), mkdir(afqFolder); end
for iif = 1:length(fascicles)
    fgWrite(fascicles1(iif),    fullfile(afqFolder,fascicles(iif).name),'mat')
end

% Save the segemented fascicles and the indices into the Mfiber
classFile2Save = fullfile(afqFolder,'tracts_classification_indices');
save(classFile2Save,'fg_classified','classification','fascicles1')

%% Build the LiFE model.
dwiFile =  '/home/frk/2t1/fp96_data_variability/dt6_run01_fliprot__4mm/run01_fliprot_aligned_trilin_aligned_trilin.nii.gz';
feFileName = ['fe_structure_',dt6_dir_name];
tic, feFolder = fullfile(dataDir,dt6_dir_name, 'LiFE');
if ~exist(feFolder,'dir'), mkdir(feFolder); end
fe = feConnectomeInit(dwiFile,fg,feFileName,feFolder,dwiFile,t1FileName);
fe = feSet(fe,'fit',feFitModel(feGet(fe,'mfiber'),feGet(fe,'dsigdemeaned'),'bbnnls'));

% Perform a virtual lesion
