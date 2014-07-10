function fibersPDB = mrtrix_track_between_rois
%
% This functions shows how to track between two ROIS using mrtrix.
% This s very helpful for ideintifying some fiber groups for example the
% optic radiation.
%
% This is how te code works.
% 1. We load two ROIs in the brai, for Shumpei's project for example we
%    will load the right-LGN and the right-Visual cortex
% 2. We create union ROI by combining these two ROIs. The union ROI is used
%    as seeding fro the fibers. mrtrix will initiate and terminate fibers only
%    within the volume defined by the Union ROI.
% 3. We create a white matter mask. THis mask is generally a large portion
%    of the white matter. A portion that contains both union ROIs. For example
%    the right hemisphere.
% 4. We use mrtrix to track between the right-LGN and righ-visual cortex.
% mrtrix will initiate fibers by seeding within the UNION ROI and it will
% only keep fibers that have paths within the white matter masks.
%
% The final result of this script is to generate lot's of canddte fibers 
% that specifically end and start from the ROI of interest. Thisis an 
% approach similar to Contrack. 
%
% INPUTS: none
% OUTPUTS: the finela name of the ROI created at each iteration
%
% Written by Franco Pestilli (c) Stanford University Vistasoft

baseDir = '/dt6/directory/';

subjectDir = '/somedirectory/';
dtFile = fullfile(baseDir, subjectDir, '/path2/dt6.mat');
refImg = fullfile(baseDir, subjectDir, '/t1/t1_acpc.nii.gz');
fibersFolder = fullfile(baseDir, subjectDir, '/save/fibers/here');

% We want to track the cortical pathway (LGN -> V1/V2 and V1/V2 -> MT)
fromRois = {'r_LGN'};
toRois   = {'V1'};

% Set upt the MRtrix trakign parameters
trackingAlgorithm = {'prob'};
lmax    = [2]; % The appropriate value depends on # of directions. For 32, use lower #'s like 4 or 6. For, 6 or 10 is good [10];
nSeeds  = 5000; % 10000; 
nFibers = 50000; %1000000;
wmMask  = [];

% Make an (include) white matter mask ROI. This mask is the smallest
% set of white matter that contains both ROIS (fromRois and toRois)
%
% We use a nifti ROi to select the portion of the White matter to use for
% seeding
wmMaskName      = fullfile(baseDir,  subjectDir, '/ROIs/front_mask'); 
[~, wmMaskName] = dtiRoiNiftiFromMat(wmMaskName,refImg,wmMaskName,1);

% Then transform the niftis into .mif
[p,f,e] = fileparts(wmMaskName );
wmMaskMifName    = fullfile(p,sprintf('%s.mif',f)); 
wmMaskNiftiName  = sprintf('%s.nii.gz',wmMaskName);
mrtrix_mrconvert(wmMaskNiftiName, wmMaskMifName);    

% This first step initializes all the files necessary for mrtrix.
% This can take a long time.
files = mrtrix_init(dtFile,lmax,fibersFolder,wmMask);

% Some of the following steps only need to be done once for each ROI,
% so we want to do some sort of unique operation on the from/toRois
individualRois = unique([fromRois, toRois]);

% Convert the ROIs from .mat or .nii.gz to .mif format.
for i_roi = 1:length(individualRois)
    if exist(fullfile(p, [individualRois{i_roi}, '.nii.gz']),'file')
        thisroi = fullfile(p, [individualRois{i_roi}, '.nii.gz']);
 
    elseif  exist(fullfile(p, [individualRois{i_roi}, '.nii']),'file')
        thisroi = fullfile(p, [individualRois{i_roi}, '.nii']);
        
    elseif   exist(fullfile(p, [individualRois{i_roi}, '.mat']),'file')
         thisroi = fullfile(p, [individualRois{i_roi}, '.mat']);
    end
    
    mrtrix_roi2mif(thisroi,refImg);
end
    
% Create joint from/to Rois to use as a mask
for nRoi = 1:length(fromRois)
    % MRTRIX tracking between 2 ROIs template.
    roi{1} = fullfile(baseDir, subjectDir, '/ROIs/', fromRois{nRoi});
    roi{2} = fullfile(baseDir, subjectDir, '/ROIs/', toRois{nRoi});
    
    roi1 = dtiRoiFromNifti([roi{1} '.nii.gz'],[],[],'.mat');
    roi2 = dtiRoiFromNifti([roi{2} '.nii.gz'],[],[],'.mat');
    
    % Make a union ROI to use as a seed mask:
    % We will generate as many seeds as requested but only inside the voume
    % defined by the Union ROI.
    %
    % The union ROI is used as seed, fibers will be generated starting ONLy
    % within this union ROI.
    roiUnion        = roi1; % seed union roi with roi1 info
    roiUnion.name   = ['union of ' roi1.name ' and ' roi2.name]; % r lgn calcarine';
    roiUnion.coords = vertcat(roiUnion.coords,roi2.coords);
    roiName         = fullfile(baseDir, subjectDir, '/ROIs/',[roi1.name '_' roi2.name '_union']);
    [~, seedMask]   = dtiRoiNiftiFromMat(roiUnion,refImg,roiName,1);
    seedRoiNiftiName= sprintf('%s.nii.gz',seedMask);
    seedRoiMifName  = sprintf('%s.mif',seedMask); 
    
    % Transform the niftis into .mif
    mrtrix_mrconvert(seedRoiNiftiName, seedRoiMifName);
        
    % We cd into the folder where we want to sae the fibers.
    cd(fibersFolder);
    
    % We geenrate and save the fibers in the current folder.
    [fibersPDB{nRoi}, status, results] = mrtrix_track_roi2roi(files, [roi{1} '.mif'], [roi{2} '.mif'], ...
        seedRoiMifName, wmMaskMifName, trackingAlgorithm{1}, ...
        nSeeds, nFibers);
    
    % fgWrite(fibersPDB,['fibername'],'pwd')
end

return
