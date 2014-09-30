%% A script showing how to create simulation ground-truth data
%
% The ground-truth data are built by the script s_make_dtiinit.m
%
% This 
%  FP took a data set from Stanford at 1.5mm 96-directions from the LiFE
%  paper and resampled to 4mm. 
%  Purpose: Smaller data size.  This is in the raw/ directory.
%
%  Then he run dtiInit on the data and created the init log and the
%  dti96trilin directory.  The additional files are produced by dtiInit as
%  well.
% 
%  Then he ran mrTrix to create a whole-brain connectome.  He requested
%  160K fibers and Lmax = 10, probabilistic.
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

%% Load the AFQ Fibers and merge them into a single

dataDir = fullfile(mrvDataRootPath,'life','data');
afqDir  = fullfile(dataDir,'AFQ');
fascicleList = dir(fullfile(afqDir,'*.mat'));

% The last one is the tracts_classification_indices, which we don't want in
% the list.  Typically, there are 20 tracts.
fascicleList = fascicleList(1:(end-1));

%% Merge all the fascicles into a single, whole-brain group

% When we merge, we may not take all the fascicles within a tract
sampF = 0.08;
for ii=1:length(fascicleList)
    fg = fgRead(fullfile(afqDir,fascicleList(ii).name));
    n = fgGet(fg,'n fibers');
    idx = randsample(1:n,ceil(sampF*n));
    fg = fgExtract(fg,idx,'keep');
    if ii==1
        fgAll = fg;
    else
        % fgMerge should be able to take a list of fullpath file names.
        fgAll = fgMerge(fgAll,fg);
    end 
end
fgAll = fgSet(fgAll,'name','AFQ Merged 4mm');

% Visualize now
% mbaDisplayConnectome(fgAll.fibers,mrvNewGraphWin);

%% Build the LiFE model

% Key files
dwiFile =  fullfile(dataDir,'run01_fliprot_aligned_trilin_aligned_trilin.nii.gz');
feName = 'afqModel'; 
feDir  = fullfile(dataDir,'LiFE');
t1File = fullfile(dataDir,'anatomy','t1.nii.gz');
feFullFname = fullfile(feDir,[feName,'.mat']);

if exist(feFullFname,'file')
    load(feFullFname);
else
    % Initialize the fascicle evaluation structure
    fe = feConnectomeInit(dwiFile,fgAll,feName,feDir,dwiFile,t1File);
    
    % Do the fit and store the answer in the fe structure
    fe = feSet(fe,'fit',feFitModel(feGet(fe,'mfiber'),feGet(fe,'dsigdemeaned'),'bbnnls'));
    
    % Save the result
    feConnectomeSave(fe)
end

%% Visualize the culled fibers, those with large weights

w = feGet(fe,'fiber weights');
% top20 = prctile(w,80);
lst = (w > 1e-6);
fgO = feGet(fe,'fibers acpc');
fgO = fgExtract(fgO,lst,'keep');

feName = 'afqModelCulled'; 
feFullFname = fullfile(feDir,[feName,'.mat']);

if exist(feFullFname,'file')
    load(feFullFname);
else
    % Initialize the fascicle evaluation structure
    fe = feConnectomeInit(dwiFile,fgAll,feName,feDir,dwiFile,t1File);
    
    % Do the fit and store the answer in the fe structure
    fe = feSet(fe,'fit',feFitModel(feGet(fe,'mfiber'),feGet(fe,'dsigdemeaned'),'bbnnls'));
    
    % Save the result
    feConnectomeSave(fe)
end

% Have a look at the most important fibers.  
% mbaDisplayConnectome(fgO.fibers,mrvNewGraphWin);

%%  Now, buid a simulation based on the culled AFQ fibers

% We will create a NIFTI that is built from the predicted diffusion of
% these fibers plus some of the mean signal in the original.  Since the
% fibers put in some mean, we will not put all of the mean back in, but
% maybe 50 or 80 percent of the mean.

% The nifti should be initialized with the header and size of the original
% DWI nifti. The anisotropic signal will be only at voxels of the
% fibers in fg. The anisotropic signal will be the contribution from each
% fiber node times the fiber weight. 
%
% We will then decide how to add in an isotropic term.  One possibility is
% to add in something proportional to the original NIFTI isotropic
% component.
%
% To do the above we will write a modified version of
%
%    feConnectomeBuildModel 
% 
% that does not remove the mean of the prediction signal, and that allows
% adding in an additional isotropic signal.
%
% Remember to check warning: 
% Warning: redundant nodes 
% > In feComputeVoxelSignal at 31
%   In feConnectomeBuildModel>(parfor body) at 83
%   In parallel_function>make_general_channel/channel_general at 870
%   In remoteParallelFunction at 30 

zeroMean = false;
fe = feConnectomeBuildModel(fe,zeroMean);
fe = feSet(fe,'name','ForwardModel');
nVoxels = feGet(fe,'n voxels');
v2fn = feGet(fe,'voxel 2 fn pair');

% Xfrom fibers in IMG then find the index in the fe structure
c = fe.fg.fibers{100};
plot3(c(1,:),c(2,:),c(3,:),'-o')

% Find the index in the fe structure
foundVoxels = feGet(fe,'find voxels',c');
sub2ind(12,31,24)

% These are the voxels in 'whichFiber'
voxList = [];
whichFiber = 375;
for ii=1:nVoxels
    if any(ismember(v2fn{ii}(:,1),whichFiber))
        voxList(end+1) = ii;
    end
end

voxCoords = feGet(fe,'roi coords');
voxCoords = voxCoords(voxList,:);
mrvNewGraphWin;
plot3(voxCoords(:,1),voxCoords(:,2),voxCoords(:,3),'-o');
axis on; grid on







