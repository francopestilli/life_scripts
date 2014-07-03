% Fiber cup directory
baseDir = '/biac4/wandell/data/diffusion/fibercup';
saveDir = '/biac4/wandell/data/diffusion/fibercup/life';
feFileName = 'fibercup_fe';
% Change to base directory
cd(baseDir);

% Read in the dwi data
dwi_path  = fullfile(baseDir,'dwi-b0650.nii.gz');
dwi_path2 = fullfile(baseDir,'dwi-b0650.nii.gz');
mask_path = fullfile(baseDir,'mask.nii.gz');

% Read the basic data fiels and fix the header
dwi     = readFileNifti(dwi_path);
b0      = dwi; 
b0.data = b0.data(:,:,:,1);
b0.ndim = 3; 
b0.pixdim = [3 3 3]; 
b0.dim    = b0.dim(1:3);

% Read in the ground truth fibers
fg_truth_path = fullfile(baseDir,'ground_truth','3mm','ground_truth_3mm.trk');
fg_truth = read_trk_to_fg(fg_truth_path);
for ii = 1:length(fg_truth.fibers)
    fg_truth.fibers{ii} = fg_truth.fibers{ii}.*repmat([-1; -1; -1],1,size(fg_truth.fibers{ii},2));
end

% Load some tracking fibers:
fg_prob_path = fullfile(baseDir,'tracking','csd_prob_track_ground_truth_endpoints.trk');
fg_prob = read_trk_to_fg(fg_prob_path);
xform = [[eye(3); 0 0 0]'; b0.qto_xyz(:,end)']';
fg_prob_acpc = dtiXformFiberCoords(fg_prob, xform);

% Add the group truth
%fg = fgMerge(fg_prob_acpc,fg_truth);

%% LiFE
if ~exist(saveDir,'dir'), mkdir(saveDir);end

% Intialize a local matlab cluster if the parallel toolbox is available.
feOpenLocalCluster;

% Initialize the Connectome
fe = feConnectomeInit(dwi_path,fg_prob_acpc,feFileName,saveDir,dwi_path2);
if cull == 1
    feC = feConnectomeCull(fe);
    w = feGet(feC,'fiberweights');
    fgCulled = feGet(feC,'fibers acpc');
    fgCulled = fgExtract(fgCulled,find(w > 10^-4));
else
    M = feGet(fe,'mfiber'); % get the life model out
    dsig = feGet(fe,'dsigdemeaned'); % Get the diffusion signal
    
    % Fit the model and cull. This will take some time...
    fefit = feFitModel(M,dsig,'bbnnls');
    fe    = feSet(fe,'fit',fefit);
    w     = feGet(fe,'fiberweights');
    fgCulled = fgExtract(fg_prob_acpc,find(w > 0));
end

% Render original fibers and data
fh = figure('name','original fibers','color','k');
mbaDisplayConnectome(fg_prob_acpc.fibers,fh,[.9 .4 .2],'single',[], [], .44);
hold on
mbaDisplayBrainSlice(b0,[0 0 1], [], [],.4);

% Render culled fibers and data
fh = figure('name','optimized fibers','color','k');
mbaDisplayConnectome(fgCulled.fibers,fh,[.2 .6 .9],'single',[], [], .44);
hold on
mbaDisplayBrainSlice(b0,[0 0 1], [], [],.4);

fh = figure('name','dwi signal','color','k');
mbaDisplayBrainSlice(b0,[0 0 2], [], [],1);
view(0,90);
box off;axis off;
set(gcf,'Color',[0 0 0]);
set(gca,'color',[0 0 0])
hold off;
set(fh, ...
    'Units','normalized', ...
    'Position', [0 0.2 0.2 0.4], ...
    'DoubleBuffer', 'on', ...
    'Color',[0 0 0]);
daspect ([ 1 1 1 ]);

save('fiber_cup_life_results.mat','fe','feC')

%% Make a plot of the weights:
figName = sprintf('Fascicle weights');
mrvNewGraphWin(figName);
[y,x] = hist(w(w > 0),logspace(-4,-.3,50));
semilogx(x,y)
set(gca,'tickdir','out','fontsize',16,'box','off')
title('fascicle weights','fontsize',16)
ylabel('number of fascicles','fontsize',16)
xlabel('fascicle weight','fontsize',16)


%% Make a plot of the RMSE:
rmse   = feGet(fe,'vox rmse');
rmsexv = feGetRep(fe,'vox rmse');
figName = sprintf('RMSE');
mrvNewGraphWin(figName);
% Non-cross-validated
[y,x] = hist(rmse,50);
plot(x,y,'k-')
hold on
% Cross-validated
[y,x] = hist(rmsexv,50);
plot(x,y,'r-')
set(gca,'tickdir','out','fontsize',16,'box','off')
title('Root-mean squared error distribution across voxels','fontsize',16)
ylabel('number of voxels','fontsize',16)
xlabel('rmse','fontsize',16)
legend({'RMSE fitted data set','RMSE cross-validated'},'fontsize',16)

%% Explained variance
r2   = feGet(fe,'r2vox');
r2xv = feGetRep(fe,'r2vox');
figName = sprintf('R2');
mrvNewGraphWin(figName);
% Non-cross-validated
[y,x] = hist(r2,50);
plot(x,y,'k-')
hold on
% Cross-validated
[y,x] = hist(r2xv,50);
plot(x,y,'r-')
set(gca,'tickdir','out','fontsize',16,'box','off')
title(sprintf('Explained variance (R2 %2.2f median, R2 %2.2f mean)',median(r2),mean(r2)),'fontsize',16)
ylabel('number of voxels','fontsize',16)
xlabel('R2','fontsize',16)
legend({'R2 fitted data set','R2 cross-validated'},'fontsize',16)



keyboard