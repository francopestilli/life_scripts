function s_fe_mt_ips_tract_hcp(hemisphere)
%
% This script performs a test of conenctivity of MT+ (Zilles t al ROI from 
% freesurfer) with the Superior parietal Cortex (Aparc Freesurfer segementation). 
% The following are the steps we perform:
%  - It loads a whole-brain tractography solution. 
%  - It loads the MTand Parietal ROI. 
%  - It finds the tract connecting the two ROIs.
%  - It builds a connectome model in the ROI defined by the tract
%  - It performs a virtual lesion on the tract connecting MT and parietal.
%
% Copyright by Franco Pestilli Stanford University, 2014

dataDir = {'2t1','2t1'};
subj_count = 0;
for iDd = 1:length(dataDir)
% Get the base directory for the data
[~,hostname] = system('hostname');
hostname = deblank(hostname);
switch dataDir{iDd}
    case {'2t2'}
        switch hostname
            case {'marcovaldo'}
                datapath = '/home/frk/2t2/HCP/';
            otherwise
                datapath = '/marcovaldo/frk/2t2/HCP/';
                
        end
        subjects = {...
            '118730', ...
            '115320', ...
            '117122', ...
            };
        
    case {'2t1'}
        switch hostname
            case {'marcovaldo'}
                datapath = '/home/frk/2t1/HCP/';
            otherwise
                datapath = '/marcovaldo/frk/2t1/HCP/';
                
        end
        subjects = {...
            '111312', ...
            '113619', ...
            '105115', ...
            '110411', ...
            };
    otherwise
        keyboard
end

if notDefined('saveDir'), savedir = fullfile('/marcovaldo/frk/Dropbox','pestilli_etal_revision',mfilename);end
if notDefined('trackingType'), trackingType = 'lmax10';end
if notDefined('hemisphere'), hemisphere = {'left','right'};end
if notDefined('plotAnatomy'), plotAnatomy = 1;end
anatomyPath     = '/dev/data/anatomy/';
parietalRoiName = 'lh_MT_label_smooth3mm.nii.gz';
mtRoiName       = 'rh_MT_label_smooth3mm.nii.gz';
parietalRoiName = 'lh_superiorparietal_label_smooth3mm.nii.gz';
mtRoiName       = 'rh_superiorparietal_label_smooth3mm.nii.gz';
            
for iSbj = 1:length(subjects)
% Load the FE structure
saveDir        = fullfile(savedir,subjects{iSbj});
fibergroupPath = fullfile(datapath,subjects{iSbj},'fibers');
fgFileToLoad   = dir(fullfile(fibergroupPath,sprintf('*%s*2000*.pdb',trackingType)));
fname          = fgFileToLoad(1).name;
fgFileToLoad   = fullfile(fibergroupPath,fname);
fgFileToLoad1   = dir(fullfile(fibergroupPath,sprintf('*%s*2000*.pdb','lmax8')));
fname1          = fgFileToLoad1(1).name;
fgFileToLoad1   = fullfile(fibergroupPath,fname1);
fgFileToLoad2   = dir(fullfile(fibergroupPath,sprintf('*%s*2000*.pdb','lmax2')));
fname2          = fgFileToLoad2(1).name;
fgFileToLoad2   = fullfile(fibergroupPath,fname2);

subj_count = subj_count+1;
for ih = 1:length(hemisphere)
    % Handle parallel computing
    if matlabpool('size') == 0
        c = parcluster;
        c.NumWorkers = 12;
        matlabpool(c);
    end
    fprintf('[%s] Loading: \n%s\n ======================================== \n\n',mfilename,fgFileToLoad)
    fg = fgRead(fgFileToLoad);
    fg = fgMerge(fg,fgRead(fgFileToLoad1));
    fg = fgMerge(fg,fgRead(fgFileToLoad2));

    % Set all the variables that depend on the hemisphere
    switch hemisphere{ih}
        case {'left'}
            parietalRoiName = 'lh_superiorparietal_label_smooth3mm.nii.gz';
            mtRoiName       = 'lh_MT_label_smooth3mm.nii.gz';
            axisLims   = [-67 -18 -110 -40 -18 80];
            vw         = [-75,30];
            slices     = {[-18 0 0],[0 -40 0],[0 0 -14 ]};
            lght       = 'left';
            SLaxLims   = [-55 2 -120 120 -20 40 ];
            histcolor{1}  = [0.4 0.4 0.4];
            histcolor{2}  = [.6 0.4 0.4];
            
        case {'right'}
            parietalRoiName = 'rh_superiorparietal_label_smooth3mm.nii.gz';
            mtRoiName       = 'rh_MT_label_smooth3mm.nii.gz';
            axisLims   = [18 67 -110 -40 -18 80];
            vw         = [75,30];
            slices     = {[18 0 0],[0 -40 0],[0 0 -14 ]};
            lght       = 'right';
            SLaxLims   = [-2 55 -120 120 -20 40 ];
            histcolor{1}  = [0 0 0];
            histcolor{2}  = [.8 0.4 0.4];
            
        otherwise
            keyboard
    end
    
    % Load the ROIs  
    FS_SUBJECT       = matchSubject2FSSUBJ(subjects{iSbj});
    roiDir           = fullfile(anatomyPath,FS_SUBJECT,'label');
    mtFileName       = fullfile(roiDir,mtRoiName);
    parietalFileName = fullfile(roiDir,parietalRoiName);
    
    % Find the fascicles in the connectome that touch both MT+ and parietal.
    mt       = dtiImportRoiFromNifti(mtFileName,[parietalFileName(1:end-7),'_ROI.mat']);
    parietal = dtiImportRoiFromNifti(parietalFileName,[parietalFileName(1:end-7),'_ROI.mat']);
    tic, fprintf('\n[%s] Segmenting tract from connectome... \n',mfilename)
    [mt2parietalTract, ~] = feSegmentFascicleFromConnectome(fg, {mt,parietal}, {'endpoints','endpoints'}, 'mt_parietal');
    
    % Clean the fibers by length, fibers that too long are likely to go far
    % frontal and not just touch MT+ and parietal.
    [~, keep]        = mbaComputeFibersOutliers(mt2parietalTract,2,2); 
    mt2parietalTract = fgExtract(mt2parietalTract,find(keep),'keep');
    toc
    
    % Find the Coordinates of the mt-parietal tract
    tic, fprintf('\n[%s] Create ROI from MT-Parietal tract... \n',mfilename)
    tractRoi = dtiCreateRoiFromFibers(mt2parietalTract);toc
    tic, fprintf('\n[%s] Removing fibers not going throught the tractROI... \n',mfilename)
    [fg,~, ~, ~] = dtiIntersectFibersWithRoi([],'and',2,tractRoi,fg);
    fg = feClipFibersToVolume(fg,tractRoi.coords,1);toc
    
    % Build LiFE model only in this volume, fit, cull
    dwiPath  = fullfile(datapath,subjects{iSbj},'diffusion_data');
    dwiFiles = dir(fullfile(dwiPath,sprintf('*2000_aligned_trilin*.gz')));
    dwiFile       = fullfile(dwiPath,dwiFiles(1).name);
    dwiFileRepeat = fullfile(dwiPath,dwiFiles(1).name);
    t1File        = fullfile(datapath,subjects{iSbj},'anatomy','T1w_acpc_dc_restore_1p25.nii.gz');
    
    % Directory where to save the fe structures
    saveDirC    = fullfile(datapath,subjects{iSbj},'connectomes');
    feFileName = [parietal.name, '_', fname(1:40), '.mat'];
    fe   = feConnectomeInit(dwiFile,fg,feFileName,saveDirC,dwiFileRepeat,t1File);
    M    = feGet(fe,'mfiber');
    dSig = feGet(fe,'dsigdemeaned');
    fit  = feFitModel(M,dSig,'bbnnls');
    fe   = feSet(fe,'fit',fit);clear fit
    fg   = feGet(fe,'fibers acpc');
    [~, keepFascicles] = feSegmentFascicleFromConnectome(fg, {mt,parietal}, {'endpoints','endpoints'}, 'mt_parietal');
    
     % Perform a virtual lesion: MT+ and parietal.
    display.tract = true;
    display.distributions = true;
    display.evidence = true;
    [SE(subj_count,ih), fig] = feVirtualLesion(fe,keepFascicles,display);    
    clear fe
    saveFig(fig(1).fh,fullfile(saveDir,[fig(1).name, '_',hemisphere{ih}]),'eps')
    saveFig(fig(2).fh,fullfile(saveDir,[fig(2).name(1:end-4), '_',hemisphere{ih}]),'eps')   
    if plotAnatomy             
        % Load the T1 file for display
        t1     = niftiRead(t1File);
        
        for ifs = 5:-1:3
        % Show te new fiber group
        figure(fig(ifs).fh); hold on
        h  = mbaDisplayBrainSlice(t1, slices{1});
        h  = mbaDisplayBrainSlice(t1, slices{2});
        h  = mbaDisplayBrainSlice(t1, slices{3});
        view(vw(1),vw(2)); axis(axisLims);
        %set(gcf,'Position',[0.0148 0.0148 .35 .87])
        saveFig(fig(ifs).fh,fullfile(saveDir,[fig(ifs).name, '_',hemisphere{ih}]),'jpg')
        close(fig(ifs).fh); drawnow
        end
    end
    close all
    matlabpool close force local
end
end
end
tic, fprintf('\n[%s] Saving results of virtual lesion... \n',mfilename)
save(fullfile(savedir,'strength_of_evidence.mat'),'SE'); toc

end % Main function

%%%%%%%%%%%%%%%%%%%%%%%
function FS_SUBJECT = matchSubject2FSSUBJ(subject)

switch subject
    case {'FP_96dirs_b2000_1p5iso'}
        FS_SUBJECT = 'pestilli_test';
        
    case {'KW_96dirs_b2000_1p5iso'}
        FS_SUBJECT = 'weiner';
        
    case {'MP_96dirs_b2000_1p5iso'}
        FS_SUBJECT = 'lmperry';
        
    case {'HT_96dirs_b2000_1p5iso'}
        FS_SUBJECT = 'takemura';
        
    case {'JW_96dirs_b2000_1p5iso'}
        FS_SUBJECT = 'winawer';
        
    case {'KK_96dirs_b2000_1p5iso'}
        FS_SUBJECT = 'knk';
        
    otherwise
        FS_SUBJECT = subject;
end
end


function [fh,sh] = makeBrainMap(fe,t1,slice,axLims,figName,saveDir)

% Make a map of the RMSE WITH and WITHOUT the fascicle:
coords  = feGet(fe,'roi coords') + 1;
xform   = feGet(fe,'xform img 2 acpc');
         
% Cross-validate RMSE
rmse = feGetRep(fe, 'vox rmse');
img  = feReplaceImageValues(nan(feGet(fe,'map size')),rmse,coords);
maxr = 50;

% Make anifti file from the rmse
ni  = niftiCreate('data',mbaNormalize(img,[0,1]), ...
    'qto_xyz',xform, ...
    'fname','rmse', ...
    'data_type',class(img));

% Open a figure
fh = mrvNewGraphWin(figName);

% Show the anatomy with the overlay
sh = mbaDisplayOverlay(t1, ni, slice, [], 'hot');

axis(axLims)

saveMap(fh,figName,saveDir,nanmean(img(:)),nanmedian(img(:)),nanstd(img(:)),maxr)
end

%---------------------------------%
function saveMap(fh,figName,saveDir,M,m,SD,maxfd)
% This helper function saves two figures for each map and eps with onlythe
% axis and a jpg with only the brain slice.
% The two can then be combined in illustrator.
%
% First we save only the slice as jpeg.
set(gca,'fontsize',16,'ztick',[-20 0 20 40], ...
    'xtick',[-50 -25 0 25 50], ...
    'tickdir','out','ticklength',[0.025 0])
axis off
saveFig(fh,fullfile(saveDir,'maps',figName),'tiff')
saveFig(fh,fullfile(saveDir,'maps',figName),'png')

% Then we save the slice with the axis as
% eps. This will only generate the axis
% that can be then combined in illustrator.
axis on
grid off

title(sprintf('mean %2.2f | median %2.2f | SD %2.2f', ...
    M,m,SD),'fontsize',16)
zlabel('Z (mm)','fontsize',16)
xlabel('X (mm)','fontsize',16)
cmap = colormap(hot(255));
colorbar('ytick',linspace(0,1,5),'yticklabel', ...    
    {linspace(0,1,5)*50}, ...
    'tickdir','out','ticklength',[0.025 0],'fontsize',16)
saveFig(fh,fullfile(saveDir,'maps',figName),1)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveFig(h,figName,type)

% MAke sure the folder to save the figure exists
[p,f,e] = fileparts(figName);
[success,message] = mkdir(p);
if ~isempty(message), disp(sprintf('%s.',message));end

% Find out which type of figure and geenerate the proper printing command.
switch type
    case {0,'jpeg','jpg'}
        printCommand = (sprintf('print(%s, ''-djpeg90'',''-r500'' , ''-noui'', ''-opengl'', ''%s'')', num2str(h),figName));
    case {1,'eps'}
        printCommand = (sprintf('print(%s, ''-cmyk'', ''-depsc2'',''-tiff'',''-r500'' , ''-noui'', ''%s'')', num2str(h),figName));
    case 'png'
        printCommand =  (sprintf('print(%s, ''-dpng'',''-r500'', ''%s'')', num2str(h),figName));
    case 'tiff'
        printCommand = (sprintf('print(%s, ''-dtiff'',''-r500'', ''%s'')', num2str(h),figName));
    case 'bmp'
        printCommand = (sprintf('print(%s, ''-dbmp256'',''-r500'', ''%s'')', num2str(h),figName));
    otherwise
        keyboard
end

% do the printing here:
fprintf('[%s] saving figure... \n%s\n',mfilename,figName);
eval(printCommand);
end