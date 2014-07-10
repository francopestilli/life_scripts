function s_fe_make_mt_parietal_tract(hemisphere,iSbj)
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

% Handle parallel computing
if matlabpool('size') == 0
    c = parcluster;
    c.NumWorkers = 12;
    matlabpool(c);
end
% Get the base directory for the data
datapath = '/marcovaldo/frk/2t1/predator/';
subjects = {...   
    'FP_96dirs_b2000_1p5iso', ...
    'KK_96dirs_b2000_1p5iso', ...
    'JW_96dirs_b2000_1p5iso', ... 
    'HT_96dirs_b2000_1p5iso', ...
    'KW_96dirs_b2000_1p5iso', ...
    'MP_96dirs_b2000_1p5iso', ...
    };

if notDefined('saveDir'), savedir = fullfile('/marcovaldo/frk/Dropbox','pestilli_etal_revision',mfilename);end
if notDefined('trackingType'), trackingType = 'lmax10';end
if notDefined('hemisphere'), hemisphere = {'left','right'};end
if notDefined('plotAnatomy'), plotAnatomy = 1;end
anatomyPath     = '/dev/data/anatomy/';
parietalRoiName = 'lh_MT_label_smooth3mm.nii.gz';
mtRoiName       = 'rh_MT_label_smooth3mm.nii.gz';
parietalRoiName = 'lh_superiorparietal_label_smooth3mm.nii.gz';
mtRoiName       = 'rh_superiorparietal_label_smooth3mm.nii.gz';
            
% Logical operation to perform ith each ROI
roi_operations   = {'and','and'};

%for iSbj = 1:length(subjects)
% Load the
t1File = fullfile(datapath,subjects{iSbj},'/anatomy/t1.nii.gz');

% Load the FE structure
saveDir        = fullfile(savedir,subjects{iSbj});
fibergroupPath = fullfile(datapath,subjects{iSbj},'fibers');
fgFileToLoad   = dir(fullfile(fibergroupPath,sprintf('*%s*.pdb',trackingType)));
fname          = fgFileToLoad(1).name;
fgFileToLoad   = fullfile(fibergroupPath,fname);
fprintf('[%s] Loading: \n%s\n ======================================== \n\n',mfilename,fgFileToLoad)
fg = fgRead(fgFileToLoad);

for ih = 1:length(hemisphere)
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
    roi      = dtiNewRoi('combined_mt_parietal');
    roi.coords = [mt.coords;parietal.coords];
    tic, fprintf('\n[%s] Segmenting tract from connectome... \n',mfilename)
    [mt2parietalTract, keepFascicles] = feSegmentFascicleFromConnectome(fg, {roi}, {'and both endpoints'}, 'mt_parietal');
    fprintf('\n[%s] DONE segmenting tract from connectome in %2.3f \n',mfilename,toc)

    % Clean the fibers by length, fibers that too long are likely to go far
    % frontal and not just touch MT+ and parietal.
    [Lnorm, Lmm]   = mbaComputeFiberLengthDistribution(mt2parietalTract, false);
    maxSD          = 1; % Max standard deviation of the fibers to keep in the group.
    fibers2delete  = Lnorm > maxSD;
    
    % Now let's get the indices of the fibers in the FE structure:
    fasIndices = find(keepFascicles);
    fasIndices = fasIndices(~fibers2delete);
     
    % Find the Coordinates of the mt-parietal tract
    fgTract     = fgExtract(mt2parietalTract,fasIndices,'keep');
    
    if plotAnatomy
        slice  = {[0 -56 0],[0 -58 0],[0 -60 0],[0 -62 0],[0 -64 0], ...
                  [0 -66 0],[0 -68 0],[0 -70 0],[0 -72 0],[0 -74 0]};
              
        % Load the T1 file for display
        t1     = niftiRead(t1File);

        % Show te new fiber group
        figName = sprintf('MT_IPS0_connection_anatomy_%s',hemisphere{ih});
        figureHandle = mrvNewGraphWin(figName);
        h  = mbaDisplayBrainSlice(t1, slices{1});
        hold on
        h  = mbaDisplayBrainSlice(t1, slices{2});
        h  = mbaDisplayBrainSlice(t1, slices{3});
        tractColor = [.8 .6 .2];%[.3 .7 .9];
        [figureHandle, lightHandle, sHandle] = mbaDisplayConnectome(fgTract.fibers,figureHandle,tractColor,'uniform');
        delete(lightHandle)
        view(vw(1),vw(2));
        axis(axisLims);
        lightHandle = camlight(lght);
        set(gcf,'Position',[0.0148 0.0148 .35 .87])
        drawnow
        keyboard
        saveFig(figureHandle,fullfile(saveDir,figName),'jpg')
        
    end
end
save(fullfile(saveDir,'strength_of_evidence.mat'),'S','DKL')
%end

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
        keyboard
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