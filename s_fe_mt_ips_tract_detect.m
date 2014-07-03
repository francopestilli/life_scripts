function s_fe_mt_ips_tract_detect(hemisphere)
%
% This script performs a test of conenctivity of MT+ (Zilles t al ROI from 
% freesurfer) with the Superior parietal Cortex (Aparc Freesurfer segementation). 
% The following are the steps we perform:
%  - It loads a whole-brain tractography solution. 
%  - It loads the MTand Parietal ROI. 
%  - It finds the tract connecting the two ROIs.
%  - It builds a connectome model in the ROI defined by the tract
%  - It performs a series of virtual lesions on fibers with different weigths.
%    we compute the quartiles of the fibers by weight and leasion different number of them by selecting them within each quartile.
%    This test is to show .
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
    'MP_96dirs_b2000_1p5iso', ...
    'KK_96dirs_b2000_1p5iso', ...
    'JW_96dirs_b2000_1p5iso', ...
    'HT_96dirs_b2000_1p5iso', ...
    'KW_96dirs_b2000_1p5iso', ...
    };

if notDefined('saveDir'), savedir = fullfile('/marcovaldo/frk/Dropbox','pestilli_etal_revision',mfilename);end
if notDefined('trackingType'), trackingType = 'lmax10';end
if notDefined('hemisphere'), hemisphere = {'left','right'};end
if notDefined('plotAnatomy'), plotAnatomy = true;end
anatomyPath     = '/marcovaldo/frk/2t1/anatomy/';
            
for iSbj = 1:length(subjects)
% Load the FE structure
saveDir         = fullfile(savedir,subjects{iSbj});
fibergroupPath  = fullfile(datapath,subjects{iSbj},'fibers');
fgFileToLoad    = dir(fullfile(fibergroupPath,sprintf('*%s*.pdb',trackingType)));
fname           = fgFileToLoad(1).name;
fgFileToLoad    = fullfile(fibergroupPath,fname);
fgFileToLoad1   = dir(fullfile(fibergroupPath,sprintf('*%s*.pdb','lmax8')));
fname1          = fgFileToLoad1(1).name;
fgFileToLoad1   = fullfile(fibergroupPath,fname1);
fgFileToLoad2   = dir(fullfile(fibergroupPath,sprintf('*%s*.pdb','lmax2')));
fname2          = fgFileToLoad2(1).name;
fgFileToLoad2   = fullfile(fibergroupPath,fname2);

for ih = 1:length(hemisphere)
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
            histcolor{1}  = [0.4 0.4 0.4];
            histcolor{2}  = [.6 0.4 0.4];
            
        case {'right'}
            parietalRoiName = 'rh_superiorparietal_label_smooth3mm.nii.gz';
            mtRoiName       = 'rh_MT_label_smooth3mm.nii.gz';
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
    [~, keep]        = mbaComputeFibersOutliers(mt2parietalTract,3,3);
    fprintf('\n[%s] Found a tract with %i fibers... \n',mfilename,sum(keep))
    mt2parietalTract = fgExtract(mt2parietalTract,find(keep),'keep');
    toc
    
    % Find the Coordinates of the mt-parietal tract
    tic, fprintf('\n[%s] Create ROI from MT-Parietal tract... \n',mfilename)
    tractRoi = dtiCreateRoiFromFibers(mt2parietalTract);
    %tractRoi = dtiRoiClean(tractRoi, 12, 'dilate');toc % We smooth the ROI to enlarge the pathneighborhod fibers for visualization
    
    tic, fprintf('\n[%s] Removing fibers not going throught the tractROI... \n',mfilename)
    [fg,~, ~, ~] = dtiIntersectFibersWithRoi([],'and',2,tractRoi,fg);
    fg = feClipFibersToVolume(fg,tractRoi.coords,1);toc
    
    % Build LiFE model only in this volume, fit, cull
    dwiPath       = fullfile(datapath,subjects{iSbj},'diffusion_data');
    dwiFiles      = dir(fullfile(dwiPath,sprintf('run*.gz')));
    dwiFile       = fullfile(dwiPath,dwiFiles(1).name);
    dwiFileRepeat = fullfile(dwiPath,dwiFiles(2).name);
    t1File        = fullfile(datapath,subjects{iSbj},'anatomy','t1.nii.gz');
    
    % Directory where to save the fe structures
    saveDirC   = fullfile(datapath,subjects{iSbj},'connectomes');
    feFileName = [parietal.name, '_', fname(1:40), '.mat'];
    fe   = feConnectomeInit(dwiFile,fg,feFileName,saveDirC,dwiFileRepeat,t1File);
    cull(iSbj,ih).name = [subjects{iSbj},'_', hemisphere{ih}];
    cull(iSbj,ih).iter = 1;
    cull(iSbj,ih).minwval = eps;
    M    = feGet(fe,'mfiber');
    dSig = feGet(fe,'dsigdemeaned');
    fe   = feSet(fe,'fit',feFitModel(M,dSig,'bbnnls'));
    sw   = sum(feGet(fe,'fiber weights') <= cull(iSbj,ih).minwval);
    cull(iSbj,ih).num2delete(cull(iSbj,ih).iter) = sw;
    cull(iSbj,ih).numtotal(cull(iSbj,ih).iter)   = size(M,2);
    cull(iSbj,ih).rmse(cull(iSbj,ih).iter)       = mean(feGet(fe,'vox rmse'));
    cull(iSbj,ih).rmse(cull(iSbj,ih).iter)       = mean(feGetRep(fe,'vox rmse'));
    cull(iSbj,ih).rrmse(cull(iSbj,ih).iter)      = mean(feGetRep(fe,'vox rmse ratio'));
    
    % Cull the connectome
    while (sw ~= 0)
        fprintf('Number of total %i and zero-weight %i fibers | culling...\n', ...
            cull(iSbj,ih).numtotal(cull(iSbj,ih).iter),cull(iSbj,ih).num2delete(cull(iSbj,ih).iter))
        M    = feGet(fe,'mfiber');
        dSig = feGet(fe,'dsigdemeaned');
        fe   = feSet(fe,'fit',feFitModel(M,dSig,'bbnnls'));
        sw   = sum(feGet(fe,'fiber weights') <= cull(iSbj,ih).minwval);
        fe   = feConnectomeReduceFibers(fe, find((feGet(fe,'fiber weights') > cull(iSbj,ih).minwval)));
        cull(iSbj,ih).iter = cull(iSbj,ih).iter + 1;
        cull(iSbj,ih).num2delete(cull(iSbj,ih).iter) = sw;
        cull(iSbj,ih).numtotal(cull(iSbj,ih).iter)   = size(M,2);
        cull(iSbj,ih).rmse(cull(iSbj,ih).iter)       = mean(feGet(fe,'vox rmse'));
        cull(iSbj,ih).rmse(cull(iSbj,ih).iter)       = mean(feGetRep(fe,'vox rmse'));
        cull(iSbj,ih).rrmse(cull(iSbj,ih).iter)      = mean(feGetRep(fe,'vox rmse ratio'));
    end
    w(iSbj,ih).name = [subjects{iSbj},'_', hemisphere{ih}];
    w(iSbj,ih).lw = log10(feGet(fe,'fiber weights'));
    w(iSbj,ih).quartiles = [.2 .8];
    w(iSbj,ih).lwvalues  = quantile(w(iSbj,ih).lw, w(iSbj,ih).quartiles);
    
    % Make aplot of the distribution of fiber weights
    fig(1).name = sprintf(['log_weights_distribution', '_',hemisphere{ih}]);
    fig(1).fh   = figure('name',fig.name,'color','w');
    w(iSbj,ih).bins = linspace(-6,0,60);
    [w(iSbj,ih).y,w(iSbj,ih).x] = hist(w(iSbj,ih).lw,w(iSbj,ih).bins);
    w(iSbj,ih).y = w(iSbj,ih).y ./ sum(w(iSbj,ih).y);
    patch([w(iSbj,ih).x,fliplr(w(iSbj,ih).x)],[w(iSbj,ih).y,zeros(size(w(iSbj,ih).y))],'k')
    set(gca,'box','off','ylim',[0 .2],'ytick',[0 .1 .2], ...
        'xlim',[-6 0],'xtick',[-6 -3 0],'fontsize',16, 'tickdir', 'out', ...
        'ticklength', [0.025 0])
    saveFig(fig(1).fh,fullfile(saveDir,fig(1).name),'eps')
    
    % Perform a virtual lesion: MT+ and parietal.
    fibersTest.num2delete = [1,4,16,64];
    for iqr = 1:length(w(iSbj,ih).quartiles)
        for ifib = 1:length(fibersTest.num2delete)
            % Now decide which fascicles to keep depending on their position in the
            % quantiles and on the number o fibers that we want to delete
            if     (iqr == 1)
                useTheseFibers = find(w(iSbj,ih).lw <= w(iSbj,ih).lwvalues(iqr));
            elseif (iqr == 2)
                useTheseFibers = find(w(iSbj,ih).lw >= w(iSbj,ih).lwvalues(iqr));
            end
            
            % Keep the length of the fibers to a resonable level
            tmpfg = feGet(fe,'fibers acpc');
%             len   = 0; c = 0;
%             while ((len < 6) || (len > 30))
%                 c = c + 1;
                keepFascicles = randsample(useTheseFibers,fibersTest.num2delete(ifib));
%                len = mean(cell2mat(cellfun(@length,tmpfg.fibers(keepFascicles),'UniformOutput',false)));
%                 if c > 30000; keyboard; end
%             end
            SE(iSbj,ih,iqr, ifib) = feVirtualLesion(fe,keepFascicles);
        end
    end
    matlabpool close force local
    tic, fprintf('\n[%s] Saving results of virtual lesion... \n',mfilename)
    save(fullfile(savedir,'strength_of_evidence_quartiles.mat'),'SE','w','cull','fibersTest'); toc
end
end

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
        printCommand = (sprintf('print(%s, ''-djpeg80'',''-r300'' , ''-noui'', ''-opengl'', ''%s'')', num2str(h),figName));
    case {1,'eps'}
        printCommand = (sprintf('print(%s, ''-cmyk'', ''-depsc2'',''-tiff'',''-r300'' , ''-noui'', ''%s'')', num2str(h),figName));
    case 'png'
        printCommand =  (sprintf('print(%s, ''-dpng'',''-r300'', ''%s'')', num2str(h),figName));
    case 'tiff'
        printCommand = (sprintf('print(%s, ''-dtiff'',''-r300'', ''%s'')', num2str(h),figName));
    case 'bmp'
        printCommand = (sprintf('print(%s, ''-dbmp256'',''-r300'', ''%s'')', num2str(h),figName));
    otherwise
        keyboard
end

% do the printing here:
fprintf('[%s] saving figure... \n%s\n',mfilename,figName);
eval(printCommand);
end