function s_fe_mt_ips_tract(hemisphere,saveDir)
%
% This script performs a test of conenctivity of MT+ (LO1 and LO2) with
% IPS0. THe following are the steps we perform:
%  - It loads a culled half hemisphere connectome (FE structure). 
%  - It loads the MT+ ROI.
%  - It loads the IPS0 ROI. 
%  - It finds the fascicles connecting MT+ and IPS0
%  - It reduces the connectome to the voxels and fibers of the conenctions
%  - It Perform a bootstrap test WITH and WITHOUT the connection between
%    MT+ and IPS0.
%
% Written by Franco Pestilli (c) Stanford University, Vista Team 2013

% Handle parallel computing
if matlabpool('size') == 0
    c = parcluster;
    c.NumWorkers = 12;
    matlabpool(c);
end
% Get the base directory for the data
datapath = '/marcovaldo/frk/2t1/predator/';
subjects = {...   
    'KK_96dirs_b2000_1p5iso', ...
    'JW_96dirs_b2000_1p5iso', ... 
    'HT_96dirs_b2000_1p5iso', ...
    'KW_96dirs_b2000_1p5iso', ...
    'FP_96dirs_b2000_1p5iso', ...
    'MP_96dirs_b2000_1p5iso', ...
    };

if notDefined('saveDir'), savedir = fullfile('/marcovaldo/frk/Dropbox','pestilli_etal_revision',mfilename);end
if notDefined('trackingType'), trackingType = 'lmax10';end
if notDefined('hemisphere'), hemisphere = {'left','right'};end
if notDefined('plotAnatomy'), plotAnatomy = 0;end

% Logical operation to perform ith each ROI
roi_operations   = {'and','and'};

for iSbj = 1
% Load the
t1File = fullfile(datapath,subjects{isbj},'/t1/t1.nii.gz');

% Load the FE structure
saveDir         = fullfile(savedir,subjects{isbj});
connectomesPath = fullfile(datapath,subjects{isbj},'connectomes');
feFileToLoad    = dir(fullfile(connectomesPath,sprintf('*%s*.mat',trackingType)));
fname           = feFileToLoad(probIndex).name(1:end-4);
feFileToLoad    = fullfile(connectomesPath,fname);
fprintf('[%s] Loading: \n%s\n ======================================== \n\n',mfilename,feFileToLoad)
load(feFileToLoad);

for ih = 1:length(hemisphere)
    % Set all the variables that depend on the hemisphere
    switch hemisphere{ih}
        case {'left'}
            ipsRoi     = 'LIPS0_cleaned.mat';
            mtRoi      = 'LTO1_O2_cleaned.mat';
            axisLims   = [-67 -18 -110 -40 -18 80];
            vw         = [-75,30];
            slices     = {[-18 0 0],[0 -40 0],[0 0 -14 ]};
            lght       = 'left';
            SLaxLims   = [-55 2 -120 120 -20 40 ];
            histcolor{1}  = [0.4 0.4 0.4];
            histcolor{2}  = [.6 0.4 0.4];
            
        case {'right'}
            ipsRoi     = 'RIPS0_cleaned.mat';
            mtRoi      = 'RTO1_O2_cleaned.mat';
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
    roiDir       = fullfile(datapath,subjects{isbj},'rois');
    mtFileName   = fullfile(roiDir,mtRoi);
    ips0FileName = fullfile(roiDir,ipsRoi);
    
    % Find the fascicles in the connectome that touch both MT+ and IPS0.
    mt   = dtiReadRoi(mtFileName);
    ips0 = dtiReadRoi(ips0FileName);
    rois = {mt,ips0};
    fg   = feGet(fe,'fg acpc');
    [mtIpsFG, keepFascicles] = feSegmentFascicleFromConnectome(fg, rois, roi_operations, 'mt_ips_zero');
    
    % Clean the fibers by length, fibers that too long are likely to go far
    % frontal and not just touch MT+ and IPS0.
    [Lnorm, Lmm]   = mbaComputeFiberLengthDistribution(mtIpsFG, false);
    maxSD          = 1; % Max standard deviation of the fibers to keep in the group.
    fibers2delete  = Lnorm > maxSD;
    
    % Now let's get the indices of the fibers in the FE structure:
    fasIndices = find(keepFascicles);
    fasIndices = fasIndices(fibers2delete);
    
    % Now let's mark the fascicles as deleted.
    keepFascicles(fasIndices) = false;
    
    % Perform a virtual lesion: MT+ and IPS0.
    [S, fh] = feVirtualLesion(fe,keepFascicles,0);
    figName = sprintf('virtual_lesion_distributions_%s',trackingType);
    saveFig(fh,fullfile(saveDir,figName),'eps')
    
    if plotAnatomy
        slice  = {[0 -56 0],[0 -58 0],[0 -60 0],[0 -62 0],[0 -64 0], ...
                  [0 -66 0],[0 -68 0],[0 -70 0],[0 -72 0],[0 -74 0]};
              
        % Load the T1 file for display
        t1     = niftiRead(t1File);
        
        for iS = 1:length(slice)
            % Make a figure of brain slice of the RMSE with the fascicle
            figName = sprintf('rmse_map_UNLESIONED_%s_slice%i_%s',feFileName(1:end-4),slice{iS}(find(slice{iS})),hemisphere{ih});
            shW = makeBrainMap(feWithFas,t1,slice{iS},SLaxLims,figName,saveDir);
            
            % Make a figure of brain slice of the RMSE without the fascicle
            figName = sprintf('rmse_map_LESIONED_%s_slice%i_%s',feFileName(1:end-4),slice{iS}(find(slice{iS})),hemisphere{ih});
            shWO = makeBrainMap(feWithoutFas,t1,slice{iS},SLaxLims,figName,saveDir);
        end
        
        % Show te new fiber group
        mtIpsFG = fgExtract(mtIpsFG,Lnorm < maxSD,'keep');
        figName = sprintf('MT_IPS0_connection_anatomy_%s_%s',hemisphere{ih},feFileName(1:end-4));
        figureHandle = mrvNewGraphWin(figName);
        h  = mbaDisplayBrainSlice(t1, slices{1});
        hold on
        h  = mbaDisplayBrainSlice(t1, slices{2});
        h  = mbaDisplayBrainSlice(t1, slices{3});
        tractColor = [.8 .6 .2];%[.3 .7 .9];
        [figureHandle, lightHandle, sHandle] = mbaDisplayConnectome(mtIpsFG.fibers,figureHandle,tractColor,'uniform');
        delete(lightHandle)
        view(vw(1),vw(2));
        axis(axisLims);
        lightHandle = camlight(lght);
        set(gcf,'Position',[0.0148 0.0148 .35 .87])
        drawnow
        saveFig(figureHandle,fullfile(saveDir,figName),'jpg')
    end
end
end

end % Main function

%%%%%%%%%%%%%%%%%%%%%%%
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