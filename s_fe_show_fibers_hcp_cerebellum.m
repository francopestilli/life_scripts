function fe = s_fe_show_fibers_hcp_cerebellum(bval,dataDir,which_subject)
%
% This function:
%  - Initializes a LIFE structure from a candidate connectome
%  - Generates an optimized connectome from a cadidate connectome using
%  LIFE method
%
%  fe = s_fe_fit_hcp()
%
% INPUTS:  none
% OUTPUTS: fe structure the optimized life structure
%
% Copyright Franco Pestilli (2014) Vistasoft Stanford University
if notDefined('bval'); bval=2000;end
if notDefined('dataDir'); dataDir='2t1';end
figVisible = 'off';

% Get the base directory for the data
[~,hostname] = system('hostname');
hostname = deblank(hostname);
switch dataDir
    case {'2t2'}
        switch hostname
            case {'marcovaldo'}
                datapath = '/home/frk/2t2/HCP/';
            otherwise
                datapath = '/marcovaldo/frk/2t2/HCP/';
                
        end
        subjects = {...
            '115320', ...
            '117122', ...
            '118730', ...
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

for isbj = which_subject
    feOpenLocalCluster;
    
    % Build the file names for the diffusion data, the anatomical MR, the fiber
    % group containing the connectome and the
    t1File        = fullfile(datapath,subjects{isbj},'anatomy','T1w_acpc_dc_restore_1p25.nii.gz');
    %cerebellar_wm = fullfile(anatomypath,subjects{isbj},'cerebellar_wm.nii.gz');
    savedir       = fullfile(datapath,subjects{isbj},'figures');
    
    % Now find all the fiber files that we will analyze
    fibersPath = fullfile(datapath,subjects{isbj},'fibers');
    fgFiles    = dir(fullfile(fibersPath,sprintf('*%s*cerebellum-optimized.mat',num2str(bval))));
    fgFileName = fullfile(fibersPath,fgFiles(1).name);
    fg         = fgRead(fgFileName);
    %  cerebellar_wm = dtiRoiFromNifti(cerebellar_wm,1,[],'mat',[],false);
    % Show the ROI:
    % plot3(cerebellar_wm.coords(:,1),cerebellar_wm.coords(:,2),cerebellar_wm.coords(:,3),'ro');
    % view(-40,40); axis equal; hold on
    
    t1     = niftiRead(t1File);
    
    % Show te new fiber group
    figName = sprintf('Cerebellum_fibers_optimized');
    figureHandle = figure('name',figName,'visible',figVisible);
    h  = mbaDisplayBrainSlice(t1, [1 0 0]);
    hold on
    tractColor = [.2  .4 .95]; subset = (randsample(1:length(fg.fibers),800));
    [figureHandle, lightHandle] = mbaDisplayConnectome(mbaFiberSplitLoops(fg.fibers(subset)),figureHandle,tractColor,'uniform');
    delete(lightHandle);
    view(90,0);
    axis([-65 65 -75 5 -62 10])
    lightHandle = camlight('right');
    drawnow
    saveFig(figureHandle,fullfile(savedir,figName),'jpeg')
    
    % Now remove the Spinal Tract
    roi = dtiCreateRoiFromFibers(fg);
    % Show the ROI:
    % plot3(roi.coords(:,1),roi.coords(:,2),roi.coords(:,3),'ro');view(-90,0)
    
    % Clip the ROI to the AC-PC plane
    roi = dtiRoiClip(roi, [], [], [-45 90]);
    % Show it after the trasformation:
    % hold on; plot3(roi.coords(:,1),roi.coords(:,2),roi.coords(:,3),'k*');
     
    % Intersect the connectome with the Spinal ROI
    tic, fprintf('[%s] Removing fibers not going through the Spinal ROI... ',mfilename)
    fg = feSegmentFascicleFromConnectome(fg, {roi}, {'not'}, 'Not_spinal');
    % Show it after the trasformation:
    % hold on; plot3(roi.coords(:,1),roi.coords(:,2),roi.coords(:,3),'y*');
       
    % Show te new fiber group
    figName = sprintf('Cerebellum_fibers_optimized_missing_spinal');
    figureHandle = figure('name',figName,'visible',figVisible);
    h  = mbaDisplayBrainSlice(t1, [0 0 -35]);
    hold on
    tractColor = [.2  .94 .75]; subset = (randsample(1:length(fg.fibers),800));
    [figureHandle, lightHandle] = mbaDisplayConnectome(mbaFiberSplitLoops(fg.fibers(subset)),figureHandle,tractColor,'uniform');
    delete(lightHandle);
    axis([-65 65 -75 5 -62 10])
    lightHandle = camlight('right');
    view(0,90);drawnow
    saveFig(figureHandle,fullfile(savedir,figName),'jpeg')
end
end


%-------------------------------%
function saveFig(h,figName,eps)
if ~exist( fileparts(figName), 'dir'), mkdir(fileparts(figName));end
fprintf('[%s] saving figure... \n%s\n',mfilename,figName);

switch eps
    case {0,'jpeg','jpg'}
        eval(sprintf('print(%s, ''-djpeg90'', ''-opengl'', ''%s'')', num2str(h),[figName,'.jpg']));
    case {1,'eps'}
        eval(sprintf('print(%s, ''-cmyk'', ''-painters'',''-depsc2'',''-tiff'',''-r500'' , ''-noui'', ''%s'')', num2str(h),[figName,'.eps']));
    case 'png'
        eval(sprintf('print(%s, ''-dpng'',''-r500'', ''%s'')', num2str(h),[figName,'.png']));
    case 'tiff'
        eval(sprintf('print(%s, ''-dtiff'',''-r500'', ''%s'')', num2str(h),[figName,'.tif']));
    case 'bmp'
        eval(sprintf('print(%s, ''-dbmp256'',''-r500'', ''%s'')', num2str(h),[figName,'.bmp']));
    otherwise
end

end
