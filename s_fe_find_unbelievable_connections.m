function s_fe_find_unbelievable_connections()
%
% This script performs a test of conenctivity between diffrerent cortical 
% ROIs that are not believable. 
%
% Copyright by Franco Pestilli Stanford University, 2014

% Handle parallel computing
if matlabpool('size') == 0
    c = parcluster;
    c.NumWorkers = 12;
    matlabpool(c);
end
% Get the base directory for the data
subjects = {...
    'FP_96dirs_b2000_1p5iso', ...
    'HT_96dirs_b2000_1p5iso', ...
    'KW_96dirs_b2000_1p5iso', ...
    'MP_96dirs_b2000_1p5iso', ...
    'KK_96dirs_b2000_1p5iso', ...
    'JW_96dirs_b2000_1p5iso', ...
    };

anatomyPath     = '/marcovaldo/frk/2t1/anatomy/';
fibersPath      = '/marcovaldo/frk/2t1/predator/';
fibs            = {'run01_fliprot_aligned_trilin_csd_lmax10_run01_fliprot_aligned_trilin_brainmask_run01_fliprot_aligned_trilin_wm_prob-500000_recomputed-rejected.mat', ...
                   'run01_fliprot_aligned_trilin_csd_lmax10_run01_fliprot_aligned_trilin_brainmask_run01_fliprot_aligned_trilin_wm_prob-500000_recomputed-optimized.mat'};
        
for iSbj = 1:length(subjects)
    % Find all the ROI computed for each subject:
    roiDir     = fullfile(anatomyPath,matchSubject2FSSUBJ(subjects{iSbj}),'label');
    allRois = dir(fullfile(roiDir,'*.mat'));
    % Compute all the combinations of the ROIs
    roisIndicesToTest = combntns(1:length(allRois),2);
    
    for ifb = 1:length(fibs)
        tic, fprintf('\n[%s] Loading the fiber group... \n',mfilename)
        fg = fgRead(fullfile(sprintf('%s',fibersPath),subjects{iSbj},'fibers',fibs{ifb}));toc
        parfor ir = 1:10%size(roisIndicesToTest,1)
            roi1 = dtiReadRoi(fullfile(roiDir,allRois(roisIndicesToTest(ir,1)).name));
            roi2 = dtiReadRoi(fullfile(roiDir,allRois(roisIndicesToTest(ir,2)).name));
      
            % Combine the two ROIs into a single ROI
            roi(ir) = roi1;
            roi(ir).coords = [roi1.coords; roi2.coords];
            roi(ir).name   = [roi1.name, 'to', roi2.name];
        end
        [fgOut{iSbj,ifb},contentiousFibers{iSbj,ifb}, keep{iSbj,ifb}, keepID{iSbj,ifb}] = dtiIntersectFibersWithRoi([],['divide','bothendpoints'], [], roi, fg);

%         keyboard
%         
%         [tract{iSbj,ifb,ir}, keepFascicles{iSbj,ifb,ir}] = feSegmentFascicleFromConnectome(fg, {roi1,roi2}, {'endpoints','endpoints'}, 'tmp');toc
%         if ~(sum(keepFascicles{iSbj,ifb,ir})==0)
%             [~, keep]        = mbaComputeFibersOutliers(tract{iSbj,ifb,ir},2,2);
%             fprintf('\n[%s] Found a tract with %i fibers touching %s and %s... ... \n',mfilename,sum(keep),allRois(roisIndicesToTest(ir,1)).name,allRois(roisIndicesToTest(ir,2)).name)
%             if ~isempty(find(keep, 1))
%                 tract{iSbj,ifb,ir} = fgExtract(tract{iSbj,ifb,ir},find(keep,1),'keep');
%                 keepFascicles{iSbj,ifb,ir}(~keep) = 0;
%             end
%         end
%           roiNames{iSbj,ifb,ir} = {allRois(roisIndicesToTest(ir,1)).name,allRois(roisIndicesToTest(ir,2)).name};
%         end
    end
end

keyboard
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

% 
% function [fh,sh] = makeBrainMap(fe,t1,slice,axLims,figName,saveDir)
% 
% % Make a map of the RMSE WITH and WITHOUT the fascicle:
% coords  = feGet(fe,'roi coords') + 1;
% xform   = feGet(fe,'xform img 2 acpc');
%          
% % Cross-validate RMSE
% rmse = feGetRep(fe, 'vox rmse');
% img  = feReplaceImageValues(nan(feGet(fe,'map size')),rmse,coords);
% maxr = 50;
% 
% % Make anifti file from the rmse
% ni  = niftiCreate('data',mbaNormalize(img,[0,1]), ...
%     'qto_xyz',xform, ...
%     'fname','rmse', ...
%     'data_type',class(img));
% 
% % Open a figure
% fh = mrvNewGraphWin(figName);
% 
% % Show the anatomy with the overlay
% sh = mbaDisplayOverlay(t1, ni, slice, [], 'hot');
% 
% axis(axLims)
% 
% saveMap(fh,figName,saveDir,nanmean(img(:)),nanmedian(img(:)),nanstd(img(:)),maxr)
% end
% 
% %---------------------------------%
% function saveMap(fh,figName,saveDir,M,m,SD,maxfd)
% % This helper function saves two figures for each map and eps with onlythe
% % axis and a jpg with only the brain slice.
% % The two can then be combined in illustrator.
% %
% % First we save only the slice as jpeg.
% set(gca,'fontsize',16,'ztick',[-20 0 20 40], ...
%     'xtick',[-50 -25 0 25 50], ...
%     'tickdir','out','ticklength',[0.025 0])
% axis off
% saveFig(fh,fullfile(saveDir,'maps',figName),'tiff')
% saveFig(fh,fullfile(saveDir,'maps',figName),'png')
% 
% % Then we save the slice with the axis as
% % eps. This will only generate the axis
% % that can be then combined in illustrator.
% axis on
% grid off
% 
% title(sprintf('mean %2.2f | median %2.2f | SD %2.2f', ...
%     M,m,SD),'fontsize',16)
% zlabel('Z (mm)','fontsize',16)
% xlabel('X (mm)','fontsize',16)
% cmap = colormap(hot(255));
% colorbar('ytick',linspace(0,1,5),'yticklabel', ...    
%     {linspace(0,1,5)*50}, ...
%     'tickdir','out','ticklength',[0.025 0],'fontsize',16)
% saveFig(fh,fullfile(saveDir,'maps',figName),1)
% end
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function saveFig(h,figName,type)
% 
% % MAke sure the folder to save the figure exists
% [p,f,e] = fileparts(figName);
% [success,message] = mkdir(p);
% if ~isempty(message), disp(sprintf('%s.',message));end
% 
% % Find out which type of figure and geenerate the proper printing command.
% switch type
%     case {0,'jpeg','jpg'}
%         printCommand = (sprintf('print(%s, ''-djpeg90'',''-r500'' , ''-noui'', ''-opengl'', ''%s'')', num2str(h),figName));
%     case {1,'eps'}
%         printCommand = (sprintf('print(%s, ''-cmyk'', ''-depsc2'',''-tiff'',''-r500'' , ''-noui'', ''%s'')', num2str(h),figName));
%     case 'png'
%         printCommand =  (sprintf('print(%s, ''-dpng'',''-r500'', ''%s'')', num2str(h),figName));
%     case 'tiff'
%         printCommand = (sprintf('print(%s, ''-dtiff'',''-r500'', ''%s'')', num2str(h),figName));
%     case 'bmp'
%         printCommand = (sprintf('print(%s, ''-dbmp256'',''-r500'', ''%s'')', num2str(h),figName));
%     otherwise
%         keyboard
% end
% 
% % do the printing here:
% fprintf('[%s] saving figure... \n%s\n',mfilename,figName);
% eval(printCommand);
% end