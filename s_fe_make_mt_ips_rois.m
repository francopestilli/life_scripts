function s_fe_make_mt_ips_rois(hemisphere)
%
% This script makes the MT and IPS rois used to find the tract reported in
% Pestilli et al., LIFE paper.
%
% Retinotopy folder: /biac3/wandell7/data/Retinotopy/
% 
% ROIs can be created by using these two function in sequence:
% (1) fs_annotationToLabelFiles(fs_subject,annotationFileName,hemisphere,labelsDir)
% (2) fs_labelFileToNiftiRoi(fs_subject,labelFileName,niftiRoiName,hemisphere,regMgzFile,smoothingKernel)
% 
% Copyright by Franco Pestilli Stanford University 2014

% Get the base directory for the data
subjects = {...   
    '105115', ...
    '110411', ...
    '111312', ...
    '113619', ...
    '115320', ...
    '117122', ...
    '118730',...
    };

if notDefined('hemisphere'), hemisphere = {'lh','rh'};end
annotationFileName = {'aparc','aparc.a2009s'};
labelFileName =  {'superiorparietal','MT'};
fsSubjectsDir = getenv('SUBJECTS_DIR');

for iSbj = 1:length(subjects)
    if ~isdir(fullfile(fsSubjectsDir))
         fs_subject = matchSubject2FSSUBJ(subjects{iSbj});
    else fs_subject = subjects{iSbj};
    end
    fsSubjectDir   = fullfile(fsSubjectsDir,fs_subject);

    % Create all the necessary label files
    for ia = 1:length(annotationFileName)
        fs_annotationToLabelFiles(fs_subject,annotationFileName{ia});
    end
    
    % Extract the label files for the ROI we want and make nifti ROIs
    for ih = 1:length(hemisphere)
        for il = 1:length(labelFileName)
            % Build a name for the nifti ROI and create it.
            labelRoiName  = sprintf('%s.%s.label',hemisphere{ih},labelFileName{il});
            niftiRoiName  = labelRoiName;
            niftiRoiName(niftiRoiName=='.') = '_';
            niftiRoiFullPath  = fullfile(fsSubjectDir,'label',  niftiRoiName);
            fs_labelFileToNiftiRoi(fs_subject,labelRoiName,niftiRoiFullPath,hemisphere{ih});
        end
    end
end

end % Main function

function FS_SUBJECT = matchSubject2FSSUBJ(subject)

switch subject
    case {'FP_96dirs_b2000_1p5iso'}
        FS_SUBJECT = 'pestilli_test';
    case {'105115'}
        FS_SUBJECT = '105115';

    case {'KW_96dirs_b2000_1p5iso'}
        
    case {'MP_96dirs_b2000_1p5iso'}
        
    case {'HT_96dirs_b2000_1p5iso'}
        
    case {'JW_96dirs_b2000_1p5iso'} 
        
    otherwise
        keyboard
end
end

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