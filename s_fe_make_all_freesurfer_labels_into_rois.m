function s_fe_make_all_freesurfer_labels_into_rois(hemisphere)
%
% This script creates all the freesurfer lables from parcellation files and 
% creates corresponding nifti and mat ROI compatible with mrDiffusion.
%
% Copyright by Franco Pestilli Stanford University 2014

% Get the base directory for the data
subjects = {...
    'KK_96dirs_b2000_1p5iso' ...
    'JW_96dirs_b2000_1p5iso', ...
    'KW_96dirs_b2000_1p5iso', ...
    'MP_96dirs_b2000_1p5iso', ...
    'HT_96dirs_b2000_1p5iso', ...
    'FP_96dirs_b2000_1p5iso', ...
    };
if notDefined('annotationFileName')
    annotationFileName = {'aparc','aparc.a2009s'};
end
fsSubjectsDir = getenv('SUBJECTS_DIR');
clobber = false;

for iSbj = 1:length(subjects)
    if ~isdir(fullfile(fsSubjectsDir,subjects{iSbj}))
         fs_subject = matchSubject2FSSUBJ(subjects{iSbj});
    else fs_subject = subjects{iSbj};
    end
    fsSubjectDir   = fullfile(fsSubjectsDir,fs_subject);

    % Create all the necessary label files
    for ia = 1:length(annotationFileName)
        fs_annotationToLabelFiles(fs_subject,annotationFileName{ia},[],fsSubjectsDir);
    end
    
    % File all the label ROIs for this subject
    labelFileNames   = dir(fullfile(fsSubjectDir,'label','*.label'));
    labelRoiName     = cell(length(labelFileNames),1);
    niftiRoiFullPath = cell(length(labelFileNames),1);
    matRoiFullPath  = cell(length(labelFileNames),1);
    for il = 1:length(labelFileNames)
        labelRoiName{il}  = labelFileNames(il).name;
        niftiRoiName      = labelRoiName{il};
        niftiRoiName(niftiRoiName=='.') = '_';
        niftiRoiFullPath{il}  = fullfile(fsSubjectDir,'label',  niftiRoiName);
        matRoiFullPath{il}   = [fullfile(fsSubjectDir,'label',  niftiRoiName),'_smooth3mm_ROI.mat'];
    end
    
    for il = 1:length(labelFileNames)
        if ~(exist([niftiRoiFullPath{il},'_smooth3mm.nii.gz'],'file')==2) || clobber
            fs_labelFileToNiftiRoi(fs_subject,labelRoiName{il},niftiRoiFullPath{il},labelFileNames(il).name(1:2),[],[],fsSubjectsDir);
        else
            fprintf('[%s] Found ROI, skipping: \n%s\n',mfilename,niftiRoiFullPath{il})
        end
        if ~(exist([matRoiFullPath{il}],'file')==2) || clobber
            dtiImportRoiFromNifti([niftiRoiFullPath{il},'_smooth3mm.nii.gz'], matRoiFullPath{il});
        else
            fprintf('[%s] Found ROI, skipping: \n%s\n',mfilename,niftiRoiFullPath{il})
        end
    end
end

end % Main function

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