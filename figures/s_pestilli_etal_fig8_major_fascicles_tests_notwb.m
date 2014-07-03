function s_pestilli_etal_fig8_major_fascicles_tests_notwb()
%
% Perform statistical tests on the strenght fo evidence for the 20 major
% white-matter pathways returned by AFQ.
%
% - load the 20 major fascicles generated with AFQ
% - load the corresponding FE structure
% - performs a virtual leasion for eahc one fo the fascicles
% - returns the statistics (strenght of evidence fo each one fo the major
%   fascicles)
%
% Copyright Franco Pestilli 2014 Stanford University
addpath(genpath('/marcovaldo/frk/git/AFQ/'))
figVisible = 'on';

% Get the base directory for the data
datapath = '/marcovaldo/frk/2t1/predator/';
subjects = {...
    'FP_96dirs_b2000_1p5iso', ...
    'JW_96dirs_b2000_1p5iso', ...
    'KK_96dirs_b2000_1p5iso', ...
    'MP_96dirs_b2000_1p5iso', ...
    'HT_96dirs_b2000_1p5iso', ...
    'KW_96dirs_b2000_1p5iso', ...
    };
if notDefined('saveDir'), savedir = fullfile('/marcovaldo/frk/Dropbox','pestilli_etal_revision',mfilename);end
if notDefined('trackingType'), trackingType = 'lmax10';end

for isbj = 1:length(subjects)
    saveDir        = fullfile(savedir,subjects{isbj});
    fibergroupPath = fullfile(datapath,subjects{isbj},'fibers');
    fgFileToLoad   = dir(fullfile(fibergroupPath,sprintf('*%s*.pdb',trackingType)));
    fname          = fgFileToLoad(1).name;
    fgFileToLoad   = fullfile(fibergroupPath,fname);
    fprintf('[%s] Loading: \n%s\n ======================================== \n\n',mfilename,fgFileToLoad)
    fg = fgRead(fgFileToLoad);
    
    % Segment the fibers using AFQ
    dtFile    = fullfile(datapath,subjects{isbj},'dtiInit','dt6.mat');
    [fg_classified,~,classification]= AFQ_SegmentFiberGroups(dtFile, fg,[],[],false);
    
    % Split the fiber groups into individual groups
    fascicles = fg2Array(fg_classified);

    % Load the fascicles
    for iFas = 1:length(fascicles)
        fprintf('[%s] Performing virtual lesion: %s...\n',mfilename, classification.names{iFas})
               
        % Find the Coordinates of the mt-parietal tract
        tic, fprintf('\n[%s] Create ROI for %s... \n',mfilename,classification.names{iFas})
        tractRoi = dtiCreateRoiFromFibers(fascicles(iFas));toc
        tic, fprintf('\n[%s] Removing fibers not going through the %s... \n',mfilename,classification.names{iFas})
        [fgClip,~, keep]  = dtiIntersectFibersWithRoi([],'and',2,tractRoi,fg);
        
        % Update the indices on the global whole-brain list
        classificationIdices = classification.index;
        classificationIdices( ~keep ) = 0;
        [fgClip, kept]   = feClipFibersToVolume(fgClip,tractRoi.coords,1);toc
        classificationIdices( ~kept ) = 0;
        fgClip.pathwayInfo = [];
        
        % Build LiFE model only in this volume, fit, cull
        dwiPath  = fullfile(datapath,subjects{isbj},'diffusion_data');
        dwiFiles = dir(fullfile(dwiPath,sprintf('run*.gz')));
        dwiFile       = fullfile(dwiPath,dwiFiles(1).name);
        dwiFileRepeat = fullfile(dwiPath,dwiFiles(2).name);
        t1File        = fullfile(datapath,subjects{isbj},'anatomy','t1.nii.gz');
        
        % Directory where to save the fe structures
        %savedir    = fullfile(datapath,subjects{isbj},'connectomes');
        fe   = feConnectomeInit(dwiFile,fgClip,[classification.names{iFas}, '.mat'],['none'],dwiFileRepeat,t1File);
        M    = feGet(fe,'mfiber');
        dSig = feGet(fe,'dsigdemeaned');
        fit  = feFitModel(M,dSig,'bbnnls');
        fe   = feSet(fe,'fit',fit);clear fit
        fgClip   = feGet(fe,'fibers acpc');
                
        % Get the fibers for the current fascicle
        keepFascicles = find(classificationIdices==iFas);
        
        % Perform a virtual lesion: MT+ and parietal.
        [S(isbj,ih), fh, KL(isbj,ih)] = feVirtualLesion(fe,keepFascicles,0);
        figName = sprintf('virtual_lesion_rmse_distributions_%s_%s',trackingType,classification.name{iFas});
        saveFig(fh(1),fullfile(saveDir,figName),'eps')
        figName = sprintf('virtual_lesion_mean_rmse_distributions_%s_%s',trackingType,classification.name{iFas});
        saveFig(fh(2),fullfile(saveDir,figName),'eps')
        close all; clear keepFascicles classificationIdices
    end
    clear fascicles clear fg_classified dtFile fe 
    save(fullfile(savedir,['average_s_customROI_',fname,'.mat']),'S','KL')
end
            

% Make a plot of the Strength of evidence for each fascicle across subjects.
figName = sprintf('StrengthEvidenceAcrossSubjects_%s',  fname);
fh  = figure('name',figName,'visible',figVisible,'color','w');
x = 1:size(S,2);
bar(x,mean(S,1),'k','EdgeColor','w'); hold on
semilogx([x; x], [mean(S,1);mean(S,1)] + ...
         [std(S,[],1)./sqrt(size(S,1));  ...
         -std(S,[],1)./sqrt(size(S,1))],'r-','linewidth',2)
ylabel('Strength of evidence (S)','FontSize',16,'FontAngle','oblique')
xlabel('Major fascicles','FontSize',16,'FontAngle','oblique')
set(gca,'fontsize',16, ...
    'xlim', [0 21], 'xtick', x,  ... 
    'ylim', [0 100],'ytick', [0 50 100], ...
    'box',  'off',  'tickdir', 'out', 'ticklength', [0.025 0])
saveFig(fh,fullfile(savedir, figName),1)

% Make a plot of the Strength of evidence for each fascicle across subjects.
figName = sprintf('StrengthEvidenceAcrossSubjectsBEST_%s',  fname);
fh  = figure('name',figName,'visible',figVisible,'color','w');
x = 1:size(S,2);
S = S([3 5 6],:);
bar(x,mean(S,1),'k','EdgeColor','w'); hold on
semilogx([x; x], [mean(S,1);mean(S,1)] + ...
         [std(S,[],1);  ...
         -std(S,[],1)],'r-','linewidth',2)
ylabel('Strength of evidence (S)','FontSize',16,'FontAngle','oblique')
xlabel('Major fascicles','FontSize',16,'FontAngle','oblique')
set(gca,'fontsize',16, ...
    'xlim', [0 21], 'xtick', x,  ... 
    'ylim', [0 160],'ytick', [0 80 160], ...
    'box',  'off',  'tickdir', 'out', 'ticklength', [0.025 0])
saveFig(fh,fullfile(savedir, figName),1)

end % Main function

%-------------------------------%
function saveFig(h,figName,eps)
if ~exist( fileparts(figName), 'dir'), mkdir(fileparts(figName));end
fprintf('[%s] saving figure... \n%s\n',mfilename,figName);

switch eps
    case {0,'jpeg'}
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