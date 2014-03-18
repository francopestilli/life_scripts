function s_pestilli_etal_fig8_major_fascicles_tests()
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
    fasciclesPath  = fullfile(datapath,subjects{isbj},'afq');
    %fasciclesFiles = dir(fullfile(fasciclesPath,'*.mat'));
    
    % Load the FE structure
    disp('Loading the FE strcuture...')
    connectomesPath = fullfile(datapath,subjects{isbj},'connectomes');
    feFileToLoad    = dir(fullfile(connectomesPath,sprintf('*%s*prob*.mat',trackingType)));
    fname           = feFileToLoad(1).name(1:end-4);
    feFileToLoad    = fullfile(connectomesPath,fname);
    load(feFileToLoad)
    disp('Loading the whole-brain fiber group...')
    if isempty(fe.fg)
        if ~isempty(strfind(fe.path.savedir,'home'))
            fe.path.savedir = fullfile('/marcovaldo/',fe.path.savedir(strfind(fe.path.savedir,'home')+length('home'):end));
        end
        fiberPath = fullfile(fileparts(fe.path.savedir),'fibers');
        fibers    = dir(fullfile(fiberPath,sprintf('*%s*.pdb',trackingType)));
        fe = feSet(fe,'fg from acpc',fgRead(fullfile(fiberPath,fibers.name)));
    end
    if isempty(fe.rep)
        if ~isempty(strfind(fe.path.dwifilerep,'home'))
            fe.path.dwifilerep = fullfile('/marcovaldo/',fe.path.dwifilerep(strfind(fe.path.dwifilerep,'home')+length('home'):end));
        end
        fe = feConnectomeSetDwi(fe,fe.path.dwifilerep,true);
    end

    tic, fprintf('[%s] Culling the whole-brain fiber group...',mfilename)
    fe = feConnectomeReduceFibers(fe,feGet(fe,'fiber weights') > 0);toc
    
    % Segment the fibers using AFQ
    dtFile    = fullfile(datapath,subjects{isbj},'dtiInit','dt6.mat');
    [fg_classified,~,classification]= AFQ_SegmentFiberGroups(dtFile, feGet(fe,'fibers acpc'),[],[],false);
    
    % Split the fiber groups into individual groups
    fascicles = fg2Array(fg_classified);

    % Load the fascicles
    for iFas = 8%1:length(fascicles)    
        fprintf('Performing visrtual lesion: %s...\n',classification.names{iFas})

        % Get the fibers for the current fascicle
        fascicles2keep = find(classification.index==iFas);
        % Make sure that the number of fascicles we are removin is similar
        % to the number of fascicles saved in the fiber group, this is just
        % a consistency check
        if ~(length(fascicles2keep)==length(fascicles(iFas).fibers)); keyboard, end

        % Reduce the connectome to the voxels and fibers of the connection between
        [S(isbj,iFas), figHandle] = feVirtualLesion(fe,fascicles2keep,0);

        classification.names{iFas}(isspace(classification.names{iFas}))='_';
        figName = fullfile(savedir,subjects{isbj},sprintf('Distributions_%s_%s',classification.names{iFas},fname));
        saveFig(figHandle,figName,true)
    end
    clear fascicles clear fg_classified dtFile fe 
end
            
save(fullfile(savedir,['average_s_',fname,'.mat']),'S')

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