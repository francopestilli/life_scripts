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
figVisible = 'off';

% Get the base directory for the data
datapath = '/marcovaldo/frk/2t1/predator/';
subjects = {...
    'HT_96dirs_b2000_1p5iso', ...
    'KW_96dirs_b2000_1p5iso', ...'JW_96dirs_b2000_1p5iso', ...
    'KK_96dirs_b2000_1p5iso', ...
    'MP_96dirs_b2000_1p5iso', ...
    'FP_96dirs_b2000_1p5iso', ...
    };
if notDefined('saveDir'), savedir = fullfile('/marcovaldo/frk/Dropbox','pestilli_etal_revision',mfilename);end
if notDefined('trackingType'), trackingType = 'lmax10';end
recompute= 0;

if recompute
    for isbj = 1:length(subjects)
        fasciclesPath  = fullfile(datapath,subjects{isbj},'afq');
        %fasciclesFiles = dir(fullfile(fasciclesPath,'*.mat'));
        
        % Load the FE structure
        disp('Loading the FE strcuture...')
        connectomesPath = fullfile(datapath,subjects{isbj},'connectomes');
        feFileToLoad    = dir(fullfile(connectomesPath,sprintf('*%s*prob*recomputed.mat',trackingType)));
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
            fe = feSet(fe,'fg from acpc',fgRead(fullfile(fiberPath,fibers(1).name)));
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
        save(fullfile(fasciclesPath,['fascicles_',fname,'.mat']),'fascicles','classification')
        
        % Load the fascicles
        for iFas = 1:length(fascicles)
            fprintf('Performing virtual lesion: %s...\n',classification.names{iFas})
            
            % Get the fibers for the current fascicle
            fascicles2keep = find(classification.index==iFas);
            % Make sure that the number of fascicles we are removin is similar
            % to the number of fascicles saved in the fiber group, this is just
            % a consistency check
            if ~(length(fascicles2keep)==length(fascicles(iFas).fibers)); keyboard, end
            if ~isempty(fascicles2keep)
                % Reduce the connectome to the voxels and fibers of the connection between
                display.distributions = true;
                display.tract         = false;
                display.evidence      = true;
                [se(isbj,iFas), fig] = feVirtualLesion(fe, fascicles2keep, display,false);
                for ifs = 1:length(fig)
                    if ~isnan(fig(ifs).h)
                        classification.names{iFas}(classification.names{iFas}==' ')='_';
                        tag = classification.names{iFas};
                        saveFig(fig(ifs).h,fullfile(fasciclesPath,fig(ifs).name,tag), fig(ifs).type);
                    end
                end
                close all
            else
                se(isbj,iFas).s.mean = nan;
                se(isbj,iFas).s.std = nan;
                se(isbj,iFas).em.mean = nan;
                se(isbj,iFas).j.mean = nan;
                se(isbj,iFas).kl.mean = nan;
            end
        end
        clear fascicles clear fg_classified dtFile fe
    end
    save(fullfile(savedir,['average_strength_of_evidence_',fname,'.mat']),'se')
else
    fname = 'average_strength_of_evidence_run01_fliprot_aligned_trilin_csd_lmax10_run01_fliprot_aligned_trilin_brainmask_run01_fliprot_aligned_trilin_wm_prob-500000_recomputed';
    load(fullfile(savedir,[fname,'.mat']),'se')
end

% Reorganize the strength of evidence and emd cross individuals and fscicles
for is = 1:size(se,1)
    for ifs = 1:size(se,2)
        s.m(is,ifs)  = se(is,ifs).s.mean;
        s.sd(is,ifs) = se(is,ifs).s.std;
        s.em(is,ifs) = se(is,ifs).em.mean;
        s.j(is,ifs)  = se(is,ifs).j.mean;
        s.kl(is,ifs) = se(is,ifs).kl.mean;
        
    end
end

% Make a plot of the Strength of evidence for each fascicle across subjects.
figName = sprintf('Strength_evidence_across_subjects_%s',  fname);
fh  = figure('name',figName,'visible',figVisible,'color','w');
x = 1:size(s.m,2);
bar(x,nanmean(s.m,1),'k','EdgeColor','w'); hold on
semilogx([x; x], [nanmean(s.m,1);nanmean(s.m,1)] + ...
         [nanstd(s.m,[],1)./sqrt(sum(~isnan(s.m)));  ...
         -nanstd(s.m,[],1)./sqrt(sum(~isnan(s.m)))],'r-','linewidth',2)
ylabel('Strength of evidence (s.d.)','FontSize',16,'FontAngle','oblique')
xlabel('Major fascicles','FontSize',16,'FontAngle','oblique')
set(gca,'fontsize',16, ...
    'xlim', [0 21], 'xtick', x,  ... 
    'ylim', [0 100],'ytick', [0 50 100], ...
    'box',  'off',  'tickdir', 'out', 'ticklength', [0.025 0])
saveFig(fh,fullfile(savedir, figName),1)

% Make a plot of the Strength of evidence for each fascicle across subjects.
figName = sprintf('Earth_movers_distance_across_individuals_%s',  fname);
fh  = figure('name',figName,'visible',figVisible,'color','w');
x = 1:size(s.em,2);
bar(x,nanmean(s.em,1),'k','EdgeColor','w'); hold on
semilogx([x; x], [nanmean(s.em,1);nanmean(s.em,1)] + ...
         [nanstd(s.em,[],1)./sqrt(sum(~isnan(s.m)));  ...
         -nanstd(s.em,[],1)./sqrt(sum(~isnan(s.m)))],'r-','linewidth',2)
ylabel('Earth movers distance (scanner units)','FontSize',16,'FontAngle','oblique')
xlabel('Major fascicles','FontSize',16,'FontAngle','oblique')
set(gca,'fontsize',16, ...
    'xlim', [0 21], 'xtick', x,  ... 
    'ylim', [0 26],'ytick', [0 14 26], ...
    'box',  'off',  'tickdir', 'out', 'ticklength', [0.025 0])
saveFig(fh,fullfile(savedir, figName),1)

% Make a plot of the Strength of evidence for each fascicle across subjects.
figName = sprintf('Jeffrey_divergence_across_individuals_%s',  fname);
fh  = figure('name',figName,'visible',figVisible,'color','w');
x = 1:size(s.j,2);
bar(x,nanmean(s.j,1),'k','EdgeColor','w'); hold on
semilogx([x; x], [nanmean(s.j,1);nanmean(s.j,1)] + ...
         [nanstd(s.j,[],1)./sqrt(sum(~isnan(s.m)));  ...
         -nanstd(s.j,[],1)./sqrt(sum(~isnan(s.m)))],'r-','linewidth',2)
ylabel('Jeffrey''s divergence (bits)','FontSize',16,'FontAngle','oblique')
xlabel('Major fascicles','FontSize',16,'FontAngle','oblique')
set(gca,'fontsize',16, ...
    'xlim', [0 21], 'xtick', x,  ... 
    'ylim', [0 2],'ytick', [0 1 2], ...
    'box',  'off',  'tickdir', 'out', 'ticklength', [0.025 0])
saveFig(fh,fullfile(savedir, figName),1)

% Make a plot of the Strength of evidence for each fascicle across subjects.
figName = sprintf('KL_divergence_across_individuals_%s',  fname);
fh  = figure('name',figName,'visible',figVisible,'color','w');
x = 1:size(s.kl,2);
bar(x,nanmean(s.kl,1),'k','EdgeColor','w'); hold on
semilogx([x; x], [nanmean(s.kl,1);nanmean(s.kl,1)] + ...
         [nanstd(s.kl,[],1)./sqrt(sum(~isnan(s.m)));  ...
         -nanstd(s.kl,[],1)./sqrt(sum(~isnan(s.m)))],'r-','linewidth',2)
ylabel('K-L divergence (bits)','FontSize',16,'FontAngle','oblique')
xlabel('Major fascicles','FontSize',16,'FontAngle','oblique')
set(gca,'fontsize',16, ...
    'xlim', [0 21], 'xtick', x,  ... 
    'ylim', [0 8],'ytick', [0 4 8], ...
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