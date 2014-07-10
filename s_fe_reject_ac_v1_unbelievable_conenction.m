function s_fe_reject_ac_v1_unbelievable_conenction()
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
    'KK_96dirs_b2000_1p5iso', ...'JW_96dirs_b2000_1p5iso', ...
    };
if notDefined('saveDir'), saveDir = fullfile('/marcovaldo/frk/Dropbox','pestilli_etal_revision',mfilename);end
anatomyPath     = '/marcovaldo/frk/2t1/anatomy/';
datapath      = '/marcovaldo/frk/2t1/predator/';
fibs            = {'run01_fliprot_aligned_trilin_csd_lmax10_run01_fliprot_aligned_trilin_brainmask_run01_fliprot_aligned_trilin_wm_prob-500000_recomputed-rejected.mat', ...
                   'run01_fliprot_aligned_trilin_csd_lmax10_run01_fliprot_aligned_trilin_brainmask_run01_fliprot_aligned_trilin_wm_prob-500000_recomputed-optimized.mat'};
hemisphere = {'left','right'}; %left and right hemisphere

for iSbj = 1:length(subjects)
    % Find all the ROI computed for each subject:
    roiDir  = fullfile(anatomyPath,matchSubject2FSSUBJ(subjects{iSbj}),'label');
    % We will load three FreeSurfer ROIs (Pericalcarine, V1 and V2) We will
    % combine them to obtain a full representation of early visual cortex
    %pericalcarine = dir(fullfile(roiDir,'*pericalcarine*.mat'));
    V1 = dir(fullfile(roiDir,'*V1_label*_smooth*.mat'));
    %V2 = dir(fullfile(roiDir,'*V2*.mat'));    
    for ih = 1:length(hemisphere)
        for ifb = 1:length(fibs)
            tic, fprintf('\n[%s] Loading the fiber group... \n',mfilename)
            fg = fgRead(fullfile(sprintf('%s',datapath),subjects{iSbj},'fibers',fibs{ifb}));toc
            
            % Load the FreeSurfer V1
            tic, fprintf('\n[%s] Building a combined V1 (%s) | Anterior Commissure... \n',mfilename,V1(ih).name)
            %roiPeri = dtiReadRoi(fullfile(roiDir,pericalcarine(ih).name));
            roiV1 = dtiReadRoi(fullfile(roiDir,V1(ih).name));
            %roiV2 = dtiReadRoi(fullfile(roiDir,V2(ih).name));
            
            % Make a spherical ROi at the location fo the Anteriro commissure
            radius = 5;
            sphere_center = [0,0,0];
            roi2 = dtiNewRoi('AC_sphere');
            roi2.coords = dtiBuildSphereCoords(sphere_center, radius);
            roi  = dtiNewRoi('AC_V1');
            roi.coords = roiV1.coords;%[roiPeri.coords; roiV1.coords; roiV2.coords; roi2.coords];toc
            
            % Extract fibers touching early visual cortex and the AC
            tic, fprintf('\n[%s] Segmenting the TRACT connecting Visual Cortex and the Anterior Commissure... \n',mfilename)
            [tract{iSbj,ifb,ih}, keepFascicles{iSbj,ifb,ih}] = feSegmentFascicleFromConnectome(fg, {roiV1, roi2}, {'and','and'}, 'tmp');toc
            fprintf('\n[%s] Found tract with %i fibers. \n',mfilename,length(tract{iSbj,ifb,ih}.fibers))
            
            % Get some statistics:
            if ~isempty(tract{iSbj,ifb,ih}.fibers)                
                tic, fprintf('\n[%s] Getting information about the TRACT... \n',mfilename)
                fgInfo{iSbj,ifb,ih}.roiname = [V1(ih).name,'_',roi2.name, num2str(radius),'mm'];
                fgInfo{iSbj,ifb,ih}.volume = size(fefgGet(tract{iSbj,ifb,ih},'uniqueimagecoords'),1) * prod([1 1 1]); % mm-cube
                fgInfo{iSbj,ifb,ih}.length = fefgGet(tract{iSbj,ifb,ih},'length'); % in mm
                dtFileName    = fullfile(datapath,subjects{iSbj},'dtiInit','dt6.mat');
                fgInfo{iSbj,ifb,ih}.fa     = fefgGet(tract{iSbj,ifb,ih},'fa',dtFileName); % fractional anysotropy
            else
                tic, fprintf('\n[%s] TRACT is emapyt information set to nan... \n',mfilename) 
                fgInfo{iSbj,ifb,ih}.roiname = [V1(ih).name,'_',roi2.name, num2str(radius),'mm'];
                fgInfo{iSbj,ifb,ih}.volume = nan; % mm-cube
                fgInfo{iSbj,ifb,ih}.length = nan; % in mm
                fgInfo{iSbj,ifb,ih}.fa     = nan; % fractional anysotropy
                toc
            end
        end
    end
end

% Show and save the fascicles
for iSbj = 1:length(subjects)
    t1File        = fullfile(datapath,subjects{iSbj},'anatomy','t1.nii.gz');
    % Display the tract and save the figure
    % Load the T1 file for display
    t1     = niftiRead(t1File);
    for ih = 1:length(hemisphere)  
        fprintf('Subject: %i Hemisphere: %i\n',iSbj,ih)
        switch hemisphere{ih}
            case {'left'}
                axisLims   = [-60 5 -120 10 0 30];
                vw         = [-55,30];
                slices     = {[-3 0 0],[0 5 0],[0 0 -17 ]};
                histcolor{1}  = [0.4 0.4 0.4];
                histcolor{2}  = [.6 0.4 0.4];
                
            case {'right'}
                axisLims   = [-20 60 -90 5 0 30];
                vw         = [55,40];
                slices     = {[-2 0 0],[0 5 0],[0 0 -17]};
                histcolor{1}  = [0 0 0];
                histcolor{2}  = [.8 0.4 0.4];
                
            otherwise
                keyboard
        end

        for ifb = 1:length(fibs)
            if ~isempty(tract{iSbj,ifb,ih}.fibers)                
                % Show te new fiber group
                if ifb == 1
                    str = 'REJECTED';
                else
                    str = 'OPTIMIZED';
                end
                figName = sprintf('V1_AC_tract_%s_HEMI_%s_%s',str,hemisphere{ih},subjects{iSbj});
                fh = figure('name',figName);
                hold on
                h  = mbaDisplayBrainSlice(t1, slices{1});
                h  = mbaDisplayBrainSlice(t1, slices{2});
                h  = mbaDisplayBrainSlice(t1, slices{3});
                view(vw(1),vw(2)); axis(axisLims);
                tract2plot   = mbaComputeFibersOutliers(tract{iSbj,ifb,ih},2,1.4);

                [fh,ligh] = mbaDisplayConnectome(tract2plot.fibers,fh,[.98 .55 .45],'single',[], [], .6,[]);
                delete(ligh);camlight(hemisphere{ih}) 
                saveFig(fh,fullfile(saveDir,figName),'jpg')
            end
        end
    end
end

% Plot the FA and fiber length 
Opt.length = []; Opt.n = [];  Opt.fa = [];
Rej.length = []; Rej.n = [];  Rej.fa = [];
for iSbj = 1:length(subjects)
    for ih = 1:length(hemisphere)
        for ifb = 1:length(fibs)
            if ~isnan(fgInfo{iSbj,ifb,ih}.length)
                if ifb == 2 % Optimized
                    Opt.length = [Opt.length; mean(fgInfo{iSbj,ifb,ih}.length)];
                    Opt.n      = [Opt.n; length(fgInfo{iSbj,ifb,ih}.length)];
                    Opt.fa     = [Opt.fa; mean((cellfun(@mean,fgInfo{iSbj,ifb,ih}.fa)))]; 
                else % Rejected
                    Rej.length = [Rej.length; nanmean(fgInfo{iSbj,ifb,ih}.length)];
                    Rej.n      = [Rej.n; length(fgInfo{iSbj,ifb,ih}.length)];
                    Rej.fa     = [Rej.fa; nanmean((cellfun(@nanmean,fgInfo{iSbj,ifb,ih}.fa)))];                     
                end
            end
        end
    end
end

% Make a plot of the length of the fibers and of the FA
figName = sprintf('length_%s_%s',hemisphere{ih},subjects{iSbj});
fh = figure('name',figName,'color','w');
plot([1 1],[mean(Rej.length),mean(Rej.length)] + [-std(Rej.length)./sqrt(length(Rej.n)), std(Rej.length)./sqrt(length(Rej.n))],'r-')
hold on
plot([1 2],[mean(Rej.length),Opt.length],'ro')
set(gca,'tickdir','out','ticklen',[.01 0],'xlim',[0 3],'ylim',[100 140],'ytick',[100 120 140],'box','off','FontSize',14)
ylabel('Tract length (mm)','FontSize',14)
saveFig(fh,fullfile(saveDir,figName),'eps')

figName = sprintf('fa_%s_%s',hemisphere{ih},subjects{iSbj});
fh = figure('name',figName,'color','w');
plot([1 1],[mean(Rej.fa),mean(Rej.fa)] + [-std(Rej.fa)./sqrt(length(Rej.n)), std(Rej.fa)./sqrt(length(Rej.n))],'r-')
hold on
plot([1 2],[mean(Rej.fa),Opt.fa],'ro')
set(gca,'tickdir','out','ticklen',[.01 0],'xlim',[0 3],'ylim',[0.2 .6],'ytick',[.2 .4 .6],'box','off','FontSize',14)
ylabel('FA','FontSize',14)
saveFig(fh,fullfile(saveDir,figName),'eps')

figName = sprintf('n_%s_%s',hemisphere{ih},subjects{iSbj});
fh = figure('name',figName,'color','w');
plot([1 1],[mean(Rej.n),mean(Rej.n)] + [-std(Rej.n)./sqrt(length(Rej.n)), std(Rej.n)./sqrt(length(Rej.n))],'r-')
hold on
plot([1 2],[mean(Rej.n),Opt.n],'ro')
set(gca,'tickdir','out','ticklen',[.01 0],'xlim',[0 3],'ylim',[0 16],'ytick',[0 8 16],'box','off','FontSize',14)
ylabel('Number of fascicles','FontSize',14)
saveFig(fh,fullfile(saveDir,figName),'eps')

save(fullfile(saveDir),'fgInfo')
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