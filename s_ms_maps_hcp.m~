function s_ms_maps_hcp(trackingType)
%
% Load FE structeres obtained by preprocessing connectomesconstrained within a
% region of interest and within the cortex and makes some basic plot of
% statistics in the connectomes
%
% See also:

% Get the base directory for the data
datapath = '/marcovaldo/frk/2t2/HCP/';
subjects = {'115320','117122','118730'};

if notDefined('saveDir'), savedir = fullfile('/marcovaldo/frk/Dropbox/','pestilli_etal_revision',mfilename);end
if notDefined('trackingType'), trackingType = 'lmax10';end

% Bins for the fiber density estimates
xBins = [1 2 4 8 16 32 64 128 256 512 1024 2048];
x     = 1:length(xBins);

% Bins for the sum of weights estimates
wxBins = [.9./(2.^[10:-1:1]) ];
wx     = 1:length(wxBins);

doFD       = 1;
figVisible = 'off';
doMAPS = false;

for isbj = 1:length(subjects)
    % High-resolution Anatomy
    t1File = fullfile(datapath,subjects{isbj},'anatomy','/T1w_acpc_dc_restore_1p25.nii.gz');
    t1     = niftiRead(t1File);
    saveDir = fullfile(savedir,subjects{isbj});
    
    % File to load   
    connectomesPath   = fullfile(datapath,subjects{isbj},'connectomes');
    feFileToLoad = dir(fullfile(connectomesPath,sprintf('*%s*right.mat',trackingType)));
    fname = feFileToLoad.name(1:end-4);
    feFileToLoad = fullfile(connectomesPath,fname);

    fprintf('[%s] Loading: \n%s\n ======================================== \n\n',mfilename,feFileToLoad)
    load(feFileToLoad);
    fprintf('[%s] Extracting info: \n%s\n ======================================== \n\n',mfilename,feFileToLoad)
    coords  = feGet(fe,'roi coords');
    xform   = feGet(fe,'xform img 2 acpc');
    mapsize = feGet(fe, 'map size');
    if isempty(fe.fg)
        fe.path.savedir = fullfile('/marcovaldo/',fe.path.savedir(strfind(fe.path.savedir,'home')+length('home'):end));
        fiberPath = fullfile(fileparts(fe.path.savedir),'fibers');
        fibers    = dir(fullfile(fiberPath,sprintf('*%s*.pdb',trackingType)));
        fg = fgRead(fullfile(fiberPath,fibers.name));
    else
        fg = feGet(fe,'fibers acpc');
    end
    w       = feGet(fe,'fiber weights');
    fgOpt   = fgExtract(fg,w > 0,'keep');
    if isempty(fe.rep)
        if ~isempty(strfind(fe.path.dwifilerep,'home'))
            fe.path.dwifilerep = fullfile('/marcovaldo/',fe.path.dwifilerep(strfind(fe.path.dwifilerep,'home')+length('home'):end));
        end
        fe = feConnectomeSetDwi(fe,fe.path.dwifilerep,true);
    end
    rmseM   = feGetRep(fe, 'vox rmse');
    rmseD   = feGetRep(fe, 'vox rmse data');
    rmseR = feGetRep(fe, 'vox rmse ratio');
    slice   = 2;%[-80:4:-2 2:4:80];
    clear fe
    
    if doFD  
        fprintf('[%s] Computing CANDIDATE fiber density: \n%s\n ======================================== \n\n',mfilename,feFileToLoad)
        % Get the fiber density
        % fd = feGet(fe,'fiber density');
        fdImg  = dtiComputeFiberDensityNoGUI(fg, xform, mapsize);  
        fprintf('[%s] Computing OPTIMAZED fiber density: \n%s\n ======================================== \n\n',mfilename,feFileToLoad)

        fdOImg = dtiComputeFiberDensityNoGUI(fgOpt, xform, mapsize);     
        fprintf('[%s] Computing WEIGHTED fiber density: \n%s\n ======================================== \n\n',mfilename,feFileToLoad)

        fdWImg = dtiComputeFiberDensityNoGUI(fgOpt, xform, mapsize,[],[],[],[],w);
    end
    
    if doMAPS
    fprintf('[%s] Making maps: \n%s\n ======================================== \n\n',mfilename,feFileToLoad)
    % Make a plot of the maps
    for is = 1:length(slice)
        
        % RMSE off the model
        map='hot';maxRmse = 90;
        figName = sprintf('RMSE_ModelCoronal_%s_slice%i',fname,slice(is));
        fh = figure('name',figName,'visible',figVisible,'color','w');
        rmseImg = feReplaceImageValues(nan(mapsize),rmseM,coords+1);
        %rmseImg(rmseImg < 1)=nan;
        rmseImg(rmseImg > maxRmse) = maxRmse;
        ni  = niftiCreate('data',rmseImg, 'fname', figName, ...
            'qto_xyz',xform, ...
            'fname','FDM', ...
            'data_type',class(rmseImg));
        
        sh = mbaDisplayOverlay(t1, ni, [0 slice(is) 0], [], map);
        saveMapCoronal(fh,figName,saveDir,nanmean(rmseImg(:)),nanmedian(rmseImg(:)),nanstd(rmseImg(:)),maxRmse,map)

        figName = sprintf('RMSE_ModelSagital_%s_slice%i',fname,slice(is));
        sh = mbaDisplayOverlay(t1, ni, [slice(is) 0 0], [], map);
        saveMapSagital(fh,figName,saveDir,nanmean(rmseImg(:)),nanmedian(rmseImg(:)),nanstd(rmseImg(:)),maxRmse,map)

        % Ratio rmse
        map = 'jet';maxRRmse = 2;minRRmse = 0.125;
        figName = sprintf('Ratio_rmseCoronal_%s_slice%i',fname,slice(is));
        fh = figure('name',figName,'visible',figVisible,'color','w');
        rrmseImg = feReplaceImageValues(nan(mapsize),rmseR,coords);
        %rrmseImg(rrmseImg > 0.73)=nan;
        rrmseImg(rrmseImg > maxRRmse)=maxRRmse;
        ni  = niftiCreate('data',rrmseImg, 'fname', figName, ...
            'qto_xyz',xform, ...
            'fname','FDM', ...
            'data_type',class(rrmseImg));
        sh = mbaDisplayOverlay(t1, ni, [0 slice(is) 0], [], map);
        saveMapCoronal(fh,figName,saveDir,nanmean(rrmseImg(:)),nanmedian(rrmseImg(:)),nanstd(rrmseImg(:)),maxRRmse,map)
                
        figName = sprintf('Ratio_rmseSagital_%s_slice%i',fname,slice(is));
        sh = mbaDisplayOverlay(t1, ni, [slice(is) 0 0], [], map);
        saveMapSagital(fh,figName,saveDir,nanmean(rrmseImg(:)),nanmedian(rrmseImg(:)),nanstd(rrmseImg(:)),maxRRmse,map)
        
        % Fiber density maps:
        % This will be used to normalize the rage of the fiber density across plots
        minfd = 2;   % Min fiber density
        maxfd = 256; % Max fiber density
        
        % Optimized connectome
        map     = 'jet';
        figName = sprintf('FibDensMapCoronalCan_%s_%s_slice%i',fname,'candidate',slice(is));
        fh      = figure('name',figName,'visible',figVisible,'color','w');
        
        % Candidate connectome
        fdImg(fdImg==0) = nan;
        fdImg(isnan(rrmseImg))=nan; 
        fdImg(fdImg > maxfd)=maxfd;
        fdImg(1,1) = maxfd;
        fdImg(1,2) = minfd;
        ni  = niftiCreate('data',(fdImg) , 'fname', figName, ...
            'qto_xyz',xform, ...
            'fname','RRMSE', ...
            'data_type',class(fdImg));
        sh = mbaDisplayOverlay(t1, ni, [0 slice(is) 0], [],map);
        saveMapCoronal(fh,figName,saveDir,nanmean(fdImg(:)),nanmedian(fdImg(:)),nanstd(fdImg(:)),maxfd ,map)
          
        figName = sprintf('FibDensMapSagitalCan_%s_%s_slice%i',fname,'candidate',slice(is));
        sh = mbaDisplayOverlay(t1, ni, [slice(is) 0 0], [],map);
        saveMapSagital(fh,figName,saveDir,nanmean(fdImg(:)),nanmedian(fdImg(:)),nanstd(fdImg(:)),maxfd ,map)

        % Optimized connectome
        figName = sprintf('FibDensMapCoronalOpt_%s_%s_slice%i',fname,'optimized',slice(is));
        fh  = figure('name',figName,'visible',figVisible,'color','w');
        fdOImg(1,1) = maxfd;
        fdOImg(1,2) = minfd;
        fdOImg(fdOImg==0) = nan; 
        fdOImg(fdOImg > maxfd)=maxfd;
        fdOImg(isnan(rrmseImg))=nan;
        ni  = niftiCreate('data',(fdOImg), ...
            'qto_xyz',xform, ...
            'fname','FDM', ...
            'data_type',class(fdOImg));
        sh = mbaDisplayOverlay(t1, ni, [0 slice(is) 0], [],map);
        saveMapCoronal(fh,figName,saveDir,nanmean(fdOImg(:)),nanmedian(fdOImg(:)),nanstd(fdOImg(:)),maxfd ,map)
         
        figName = sprintf('FibDensMapSagitalOpt_%s_%s_slice%i',fname,'optimized',slice(is));
        sh = mbaDisplayOverlay(t1, ni, [slice(is) 0 0], [],map);
        saveMapSagital(fh,figName,saveDir,nanmean(fdOImg(:)),nanmedian(fdOImg(:)),nanstd(fdOImg(:)),maxfd ,map)
   
        % Weight density (sum of weights)
        map = 'hsv';
        figName = sprintf('WeightMapCoronal_%s_%s_slice%i',fname,'optimized',slice(is));
        fh  = figure('name',figName,'visible',figVisible,'color','w');
        fdWImg(isnan(rrmseImg))=nan;
        ni  = niftiCreate('data',fdWImg, 'fname', figName, ...
            'qto_xyz',xform, ...
            'fname','FDM', ...
            'data_type',class(fdWImg));
        sh = mbaDisplayOverlay(t1, ni, [0 slice(is) 0], [],map);
        saveMapCoronal(fh,figName,saveDir,nanmean(fdWImg(:)),nanmedian(fdWImg(:)),nanstd(fdWImg(:)),max(fdWImg(:)) ,map)
        

        figName = sprintf('WeightMapSagital_%s_%s_slice%i',fname,'optimized',slice(is));
        sh = mbaDisplayOverlay(t1, ni, [slice(is) 0 0], [],map);
        saveMapSagital(fh,figName,saveDir,nanmean(fdWImg(:)),nanmedian(fdWImg(:)),nanstd(fdWImg(:)),max(fdWImg(:)) ,map)
       
        % RMSE of the data
        figName = sprintf('RMSE_DataCoronal_%s_slice%i',fname,slice(is));
        fh = figure('name',figName,'visible',figVisible,'color','w');
        rmseImg = feReplaceImageValues(nan(mapsize),rmseD,coords);
        rmseImg(isnan(rrmseImg))=nan;
        rmseImg(rmseImg>maxRmse) = maxRmse;
        ni  = niftiCreate('data',rmseImg, 'fname', figName, ...
            'qto_xyz',xform, ...
            'fname','FDM', ...
            'data_type',class(rmseImg));
        sh = mbaDisplayOverlay(t1, ni, [0 slice(is) 0], [], map);
        saveMapCoronal(fh,figName,saveDir,nanmean(rmseImg(:)),nanmedian(rmseImg(:)),nanstd(rmseImg(:)),maxRmse,map)
                  
        figName = sprintf('RMSE_DataSagital_%s_slice%i',fname,slice(is));
        sh = mbaDisplayOverlay(t1, ni, [slice(is) 0 0], [], map);
        saveMapSagital(fh,figName,saveDir,nanmean(rmseImg(:)),nanmedian(rmseImg(:)),nanstd(rmseImg(:)),maxRmse,map)
  
        close all
        drawnow
    end
    end
    
    % Histogram plots
    figName = sprintf('FibDensHistCandVSOpt_%s',fname);
    fh  = figure('name',figName,'visible',figVisible,'color','w');
    x = [1 2.^[1 2 3 4 5 6 7 8 9 10]];
    fdOImg(fdOImg == 0) = nan;
    fdImg(fdImg == 0) = nan;
    [yFD(isbj,:),xFD(isbj,:)]  = hist(fdImg(:),x);
    yFD(isbj,:) = 100*yFD(isbj,:)./sum(yFD(isbj,:));
    [yoFD(isbj,:),xoFD(isbj,:)]= hist(fdOImg(:),x);
     yoFD(isbj,:) = 100*yoFD(isbj,:)./sum(yoFD(isbj,:));
    semilogx(xFD(isbj,:),yFD(isbj,:),'k-',xoFD(isbj,:),yoFD(isbj,:),'r-','linewidth',2)
    ylabel('Percent voxels','FontSize',16,'FontAngle','oblique')
    xlabel('Fascicles per voxel','FontSize',16,'FontAngle','oblique')
    legend(gca,{'Candidated','Optimized'},'box','off')
    set(gca,'fontsize',16, ...
        'ylim', [0 30], ...
        'ytick',[0 15 30], ...
        'xlim', [0.5 2^10],'xtick',x,...
        'box','off','tickdir','out','ticklength',[0.025 0])
    saveFig(fh,fullfile(saveDir, figName),1)
    
    figName = sprintf('RMSE_HistDataVSOpt_%s',fname);
    fh  = figure('name',figName,'visible',figVisible,'color','w');
    x   = 0:10:700;
   [yoRMSE(isbj,:),xoRMSE(isbj,:)]= hist(rmseM(:),x);    
    yoRMSE(isbj,:) = 100*yoRMSE(isbj,:)./sum(yoRMSE(isbj,:));
    plot(xoRMSE(isbj,:),yoRMSE(isbj,:),'r-','linewidth',2)
    ylabel('Percent voxels','FontSize',16,'FontAngle','oblique')
    xlabel('RMSE (raw scanner units)','FontSize',16,'FontAngle','oblique')
    legend(gca,{'Model'},'box','off')
    set(gca,'fontsize',16, ...
        'ylim', [0 14], ...
        'ytick',[0 7 14], ...
        'xlim', [0 700],'xtick',[0 350 700],...
        'box','off','tickdir','out','ticklength',[0.025 0])
    saveFig(fh,fullfile(saveDir, figName),1)
end

% Average histograms
saveDir = fullfile(savedir,'average_hcp_1p25mm');

% Save the results to file, it takes along time to load all these FE strctures...
m.density.candidatey  = mean(yFD,1);
m.density.candidateSte= std(yFD,[],1)./sqrt(size(yFD,1));
m.density.optimaly    = mean(yoFD,1);
m.density.optimalSte  = std(yoFD,[],1)./sqrt(size(yoFD,1));
m.density.x     = xFD(isbj,:);
m.density.units = {'x=Fascicles per voxel','y=percent voxels'};
m.density.yFD   = yFD;
m.density.yoFD  = yoFD;

% rmse data vs. model
m.rmse.model   = mean(yoRMSE,1);
m.rmse.modelSte= std(yoRMSE,[],1)./sqrt(size(yoRMSE,1));
m.rmse.x       = xoRMSE(isbj,:);
m.rmse.units = {'x=rmse (raw scanner units)','y=percent voxels'};
m.rmse.yoRMSE=yoRMSE;

mkdir(saveDir)
save(fullfile(saveDir,'mean_histograms.mat'),'m')

% Histogram plots
figName = sprintf('FibDensHistCandVSOpt_%s',fname);
fh  = figure('name',figName,'visible',figVisible,'color','w');
semilogx(m.density.x,m.density.candidatey,'k-',m.density.x,m.density.optimaly,'r-','linewidth',2)
hold on
semilogx([m.density.x;m.density.x], [m.density.candidatey-m.density.candidateSte; ...
                                     m.density.candidatey+m.density.candidateSte],'k-', ...
         [m.density.x;m.density.x], [m.density.optimaly-m.density.optimalSte; ...
                                    m.density.optimaly+m.density.optimalSte],'r-','linewidth',2)

ylabel('Percent voxels','FontSize',16,'FontAngle','oblique')
xlabel('Fascicles per voxel','FontSize',16,'FontAngle','oblique')
legend(gca,{'Candidated','Optimized'},'box','off')
set(gca,'fontsize',16, ...
    'ylim', [0 30], ...
    'ytick',[0 15 30], ...
    'xlim', [0.5 2^10],'xtick',[m.density.x],...
    'box','off','tickdir','out','ticklength',[0.025 0])
saveFig(fh,fullfile(saveDir, figName),1)

figName = sprintf('RMSE_mean_HistDataVSOpt_%s',fname);
fh  = figure('name',figName,'visible',figVisible,'color','w');
plot(m.rmse.x,m.rmse.model,'r-','linewidth',2)
hold on
plot([m.rmse.x;m.rmse.x],[m.rmse.model-m.rmse.modelSte;m.rmse.model+m.rmse.modelSte],'r-','linewidth',2)

ylabel('Percent voxels','FontSize',16,'FontAngle','oblique')
xlabel('RMSE (raw scanner units)','FontSize',16,'FontAngle','oblique')
    legend(gca,{'Model'},'box','off')
    set(gca,'fontsize',16, ...
        'ylim', [0 14], ...
        'ytick',[0 7 14], ...
        'xlim', [0 700],'xtick',[0 350 700],...
        'box','off','tickdir','out','ticklength',[0.025 0])
saveFig(fh,fullfile(saveDir, figName),1)

end  % Main function


%---------------------------------%
function saveMapSagital(fh,figName,saveDir,M,m,SD,maxfd,map)
% This helper function saves two figures for each map and eps with onlythe
% axis and a jpg with only the brain slice.
% The two can then be combined in illustrator.
%
% First we save only the slice as jpeg.
set(gca,'fontsize',16,'ytick',[-80 -40 0 40 80], ...
    'ztick',[-40  0  40  80], ...
    'xlim',[-80 80],'ylim',[-110 100],'zlim',[-60 80],'tickdir','out','ticklength',[0.025 0])
axis off
saveFig(fh,fullfile(saveDir,figName),'tiff')
saveFig(fh,fullfile(saveDir,figName),'png')

% Then we save the slice with the axis as
% eps. This will only generate the axis
% that can be then combined in illustrator.
axis on
grid off

title(sprintf('mean %2.2f | median %2.2f | SD %2.2f', ...
    M,m,SD),'fontsize',16,'FontAngle','oblique')
zlabel('Z (mm)','fontsize',16,'FontAngle','oblique')
xlabel('X (mm)','fontsize',16,'FontAngle','oblique')
cmap = colormap(eval(sprintf('%s(255)',map)));
colorbar('ytick',linspace(0,1,5),'yticklabel', ...
    {1, num2str(ceil(maxfd/8)), num2str(ceil(maxfd/4)), ...
    num2str(ceil(maxfd/2)), num2str(ceil(maxfd))}, ...
    'tickdir','out','ticklength',[0.025 0],'fontsize',16)
saveFig(fh,fullfile(saveDir,figName),1)
end

%---------------------------------%
function saveMapCoronal(fh,figName,saveDir,M,m,SD,maxfd,map)
% This helper function saves two figures for each map and eps with onlythe
% axis and a jpg with only the brain slice.
% The two can then be combined in illustrator.
%
% First we save only the slice as jpeg.
set(gca,'fontsize',16,'ztick',[-20 -10 0 10 20], ...
    'xtick',[0 10 20 30 40 50], ...
    'xlim',[-5 70],'zlim',[-30 40],'tickdir','out','ticklength',[0.025 0])
axis off
saveFig(fh,fullfile(saveDir, figName),'tiff')
saveFig(fh,fullfile(saveDir, figName),'png')

% Then we save the slice with the axis as
% eps. This will only generate the axis
% that can be then combined in illustrator.
axis on
grid off

title(sprintf('mean %2.2f | median %2.2f | SD %2.2f', ...
    M,m,SD),'fontsize',16,'FontAngle','oblique')
zlabel('Z (mm)','fontsize',16,'FontAngle','oblique')
xlabel('X (mm)','fontsize',16,'FontAngle','oblique')
cmap = colormap(eval(sprintf('%s(255)',map)));
colorbar('ytick',linspace(0,1,5),'yticklabel', ...
    {1, num2str(ceil(maxfd/8)), num2str(ceil(maxfd/4)), ...
    num2str(ceil(maxfd/2)), num2str(ceil(maxfd))}, ...
    'tickdir','out','ticklength',[0.025 0],'fontsize',16)
saveFig(fh,fullfile(saveDir, figName),1)
end

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