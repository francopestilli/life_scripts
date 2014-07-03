function fe = s_fe_probabilistic_vs_deterministic_96()
%
% This function llustrates how to:
%  - initialize a LIFE structure from a candidate connectome
%  - Generate an optimized connectome from a cadidate connectome using the
%  LIFE strustrue
%
%  fe = s_fe_fit()
% 
% INPUTS:  none
% OUTPUTS: fe structure the optimized life structure
%
% Copyright Franco Pestilli (2013) Vistasoft Stanford University.
% Get the base directory for the data
datapath = '/marcovaldo/frk/2t1/predator/';
subjects = {...
    'HT_96dirs_b2000_1p5iso', ...    'MP_96dirs_b2000_1p5iso', ...
    'KW_96dirs_b2000_1p5iso', ...
    'KK_96dirs_b2000_1p5iso', ...'JW_96dirs_b2000_1p5iso', ...
    'FP_96dirs_b2000_1p5iso', ...
    };

if notDefined('saveDir'), savedir = fullfile('/marcovaldo/frk/Dropbox','pestilli_etal_revision',mfilename);end
if notDefined('trackingType'), trackingType = 'lmax10';end
display.evidence = 1;
display.distributions = 1;

for isbj = 1:length(subjects)
    saveDir = fullfile(savedir,subjects{isbj});
    
    %% Load two pre-culled connectomes
    feProbFileName = sprintf('*%s*prob*recomputed.mat',trackingType);
    feDetFileName  = '*tensor*recomputed.mat';
    connectomesPath   = fullfile(datapath,subjects{isbj},'connectomes');
    feFileToLoad = dir(fullfile(connectomesPath,feProbFileName));
    fname = feFileToLoad(1).name(1:end-4);
    feProbFileToLoad = fullfile(connectomesPath,fname);
    
    feFileToLoad = dir(fullfile(connectomesPath,feDetFileName));
    fname = feFileToLoad.name(1:end-4);
    feDetFileToLoad = fullfile(connectomesPath,fname);
     
    % Extract the RMSE, R_rmse and the coordinates of the white matter
    fprintf('[%s] Loading: \n%s\n ======================================== \n\n',mfilename,feProbFileToLoad)
    load(feProbFileToLoad);
    if isempty(fe.rep)
        if ~isempty(strfind(fe.path.dwifilerep,'home'))
            fe.path.dwifilerep = fullfile('/marcovaldo/',fe.path.dwifilerep(strfind(fe.path.dwifilerep,'home')+length('home'):end));
        end
        fe = feConnectomeSetDwi(fe,fe.path.dwifilerep,true);
    end
    p.rmse   = feGetRep(fe,'vox rmse');
    p.rrmse  = feGetRep(fe,'vox rmse ratio');
    p.coords = feGet(   fe,'roi coords');
    clear fe
        
    % Extract the RMSE, R_rmse and the coordinates of the white matter
    fprintf('[%s] Loading: \n%s\n ======================================== \n\n',mfilename,feDetFileToLoad)
    load(feDetFileToLoad);
    if isempty(fe.rep)
        if ~isempty(strfind(fe.path.dwifilerep,'home'))
            fe.path.dwifilerep = fullfile('/marcovaldo/',fe.path.dwifilerep(strfind(fe.path.dwifilerep,'home')+length('home'):end));
        end
        fe = feConnectomeSetDwi(fe,fe.path.dwifilerep,true);
    end
    d.rmse   = feGetRep(fe, 'vox rmse');
    d.rrmse  = feGetRep(fe, 'vox rmse ratio');
    d.coords = feGet(   fe, 'roi coords');
    clear fe
    
    %% Find the common coordinates between the two connectomes
    %
    % There are more coordinates in the Prob conectome, because the tracking
    % fills up more White-matter.
    %
    % So, first we find the indices in the probabilistic connectome of the
    % coordinate in the deterministic conenctome.
    %
    % But there are some of the coordinates in the Deterministic conectome that
    % are NOT in the Probabilistic connectome.
    %
    % So, second we find the indices in the Deterministic connectome of the
    % subset of coordinates in the Probabilistic connectome found in the
    % previous step.
    
    % First we find the coordinates in the Probabilistic conectome that are
    % also in the Deterministic connectome.
    prob.coordsIdx = ismember(p.coords,d.coords,'rows');
    
    % Second we find the coordinates in the Deterministic connectome that are
    % also in the Probabilistic connectome.
    prob.coords   = p.coords(prob.coordsIdx,:);
    det.coordsIdx = ismember(d.coords,prob.coords,'rows');
    det.coords    = d.coords(det.coordsIdx,:);
    
    % What we really need is detCoordsIdx and probCoordsIdx. These allow us to
    % find the common voxel indices in rmse and rrmse, etc.
    prob.rmse  = p.rmse( prob.coordsIdx);
    prob.rrmse = p.rrmse(prob.coordsIdx);
    
    det.rmse  = d.rmse( det.coordsIdx);
    det.rrmse = d.rrmse(det.coordsIdx);
    
    %% Strength of evidence
    se = feComputeEvidence(prob.rmse,det.rmse);
    if display.distributions
        % Raw RMSE distirbutions
        fig(1).name = sprintf('rmse_distributions_%s',mfilename);
        fig(1).fh   = figure('name',fig(1).name,'color','w');
        set(fig(1).fh,'Units','normalized','Position',[0.007 0.55  0.28 0.36]);
        plot(se.lesion.xhist,se.lesion.hist,'-','color', [.95 .45 .1],'linewidth',2); hold on
        plot(se.nolesion.xhist,se.nolesion.hist,'-','linewidth',2, 'color', [.1 .45 .95])
        plot([se.nolesion.rmse.mean,se.nolesion.rmse.mean], [0,0.2],'-','color',[.1 .45 .95] )
        plot([se.lesion.rmse.mean,se.lesion.rmse.mean], [0,0.2], '-', 'color',[.95 .45 .1])
        title(sprintf('mean RMSE\nno-lesion %2.3f | lesion %2.2f', ...
            se.nolesion.rmse.mean,se.lesion.rmse.mean),'fontsize',16)
        ylabel('Probability', 'fontsize',14);xlabel('RMSE', 'fontsize',14)
        legend({'Lesion','No lesion'},'fontsize',14);
        set(gca,'box','off','xtick',[0 round(se.xrange(2)/2) se.xrange(2)], ...
            'ytick',[0 .1 .2],'xlim',[0 se.xrange(2)],'ylim',[0 .25], ...
            'tickdir', 'out', 'ticklength', [0.025 0]);
        
        % Plot the null distribution and the empirical difference
        ywo_e  = se.s.lesioned_e;
        y_e    = se.s.unlesioned_e;
        woxhis = se.s.lesioned.xbins;
        xhis   = se.s.unlesioned.xbins;
        min_x  = se.s.min_x;
        max_x  = se.s.max_x;
        fig(2).name = sprintf('virtual_lesion_test_mean_rmse_hist_%s_%s',mfilename,'prob_det');
        fig(2).fh   = figure('name',fig(2).name,'color','w');
        set(fig(2).fh,'Units','normalized','Position',[0.007 0.55  0.28 0.36]);
        patch([xhis,xhis],y_e(:),[.1 .45 .95],'FaceColor',[.1 .45 .95],'EdgeColor',[.1 .45 .95]); % Distribution as the +/- 2SD
        hold on
        patch([woxhis,woxhis],ywo_e(:),[.95 .45 .1],'FaceColor',[.95 .45 .1],'EdgeColor',[.95 .45 .1]); % Distribution as the +/- 2SD
        set(gca,'tickdir','out', ...
            'box','off', ...
            'ylim',[0 0.6], ...
            'xlim',[min_x,max_x], ...
            'ytick',[0 0.25 0.5], ...
            'xtick',round(linspace(min_x,max_x,4)), ...
            'fontsize',16)
        ylabel('Probability','fontsize',16)
        xlabel('rmse','fontsize',16')
        title(sprintf('Strength of connection evidence %2.3f',(se.s.mean)), ...
            'FontSize',16)    
        saveFig(fig(1).fh,fullfile(saveDir, fig(1).name),1)
    end
    
    % RMSE scatter-density plot of  Probabilistic and Deterministic
    figNameRmse = sprintf('prob_vs_det_rmse_common_voxels_ma_%s',fname);
    fh = scatterPlotRMSE(det,prob,figNameRmse);
    saveFig(fh,fullfile(saveDir, figNameRmse),1)
    saveFig(fh,fullfile(saveDir, figNameRmse),0)

    if display.evidence
        % Binned RMSE distributions
        fig(3).name = sprintf('Size_of_effect_of_the_lesion_%s',mfilename);
        fig(3).fh   = figure('name',fig(3).name,'color','w');
        set(fig(3).fh,'Units','normalized','Position',[0.007 0.55  0.28 0.36]);
        subplot(1,4,1)
        plot(1,se.s.mean,'-o','color', [.95 .45 .1],'linewidth',2); hold on
        plot([1,1], [se.s.mean,se.s.mean] + [-se.s.std,se.s.std], '-','color',[.95 .45 .1] )
        ylabel('S (s.d.)', 'fontsize',14);
        set(gca,'box','off','xlim',[0 2], 'ylim',[0 ceil(se.s.mean + se.s.std)], ...
            'tickdir', 'out', 'ticklength', [0.025 0])
        subplot(1,4,2)
        plot(1,se.em.mean,'-o','color', [.95 .45 .1],'linewidth',2); hold on
        ylabel('Earth mover''s distance', 'fontsize',14);
        set(gca,'box','off','xlim',[0 2], 'ylim',[0 ceil(se.em.mean)], ...
            'tickdir', 'out', 'ticklength', [0.025 0])
        subplot(1,4,3)
        plot(1,se.kl.mean,'-o','color', [.95 .45 .1],'linewidth',2); hold on
        ylabel('K-L divergence (bits)', 'fontsize',14);
        set(gca,'box','off','xlim',[0 2], 'ylim',[0 ceil(se.kl.mean)], ...
            'tickdir', 'out', 'ticklength', [0.025 0])
        subplot(1,4,4)
        plot(1,se.j.mean,'-o','color', [.95 .45 .1],'linewidth',2); hold on
        ylabel('Jeffrey''s divergence (bits)', 'fontsize',14);
        set(gca,'box','off','xlim',[0 2], 'ylim',[0 ceil(se.j.mean)], ...
            'tickdir', 'out', 'ticklength', [0.025 0])    
        saveFig(fig(2).fh,fullfile(saveDir, fig(2).name),1)
        saveFig(fig(3).fh,fullfile(saveDir, fig(3).name),1)
    end

    close all
    save(fullfile(saveDir,sprintf('average_results_%s.mat',trackingType)),'se')
end


end

%------------------------------------%
function fhRmseMap = scatterPlotRMSE(det,prob,figNameRmse)
fhRmseMap = mrvNewGraphWin(figNameRmse);
[ymap,x]  = hist3([det.rmse;prob.rmse]',{[10:1:90], [10:1:90]});
ymap = ymap./length(prob.rmse);
sh   = imagesc(flipud(log10(ymap)));
cm   = colormap(flipud(hot)); view(0,90);
axis('square')      
set(gca, ...
    'xlim',[1 length(x{1})],...
    'ylim',[1 length(x{1})], ...
    'ytick',[1 (length(x{1})/2) length(x{1})], ...
    'xtick',[1 (length(x{1})/2) length(x{1})], ...
    'yticklabel',[x{1}(end) x{1}(round(end/2)) x{1}(1)], ...
    'xticklabel',[x{1}(1)   x{1}(round(end/2)) x{1}(end)], ...
    'tickdir','out','ticklen',[.025 .05],'box','off', ...
    'fontsize',16,'visible','on')
hold on
plot3([1 length(x{1})],[length(x{1}) 1],[max(ymap(:)) max(ymap(:))],'k-','linewidth',1)
ylabel('Deterministic_{rmse}','fontsize',12)
xlabel('Probabilistic_{rmse}','fontsize',12)
cb = colorbar;
tck = get(cb,'ytick');
set(cb,'yTick',[min(tck)  mean(tck) max(tck)], ... 'yTickLabel',10.^[min(tck), mean(tck), max(tck)], ...
    'tickdir','out','ticklen',[.025 .05],'box','on', ...
    'fontsize',16,'visible','on')
end

%---------------------------------------%
function fh = distributionPlotRMSE(y_e,ywo_e,dprime,xhis,woxhis,figName)

h1.ylim  = [0 0.6];
h1.xlim  = [22,34];
h1.ytick = [0 0.3 0.6];
h1.xtick = [28 30 32 34];
h2.ylim  = [0 0.4];
h2.xlim  = [28,32];
h2.ytick = [0 0.2 0.4];
h2.xtick = [28 30 32];
histcolor{1} = [0 0 0];
histcolor{2} = [.95 .6 .5];

fh = mrvNewGraphWin(figName);
patch([xhis,xhis],y_e(:),histcolor{1},'FaceColor',histcolor{1},'EdgeColor',histcolor{1});
hold on
patch([woxhis,woxhis],ywo_e(:),histcolor{2},'FaceColor',histcolor{2},'EdgeColor',histcolor{2}); 
set(gca,'tickdir','out', ...
        'box','off', ...
        'ylim',[0 .6], ... 
        'xlim',[28 34], ...
        'xtick',[28 30 32 34], ...
        'ytick',[0 .3 .6], ...
        'fontsize',16)
ylabel('Probability','fontsize',16)
xlabel('rmse','fontsize',16')

title(sprintf('Strength of evidence:\n mean %2.3f - std %2.3f',mean(dprime),std(dprime)), ...
    'FontSize',16)
legend({'Probabilistic','Deterministic'})
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