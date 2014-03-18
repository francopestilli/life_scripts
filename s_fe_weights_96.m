function s_fe_weights_96(bval,tractographyType)
%
% This function:
%  - Load a series of results from precomputed connectomes
%  - Extract the weights and fiber lenths.
%  - Make a plot of the distribution of weights and fiber lengths.
%
%  s_fe_weights(bval,tractographyType)
% 
% Copyright Franco Pestilli (2014) Stanford University

% Get the base directory for the data
datapath = '/home/frk/2t1/predator/';
subjects = {...
            'KK_96dirs_b2000_1p5iso', ...           
            'MP_96dirs_b2000_1p5iso', ...
            'JW_96dirs_b2000_1p5iso', ...
            'HT_96dirs_b2000_1p5iso', ...
            'KW_96dirs_b2000_1p5iso', ...
            'FP_96dirs_b2000_1p5iso', ...
            };
        
if notDefined('tractographyType'), tractographyType = 'lmax10'; end
if notDefined('bval'), bval = []; end
if notDefined('saveDir'), saveDir = fullfile('/home/frk/Dropbox','pestilli_etal_revision',mfilename);end

% These are default parameters for plotting
nBins = 11; % first bin is zero weight second is
bins_w = linspace(-10,-1,nBins);
proportionDeleted = nan(size(subjects,2),2);
w = cell(1,length(subjects));
len = w;
for isbj = 1:length(subjects)
    % Directory where to load the results
    resultsPath  = fullfile(datapath,subjects{isbj},'results');
    savedir      = fullfile(saveDir,subjects{isbj});

    if ~isempty(bval)
       resFiles       = dir(fullfile(resultsPath,sprintf('*%s*.mat',num2str(bval))));
    else
       resFiles       = dir(fullfile(resultsPath,sprintf('*%s*.mat',tractographyType)));
    end
            
    % We build one modelper fiber group, whole brain fiber group
    for iRes = 1:length(resFiles)
        % The final connectome and dat astructure will be saved with this name:
        [~,resFileName] = fileparts(resFiles(iRes).name);
        
        % Buil a full-file of the fibers and the RESULTS structure to load
        resFileName2Load = fullfile(resultsPath,[resFileName,'.mat']);

        % Initialize the Connectome
        fprintf('[%s] Loading results: \n%s\n',mfilename,resFileName2Load)
        load(resFileName2Load);
        
        % Reorganize the weights and fiber-lengths in structures divided by
        % subjects
        w{isbj}   = results.weights(results.weights > 0);
        olen{isbj} = results.optimized.length;
        clen{isbj} = results.candidate.length;

        % Distribution of length
        [f, ol(isbj,:), cl(isbj,:), lbins] = plotLengthHist(tractographyType,olen{isbj},clen{isbj});
        feSavefig(f,'verbose','yes','figName',[resFileName, '_optimized_fiberlength_hist'],'figDir',savedir,'figType','eps');
        close(f)
            
        % Distribution of weights
        [f, ow(isbj,:), wbins] = plotWeightsHist(tractographyType,w{isbj});
        feSavefig(f,'verbose','yes','figName',[resFileName, '_optimized_weights_hist'],'figDir',savedir,'figType','eps');
        close(f)

        % 2d histogram
        %[f, h] = plotLengthWeightScatter(tractographyType,w{isbj},len{isbj});
        %feSavefig(f,'verbose','yes','figName',[resFileName, '_weight_vs_fiberlength_2d'],'figDir',savedir,'figType','eps');
        %set(gca, 'visible', 'off')
        %feSavefig(f,'verbose','yes','figName',[resFileName, '_weight_vs_fiberlength_2d'],'figDir',savedir,'figType','jpg');
        %close(f)

        % Just reorganize the vectors
        %len{isbj} = len{isbj}';
        %w{isbj}   = w{isbj}';

        % Save the proportion of zero-weight and non-zero-weight fibers
        proportionDeleted(isbj,:) = [results.candidate.n - results.optimized.n  results.optimized.n]./results.candidate.n;
        [f, h] = plotFiberPie(resFileName,proportionDeleted(isbj,:));
        feSavefig(f,'verbose','yes','figName',[resFileName, '_pie'],'figDir',savedir,'figType','eps');
        close(f)

     end
end

mow = mean(ow,1);
owerr = [mow;mow]+[-std(ow)./sqrt(size(ow,1));std(ow)./sqrt(size(ow,1))];
f = plotMeanWeightsHist(tractographyType,mow,owerr,wbins);
feSavefig(f,'verbose','yes','figName',[tractographyType, '_weight_hist_across_subjects'],'figDir',fullfile(saveDir,'averages'),'figType','eps');
close(f)

mol = mean(ol,1);
olerr = [mol;mol]+[-std(ol)./sqrt(size(ow,1));std(ol)./sqrt(size(ow,1))];
mcl = mean(cl,1);
clerr = [mcl;mcl]+[-std(cl)./sqrt(size(ow,1));std(cl)./sqrt(size(ow,1))];
f = plotMeanLengthHist(tractographyType,mol,olerr,mcl,clerr,lbins);
feSavefig(f,'verbose','yes','figName',[tractographyType, '_length_hist_across_subjects'],'figDir',fullfile(saveDir,'averages'),'figType','eps');
close(f)

ratioLength = ol./cl;
drl = mean(ratioLength,1);
drlerr = [drl;drl]+[-std(ratioLength)./sqrt(size(ratioLength,1));std(ratioLength)./sqrt(size(ratioLength,1))];
f = plotMeanLengthRatioHist(tractographyType,drl,drlerr,lbins);
feSavefig(f,'verbose','yes','figName',[tractographyType, '_length_ratio_hist_across_subjects'],'figDir',fullfile(saveDir,'averages'),'figType','eps');
close(f)

[f, h] = plotProportionDeleted(tractographyType,proportionDeleted);
feSavefig(f,'verbose','yes','figName',[tractographyType, '_deleted_fibers_across_subjects'],'figDir',fullfile(saveDir,'averages'),'figType','eps');
close(f)

% [f, h] = plotLengthWeightScatter(['_ACROSS_FIVE_SUBJECTS_',tractographyType],cell2mat(w)',cell2mat(len)');
% feSavefig(f,'verbose','yes','figName',[tractographyType, '_weight_vs_fiberlength_2d_across_subjects'],'figDir','~/2t1/predator/average_figures','figType','eps');
% set(gca, 'visible', 'off')
% feSavefig(f,'verbose','yes','figName',[tractographyType, '_weight_vs_fiberlength_2d_across_subjects'],'figDir','~/2t1/predator/average_figures','figType','jpg');
% close(f)
% 
end % End main function

function f = plotMeanLengthRatioHist(tractographyType,drl,drlerr,bins)
fontSiz = 15;
% Make a 2D histogram (scatter plot) of the percent deleted and kept fibers
f = figure('name',sprintf('mean_length_ratio_hist_%s',tractographyType),'color','w');
sh = plot([bins;bins],drlerr,'-','color','r','linewidth',2);
hold on
sh = plot(bins,(drl),'-','color','k','markeredgecolor','w','markerfacecolor','k','markersize',12);
ylabel('Proportion deleted fascicles','fontsize',fontSiz)
xlabel('Fascicle length','fontsize',fontSiz)
set(gca,'xlim',[0 200],...
    'xtick',[0 100 200],...
    'ylim',    [0.25 1], ...
    'ytick',[0.25 0.5 1], ...
    'yscale','log', ...
    'tickdir','out','box','off', ...
    'fontsize',fontSiz,'visible','on')

end

function f = plotMeanLengthHist(tractographyType,ol,olerr,cl,clerr,bins)
fontSiz = 15;
% Make a 2D histogram (scatter plot) of the percent deleted and kept fibers
f = figure('name',sprintf('mean_length_hist_%s',tractographyType),'color','w');
sh = plot([bins;bins],olerr,'-','color','r','linewidth',2);
hold on
sh = plot([bins;bins],clerr,'-','color','r','linewidth',2);
sh = plot(bins,(ol),'-','color','k','markeredgecolor','w','markerfacecolor','k','markersize',12);
sh = plot(bins,(cl),'-','color',[.6 .6 .6],'markeredgecolor','w','markerfacecolor',[.6 .6 .6],'markersize',12);
ylabel('Number of fascicles','fontsize',fontSiz)
xlabel('Fascicle length','fontsize',fontSiz)
set(gca,'xlim',[0 200],...
    'xtick',[0 100 200],...
    'ylim',    [0 2.4*10^5], ...
    'ytick',[0 1.2*10^5 2.4*10^5], ...
    'tickdir','out','box','off', ...
    'fontsize',fontSiz,'visible','on')

end

function f = plotMeanWeightsHist(tractographyType,ow,owerr,bins)
fontSiz = 15;
% Make a 2D histogram (scatter plot) of the percent deleted and kept fibers
f = figure('name',sprintf('mean_weights_hist_%s',tractographyType),'color','w');
sh = plot([bins;bins],owerr,'-','color','r','linewidth',2);
hold on
sh = plot(bins,(ow),'-','color','k','markeredgecolor','w','markerfacecolor','k','markersize',12);
ylabel('Number of fascicles','fontsize',fontSiz)
xlabel('log_{10}(Fascicle weight)','fontsize',fontSiz)
set(gca,'xlim',[-6 0],...
    'xtick',[-6 -3 0],...
    'ylim',    [0 80000], ...
    'ytick',[0 40000 80000], ...
    'tickdir','out','box','off', ...
    'fontsize',fontSiz,'visible','on')

end

function [f, oy,bins] = plotWeightsHist(tractographyType,ow)
fontSiz = 15;
% Make a 2D histogram (scatter plot) of the percent deleted and kept fibers
f = figure('name',sprintf('weights_hist_%s',tractographyType),'color','w');
bins = linspace(-6,0,22);
[oy,x]  = hist(log10(ow),bins);
sh = plot(bins,oy,'-','color','k','markeredgecolor','w','markerfacecolor','k','markersize',12);
ylabel('Number of fascicles','fontsize',fontSiz)
xlabel('log_{10}(Fiber weight)','fontsize',fontSiz)
set(gca,'xlim',[-6 0],...
    'xtick',[-6 -3 0],...
    'ylim',    [0 80000], ...
    'ytick',[0 40000 80000], ...
    'tickdir','out','box','off', ...
    'fontsize',fontSiz,'visible','on')

end

function [f, oy, cy, bins] = plotLengthHist(tractographyType,olen,clen)
fontSiz = 15;
% Make a 2D histogram (scatter plot) of the percent deleted and kept fibers
f = figure('name',sprintf('fiberlength_hist_%s',tractographyType),'color','w');
bins = [10:36:200];
[oy,x]  = hist(olen,bins);
[cy,x]  = hist(clen,bins);
sh = plot(bins,(cy),'-','color',[.6 .6 .6],'markeredgecolor','w','markerfacecolor',[.6 .6 .6],'markersize',12);
hold on
sh = plot(bins,(oy),'-','color','k','markeredgecolor','w','markerfacecolor','k','markersize',12);
ylabel('Number of fascicles','fontsize',fontSiz)
xlabel('Fascicle length (mm)','fontsize',fontSiz)
set(gca,'xlim',[0 200],...
    'xtick',[10 100 200],...
    'ylim',    [0 2.4*10^5], ...
    'ytick',[0 1.2*10^5 2.4*10^5], ...
    'tickdir','out','box','off', ...
    'fontsize',fontSiz,'visible','on')

end

function [f, ymap] = plotLengthWeightScatter(tractographyType,w,len)
fontSiz = 15;
% Make a 2D histogram (scatter plot) of the percent deleted and kept fibers
f = figure('name',sprintf('weight_vs_fiberlength_%s',tractographyType),'color','w');
nBins = 46;
ybins = linspace(0,nBins,nBins);
ybins = logspace(log10(5),log10(400),nBins);
xbins = linspace(-9,-1,nBins);
[ymap,x]  = hist3([(len), log10(w)],{ ybins, xbins});
ymap = ymap./numel(len);
sh = imagesc(flipud(ymap));
cm = colormap(flipud(hot)); 
view(0,90); axis('square');
xlabel('log_1_0(Fascicle weight)','fontsize',fontSiz)
ylabel('Fascicle length (mm)','fontsize',fontSiz)

set(gca,'xlim',[nBins/2 nBins],...
    'ylim',    [5 nBins], ...
    'ytick',[5 nBins/2 nBins], ...
    'yticklabel',[round(ybins(nBins)) round(ybins(round(mean([nBins, nBins/2])))) round(ybins(nBins/2)) ], ...
    'xtick',[nBins/2 mean([nBins, nBins/2]) nBins], ...
    'xticklabel',[round(xbins(nBins/2)) round(xbins(round(mean([nBins, nBins/2])))) round(xbins(nBins)) ], ...
    'tickdir','out','box','off', ...
    'fontsize',fontSiz,'visible','on')

end

function [f, h] = plotProportionDeleted(tractographyType,proportionDeleted) 
% Make aplot of the percent deleted and kept fibers
fontSiz = 15;
f = figure('name',sprintf('mean_proportion_deleted_fibers_%s',tractographyType),'color','w');
m  = mean(proportionDeleted,1);
sd = [m - std(proportionDeleted)./sqrt(size(proportionDeleted,1)); m + std(proportionDeleted./sqrt(size(proportionDeleted,1)))];
h(1) = bar(m,'facecolor','k'); hold on
plot([1 2; 1 2],sd,'r-','linewidth',2)
ylabel('Proportion','fontsize',fontSiz)
set(gca,'tickdir','out','color','w','box','off', ...
    'ylim',[0 .75], ...
    'ytick',[0 .25 .5 .75], ...
    'xlim',[0 3], 'fontsize',fontSiz,...
    'xticklabel',{'Deleted','Kept'})
end

function [f, h] = plotFiberPie(resFileName,proportionDeleted)
% Make a pie graphshwing the percent of values below and above zero
fontSiz = 15;
f = figure('name',sprintf('pie_of_deleted_fibers_%s',resFileName),'color','w');
h = pie( proportionDeleted,[1 0]);
colormap gray
patchObjs = findobj(h,'Type','patch');
textObjs  = findobj(h,'Type','text');
oldStr    = get(textObjs,{'String'});
val    = get(textObjs,{'Extent'});
oldExt = cat(1,val{:});
newStr = {sprintf('%s deleted',oldStr{1});sprintf('%s kept',oldStr{2})};
set(textObjs,{'String'},newStr)
end

        
% load run01_fliprot_aligned_trilin_csd_lmax10_run01_fliprot_aligned_trilin_brainmask_run01_fliprot_aligned_trilin_wm_prob-800000-fiberStatsResults.mat
% w = results.weights(results.weights>0);
% [f,x_vals] = ecdf(log10(w));
% F = plot(x_vals,f);
% hold on
% G = plot(x_vals,normcdf(x_vals,mean(log10(w)),std(log10(w))),'r-');
% [H, pValue, KSstatistic, criticalValue] = kstest(f)