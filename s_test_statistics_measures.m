function s_test_statistics_measures
% Test the variability of S and K-L or Jeffrey as function of:
% - Effect size, namely he difference between the distribution of errors
% with and without the lesion
% - number of voxels, the numer of observations sustaining the hypothesis
%
% Copyright Franco Pestilli Stanford University 2014
if notDefined('saveDir'), savedir = fullfile('/marcovaldo/frk/Dropbox','pestilli_etal_revision',mfilename);end

nBoots = 500;
nvox= [ 2^13 2^12 2^11 2^10 2^9 2^7 2^6 2^5 2^4];
nboots = 10000; 
nmontecarlo = 4;

% We set up two idealized RMSE distributions one for the lesioned and one for the ulesioned connectome
% We change the number of voxels in each distribution and measure the
% resulting indices (KL, J and S).
%
% The we repeate the same operation for a different problem. Tis time the
% distributions have a different effect size (meaning lesioning has a
% different impact on the rmse.

% RMSE distributions ar emodelled as Rician distributions.
% Lesioned RMSE distrbution
vu = 80;% 20 
su = 100;% 20

% Unlsioned RMSE distrbution
vl = 80;% 100
sl = [300, 250, 200, 150, 100];% Try s = 20 vs 100 for an example of how two similar means generate very different K-

colors = copper(8);
colors = colors(4:end,:);

for is = 1:length(sl)    
for invx=1:length(nvox)
    fprintf('[%s] nvoxel: %i\n',mfilename,nvox(invx))
    for ibt = 1:nBoots 
        fprintf('[%s] boot num %i\n',mfilename,ibt)
        % Simulate a couple of RMSE distributions
        unlesioned_rmseall = ricernd(vu*ones(1, nvox(invx)), su);
        lesioned_rmseall   = ricernd(vl*ones(1, nvox(invx)), sl(is));
        urmse(ibt,is) = mean(unlesioned_rmseall);
        lrmse(ibt,is) = mean(lesioned_rmseall);

        se.xrange(1) = 0;  
        se.xrange(2) = 1000;
        se.nbins     = 60;
        se.bins      = linspace(se.xrange(1),se.xrange(2),se.nbins);
        [se.lesion.hist, se.lesion.xhist] = hist(lesioned_rmseall,se.bins);
        se.lesion.hist = se.lesion.hist ./ sum(se.lesion.hist);
        [se.nolesion.hist, se.nolesion.xhist] = hist(unlesioned_rmseall,se.bins);
        se.nolesion.hist   = se.nolesion.hist ./ sum(se.nolesion.hist);
        
        % KL
        tmp = se.nolesion.hist .* log2( (se.nolesion.hist) ./ (se.lesion.hist + eps) );
        tmp(isnan(tmp)) = 0;
        se.value(ibt) = nansum(tmp);clear tmp
        
        % Jeffrey
        tmp = se.nolesion.hist .* log2( (se.nolesion.hist) ./ ((se.lesion.hist + se.lesion.hist + eps)./2)  ) + ...
              se.lesion.hist .* log2( (se.lesion.hist) ./ ((se.lesion.hist + se.lesion.hist + eps)./2) );
        tmp(isnan(tmp)) = 0;
        se.jvalue(ibt) = nansum(tmp);clear tmp     
        
        % EMD
        [~,tmp_em(ibt)] = emd(se.nolesion.xhist',se.lesion.xhist',se.nolesion.hist',se.lesion.hist',@gdf);    
    end
    
    se.em.name = sprintf('Earth Mover''s distance: http://en.wikipedia.org/wiki/Earth_mover''s_distance');
    se.em.mean(invx,is) = mean(tmp_em);
    se.em.std(invx,is)  = std(tmp_em);
    
    se.kl.name = sprintf('Kullback–Leibler divergence: http://en.wikipedia.org/wiki/Kullback-Leibler_divergence');
    se.kl.mean(invx,is) = mean(se.value);
    se.kl.std(invx,is)  = std(se.value);
    
    se.j.name = sprintf('Jeffrey''s divergence: http://en.wikipedia.org/wiki/Divergence_(statistics)');
    se.j.mean(invx,is)  = mean(se.jvalue);
    se.j.std(invx,is)   = std(se.jvalue);
    
    s = compute_s(lesioned_rmseall, unlesioned_rmseall,nboots,nmontecarlo,[]);
    se.s.name = sprintf('strength of evidence, d-prime: http://en.wikipedia.org/wiki/Effect_size');
    se.s.mean(invx,is) = s.mean; 
    se.s.std(invx,is)  = s.std; 

    if invx == 1
        % Make a plot
        fh(1) = figure('name',sprintf('%s_distributions_effect_size%i',mfilename, is),'color','w');
        plot(se.lesion.xhist,se.lesion.hist,'-','color',colors(is,:),'linewidth',2); hold on
        plot(se.nolesion.xhist,se.nolesion.hist,'k-','linewidth',2)
        title(sprintf('mean RMSE \n no-lesion %2.3f\ | lesion %2.2f',mean(unlesioned_rmseall),mean(lesioned_rmseall)),'fontsize',16)
        ylabel('Probability', 'fontsize',14);xlabel('RMSE', 'fontsize',14)
        if is == 1, legend({'Lesion','nvoxolesion'},'fontsize',14);end
        axis square
        set(gca,'box','off','xtick',[0 500 1000],'ytick',[0 .06 .12],'xlim',[0 1000],'ylim',[0 .125], ...
            'tickdir', 'out', 'ticklength', [0.025 0])
        drawnow
        saveFig(fh(1),fullfile(savedir, sprintf('%s_distributions_effect_size%i',mfilename, is)),1)
    end

end
end

fh(2) = figure('name',sprintf('%s_S',mfilename),'color','w');
hold on
for is = 1:size(se.s.mean,2)
    plot(nvox',se.s.mean(:,is),'-','color',colors(is,:),'linewidth',2);
    hold on
    plot([nvox',nvox']',[se.s.mean(:,is),se.s.mean(:,is)]' + [-se.s.std(:,is),se.s.std(:,is)]','-','color',colors(is,:),'linewidth',2);
end
ylabel('S', 'fontsize',14);xlabel('number of voxels', 'fontsize',14)
set(gca,'box','off', 'fontsize',14,'tickdir','out', 'ticklength', [0.025 0],'xscale','log','xtick',fliplr(nvox))
axis square
saveFig(fh(2),fullfile(savedir, sprintf('%s_S',mfilename)),1)

fh(3) = figure('name',sprintf('%s_KL',mfilename),'color','w');
hold on
for is = 1:size(se.s.mean,2)
    plot(nvox',se.kl.mean(:,is),'-','color',colors(is,:),'linewidth',2);
    hold on
    plot([nvox',nvox']',[se.kl.mean(:,is),se.kl.mean(:,is)]' + [-se.kl.std(:,is),se.kl.std(:,is)]','-','color',colors(is,:),'linewidth',2);
end
ylabel('K-L', 'fontsize',14);xlabel('number of voxels', 'fontsize',14)
set(gca,'box','off', 'fontsize',14,'tickdir','out', 'ticklength', [0.025 0],'xscale','log','xtick',fliplr(nvox))
axis square
saveFig(fh(3),fullfile(savedir, sprintf('%s_KL',mfilename)),1)

fh(4) = figure('name',sprintf('%s_JEFFREY',mfilename),'color','w');
hold on
for is = 1:size(se.s.mean,2)
    plot(nvox',se.j.mean(:,is),'-','color',colors(is,:),'linewidth',2);
    hold on
    plot([nvox',nvox']',[se.j.mean(:,is),se.j.mean(:,is)]' + [-se.j.std(:,is),se.j.std(:,is)]','-','color',colors(is,:),'linewidth',2);
end
ylabel('Jeffrey', 'fontsize',14);xlabel('number of voxels', 'fontsize',14)
set(gca,'box','off', 'fontsize',14,'tickdir','out', 'ticklength', [0.025 0],'xscale','log','xtick',fliplr(nvox))
axis square
saveFig(fh(4),fullfile(savedir, sprintf('%s_JEFFREY',mfilename)),1)

fh(5) = figure('name',sprintf('%s_EMD',mfilename),'color','w');
hold on
for is = 1:size(se.s.mean,2)
    plot(nvox',se.em.mean(:,is),'-','color',colors(is,:),'linewidth',2);
    hold on
    plot([nvox',nvox']',[se.em.mean(:,is),se.em.mean(:,is)]' + [-se.em.std(:,is),se.em.std(:,is)]','-','color',colors(is,:),'linewidth',2);
end
ylabel('EMD', 'fontsize',14);xlabel('number of voxels', 'fontsize',14)
set(gca,'box','off', 'fontsize',14,'tickdir','out', 'ticklength', [0.025 0],'xscale','log','xtick',fliplr(nvox))
axis square
saveFig(fh(5),fullfile(savedir, sprintf('%s_EMD',mfilename)),1)

fh(6) = figure('name',sprintf('%s_COLORBAR',mfilename),'color','w');
colormap(colors)
colorbar('ytick',[1:size(colors,1)]+.5,'yticklabel',round([mean(lrmse,1) - mean(urmse,1) ]))
set(gca, 'fontsize',14,'tickdir','out', 'ticklength', [0.025 0])
saveFig(fh(6),fullfile(savedir, sprintf('%s_COLORBAR',mfilename)),1)

% Save the results
save(fullfile(savedir, 'results_1'),'se')

keyboard
end % Main function

function s = compute_s(lesioned_rmseall, unlesioned_rmseall,nboots,nmontecarlo,sbpl)
disp('computing S ...')
% The following is the code for the bootstrap test on the MEAN rmse
sizeunlesioned    = length(unlesioned_rmseall);
nullDistributionW = nan(nboots,nmontecarlo);
nullDistributionWO = nan(nboots,nmontecarlo);
min_x = floor(mean([unlesioned_rmseall]) - mean([unlesioned_rmseall])*.05);
max_x =  ceil(mean([lesioned_rmseall])   + mean([lesioned_rmseall])*.05);
min_x = min([min_x,max_x]);
max_x = max([max_x,min_x]);
for inm = 1:nmontecarlo
    parfor ibt = 1:nboots
        nullDistributionW(ibt,inm) = mean(randsample(unlesioned_rmseall,   sizeunlesioned,true));      
        nullDistributionWO(ibt,inm) = mean(randsample(lesioned_rmseall,sizeunlesioned,true));
    end
    
    % Distribution unlesioned
    [y(:,inm),xhis] = hist(nullDistributionW(:,inm),linspace(min_x,max_x,200));
    y(:,inm) = y(:,inm)./sum(y(:,inm));
    
    % Distribution lesioned
    [woy(:,inm),woxhis] = hist(nullDistributionWO(:,inm),linspace(min_x,max_x,200));
    woy(:,inm) = woy(:,inm)./sum(woy(:,inm));
end
y_m = mean(y,2);
y_e = [y_m, y_m] + 2*[-std(y,[],2),std(y,[],2)];

ywo_m = mean(woy,2);
ywo_e = [ywo_m, ywo_m] + 2*[-std(woy,[],2),std(woy,[],2)];

if ~isempty(sbpl)
% Plot the null distribution and the empirical difference
subplot(2,3,sbpl);
patch([xhis,xhis],y_e(:),'b','FaceColor',[.67 .86 .96],'EdgeColor',[.67 .86 .96]); % Distribution as the +/- 2SD
hold on
patch([woxhis,woxhis],ywo_e(:),[.97 .66 .76],'FaceColor',[.97 .66 .76],'EdgeColor',[.97 .66 .76]); % Distribution as the +/- 2SD
set(gca,'tickdir','out', ...
        'box','off', ...
        'ylim',[ 0.2], ... 
        'xlim',[min_x,max_x], ...
        'ytick',[0 0.1 0.2], ...
        'xtick',round(linspace(min_x,max_x,4)), ...
        'fontsize',16)
ylabel('Probability','fontsize',16)
xlabel('rmse','fontsize',16')
end

% (3) Compute the probability that the empirical difference (1) was
%     observed by chance given th data, by looking at the percentile of the
%     empirical difference in the Nul distribution (2).
s.mean = mean(diff([mean(nullDistributionW,1); ...
                          mean(nullDistributionWO,1)])./sqrt(sum([std(nullDistributionW,[],1);std(nullDistributionWO,[],1)].^2,1)));
s.std  = std(diff([mean(nullDistributionW,1); ...
                          mean(nullDistributionWO,1)])./sqrt(sum([std(nullDistributionW,[],1);std(nullDistributionWO,[],1)].^2,1)));

if ~isempty(sbpl)
title(sprintf('S %2.3f',(s)), 'FontSize',16);
end

end

function r = ricernd(v, s)
%RICERND Random samples from the Rice/Rician probability distribution.
%   r = ricernd(v, s) returns random sample(s) from the Rice (aka Rician) 
%   distribution with parameters v and s.
%   (either v or s may be arrays, if both are, they must match in size)
%
%   R ~ Rice(v, s) if R = sqrt(X^2 + Y^2), where X ~ nvox(v*cos(a), s^2) and
%   Y ~ nvox(v*sin(a), s^2) are independent normal distributions (any real a).
%   nvoxote that v and s are *not* the mean and standard deviation of R!
%
%   The size of Y is the common size of the input arguments.  A scalar
%   input functions as a constant matrix of the same size as the other
%   inputs.
%
%   nvoxote, to add Rician noise to data, with given s and data-dependent v:
%     new = ricernd(old, s);
%
%   Reference: http://en.wikipedia.org/wiki/Rice_distribution (!)
%
%   Example:
%
%     % Compare histogram of random samples with theoretical PDF:
%     v = 4; s = 3; nvox= 1000;
%     r = ricernd(v*ones(1, nvox), s);
%     c = linspace(0, ceil(max(r)), 20);
%     w = c(2); % histogram bin-width
%     h = histc(r, c); bar(c, h, 'histc'); hold on
%     xl = xlim; x = linspace(xl(1), xl(2), 100);
%     plot(x, nvox*w*ricepdf(x, v, s), 'r');
%     
%   See also RICEPDF, RICESTAT, RICEFIT

%   Missing (?) 'See also's RICECDF, RICEInvoxV, RICELIKE

%   Inspired by normpdf from the MATLAB statistics toolbox
%   Copyright 2008 Ged Ridgway (Ged at cantab dot net)

if isscalar(v)
    dim = size(s);
elseif isscalar(s)
    dim = size(v);
elseif all(isequal(size(v), size(s)))
    % (both non-scalar, matching)
    dim = size(v); % == size(s)
else
    error('ricernd:InputSizeMismatch','Sizes of s and v inconsistent.')
end

x = s .* randn(dim) + v;
y = s .* randn(dim);
r = sqrt(x.^2 + y.^2);

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
