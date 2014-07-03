function s_fe_mt_ips_tract_detect_plot()
%
% This script performs a test of conenctivity of MT+ (Zilles t al ROI from 
% freesurfer) with the Superior parietal Cortex (Aparc Freesurfer segementation). 
% The following are the steps we perform:
%  - It loads a whole-brain tractography solution. 
%  - It loads the MTand Parietal ROI. 
%  - It finds the tract connecting the two ROIs.
%  - It builds a connectome model in the ROI defined by the tract
%  - It performs a series of virtual lesions on fibers with different weigths.
%    we compute the quartiles of the fibers by weight and leasion different number of them by selecting them within each quartile.
%    This test is to show .
%
% Copyright by Franco Pestilli Stanford University, 2014

% Handle parallel computing
if matlabpool('size') == 0
    c = parcluster;
    c.NumWorkers = 12;
    matlabpool(c);
end


if notDefined('saveDir'), saveDir = fullfile('/marcovaldo/frk/Dropbox','pestilli_etal_revision','s_fe_mt_ips_tract_detect',mfilename);end
            
tic, fprintf('\n[%s] loading results of s_fe_mt_ips_tract_detect... \n',mfilename)
load(fullfile('/marcovaldo/frk/Dropbox','pestilli_etal_revision','s_fe_mt_ips_tract_detect','strength_of_evidence_quartiles.mat'),'SE','w','cull'); toc

% Reorganzie the results of the virtual lesions: SE(iSbj,ih,iqr, ifib)
se = reshape(SE,10,2,4);
for isbj = 1:size(se,1)
    for iqrt = 1:size(se,2)
        for infbs = 1:size(se,3)    
            res.s(isbj,iqrt,infbs) = se(isbj,iqrt,infbs).s.mean;
            res.e(isbj,iqrt,infbs) = se(isbj,iqrt,infbs).em.mean;
            res.j(isbj,iqrt,infbs) = se(isbj,iqrt,infbs).j.mean;
            res.k(isbj,iqrt,infbs) = se(isbj,iqrt,infbs).kl.mean;
        end
    end
end

m.s = squeeze(median(res.s,  1));
sd.s = squeeze(std(res.s,[],1) ./ sqrt(size(res.s,1)));
sd.s(sd.s > 20) = 2;

m.e = squeeze(median(res.e,  1));
sd.e = squeeze(std(res.e,[],1) ./ sqrt(size(res.e,1)));

m.j = squeeze(median(res.j,  1));
sd.j = squeeze(std(res.j,[],1) ./ sqrt(size(res.j,1)));

m.k = squeeze(median(res.k,  1));
sd.k = squeeze(std(res.k,[],1) ./ sqrt(size(res.e,1)));

% Plot the results.
close all
whichquartiles = [1,2,3];
colors ={'r','g','b'};
x = [1,2,4,8,16,32,64];

figName = sprintf('s_fe_mt_ips_deletion_tests_strength_of_evidence_across_subjects');
fh =figure('name',figName,'color','w');
for ip = whichquartiles
    semilogx(x,m.s(ip,:),'-','color',colors{ip})
    hold on
    semilogx([x;x],[m.s(ip,:);m.s(ip,:)] + [-sd.s(ip,:);sd.s(ip,:)],'-','color',colors{ip})
end
ylabel('EMD');xlabel('Number of lesioned fibers')
set(gca,'box','off','tickdir','out','xtick',x)
saveFig(fh,fullfile(saveDir,figName),'eps')

figName = sprintf('s_fe_mt_ips_deletion_tests_EMD_across_subjects');
fh =figure('name',figName,'color','w');
for ip = whichquartiles
    semilogx(x,m.e(ip,:),'-','color',colors{ip});
    hold on
    semilogx([x;x],[m.e(ip,:);m.e(ip,:)] + [-sd.e(ip,:);sd.e(ip,:)],'-','color',colors{ip})
end
ylabel('EMD');xlabel('Number of lesioned fibers')
set(gca,'box','off','tickdir','out','xtick',x)
saveFig(fh,fullfile(saveDir,figName),'eps')

figName = sprintf('s_fe_mt_ips_deletion_tests_jeffrey_across_subjects');
fh =figure('name',figName,'color','w');
for ip = whichquartiles
    semilogx(x,m.j(ip,:),'-','color',colors{ip});
    hold on
    semilogx([x;x],[m.j(ip,:);m.j(ip,:)] + [-sd.j(ip,:);sd.j(ip,:)],'-','color',colors{ip})
end
ylabel('Jeffrey');xlabel('Number of lesioned fibers')
set(gca,'box','off','tickdir','out','xtick',x)
saveFig(fh,fullfile(saveDir,figName),'eps')

figName = sprintf('s_fe_mt_ips_deletion_tests_KL_across_subjects');
fh =figure('name',figName,'color','w');
for ip = whichquartiles
    semilogx(x,m.k(ip,:),'-','color',colors{ip});
    hold on
    semilogx([x;x],[m.k(ip,:);m.k(ip,:)] + [-sd.k(ip,:);sd.k(ip,:)],'-','color',colors{ip})
end
ylabel('K-L');xlabel('Number of lesioned fibers')
set(gca,'box','off','tickdir','out','xtick',x)
saveFig(fh,fullfile(saveDir,figName),'eps')

% Plot the culling procedure
figName = sprintf('s_fe_mt_ips_deletion_tests_fiber_number');
fh =figure('name',figName,'color','w');
for isbj = 1:size(cull,1), hold on
    for ih = 1:size(cull,2)
        plot(cull(isbj,ih).numtotal(2:end)./cull(isbj,ih).numtotal(1),'-o','color',[.6 .6 .6]./isbj)
    end
end
set(gca,'box','off','tickdir','out','ylim',[0 1])
saveFig(fh,fullfile(saveDir,figName),'eps')

figName = sprintf('s_fe_mt_ips_deletion_tests_rmse');
fh =figure('name',figName,'color','w');
for isbj = 1:size(cull,1), hold on
    for ih = 1:size(cull,2)
        plot(cull(isbj,ih).rmse(1:end),'-o','color',[.6 .6 .6]./isbj)
    end
end
ylabel('RMSE');xlabel('Culling iteration number')
set(gca,'box','off','tickdir','out')
saveFig(fh,fullfile(saveDir,figName),'eps')
       
figName = sprintf('s_fe_mt_ips_deletion_tests_r_rmse');
fh =figure('name',figName,'color','w');
for isbj = 1:size(cull,1), hold on
    for ih = 1:size(cull,2)
        plot(cull(isbj,ih).rrmse(1:end),'-o','color',[.6 .6 .6]./isbj)
    end
end
ylabel('R_{RMSE}');xlabel('Culling iteration number')
set(gca,'box','off','tickdir','out','ylim',[.75 1])
saveFig(fh,fullfile(saveDir,figName),'eps')


end % Main function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveFig(h,figName,type)

% MAke sure the folder to save the figure exists
[p,~,~] = fileparts(figName);
[~,message] = mkdir(p);
if ~isempty(message), disp(sprintf('s_fe_mt_ips_deletion_tests_s_fe_mt_ips_deletion_tests_%s.',message));end

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