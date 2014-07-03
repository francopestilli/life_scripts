function s_fe_mt_ips_plot_S()

if notDefined('saveDir'), saveDir = fullfile('/marcovaldo/frk/Dropbox','pestilli_etal_revision',mfilename);end

dataset = {'96','hcp'};
for id = 1:length(dataset)
    switch dataset{id}
        case '96'
            cd /home/frk/Dropbox/pestilli_etal_revision/s_fe_mt_ips_tract
            load strength_of_evidence.mat
            eylims = [0 26];
            eyticks = [0 13 26];
        case 'hcp'
            cd /home/frk/Dropbox/pestilli_etal_revision/s_fe_mt_ips_tract_hcp/
            load strength_of_evidence.mat
            eylims = [0 500];
            eyticks = [0 250 500];
        otherwise
            keyboard
    end

% Extract all the metrics in a easy format for plotitng
for is = 1:size(SE,1)
    for ih = 1:size(SE,2)
        kl(is,ih) = SE(is,ih).kl.mean;
        S(is,ih) = SE(is,ih).s.mean;      
        J(is,ih) = SE(is,ih).j.mean;
        E(is,ih) = SE(is,ih).em.mean;
    end
end

mKL = mean(kl(:));
seKL = std(kl(:))./sqrt(numel(kl));
figName = sprintf('MT_parietal_connection_KL_%s_%s',dataset{id},mfilename);
fh(1) = mrvNewGraphWin(figName);
bar(1,mKL,'k');
hold on
plot([1; 1]',[mKL; mKL] + [seKL; -seKL],'r', 'linewidth',4)
set(gca,'tickdir','out', ...
        'box','off', ...
        'ylim',[0 6], ... 
        'xlim',[0 2], ...
        'ytick',[0 3 6], ...
        'xtick',[1], ...  
        'xticklabel',{''}, ...
        'fontsize',16)
ylabel(sprintf('Kullback-Leibler divergence'),'fontsize',16)
saveFig(fh(1),fullfile(saveDir,figName),'eps')
 
mS  = mean(S(:));
seS = std(S(:)) ./ sqrt(numel(S(:)));

figName = sprintf('MT_parietal_connection_S_%s_%s',dataset{id},mfilename);
fh(1) = mrvNewGraphWin(figName);
bar([1],mS,'k');
hold on
plot([1; 1]',[mS; mS] + [seS; -seS],'r', 'linewidth',4)
set(gca,'tickdir','out', ...
        'box','off', ...
        'ylim',[0 50], ... 
        'xlim',[0 2], ...
        'ytick',[0 25 50], ...
        'xtick',[1 2], ...  
        'xticklabel',{''}, ...
        'fontsize',16)
ylabel('Strength of evidence (S)','fontsize',16)
saveFig(fh(1),fullfile(saveDir,figName),'eps')

mJ  = mean(J(:));
seJ = std(J(:)) ./ sqrt(numel(J(:)));

figName = sprintf('MT_parietal_connection_JEFFREY_%s_%s',dataset{id},mfilename); 
fh(1) = mrvNewGraphWin(figName);
bar([1],mJ,'k');
hold on
plot([1; 1]',[mJ; mJ] + [seJ; -seJ],'r', 'linewidth',4)
set(gca,'tickdir','out', ...
        'box','off', ...
        'ylim',[0 2], ... 
        'xlim',[0 2], ...
        'ytick',[0 1 2], ...
        'xtick',[1 2], ...  
        'xticklabel',{''}, ...
        'fontsize',16)
ylabel('Jeffrey''s divergence (bits)','fontsize',16)
saveFig(fh(1),fullfile(saveDir,figName),'eps')

mE  = mean(E(:));
seE = std(E(:)) ./ sqrt(numel(E(:)));

figName = sprintf('MT_parietal_connection_EMD_%s_%s',dataset{id},mfilename); 
fh(1) = mrvNewGraphWin(figName);
bar([1],mE,'k');
hold on
plot([1; 1]',[mE; mE] + [seE; -seE],'r', 'linewidth',4)
set(gca,'tickdir','out', ...
        'box','off', ...
        'ylim',eylims, ... 
        'xlim',[0 2], ...
        'ytick',eyticks, ...
        'xtick',[1 2], ...  
        'xticklabel',{''}, ...
        'fontsize',16)
ylabel('Earth movers'' distance','fontsize',16)
saveFig(fh(1),fullfile(saveDir,figName),'eps')

end


clear kl mKL seKL mS seS S mS seS S

% %% Human Connectome data
% 
% 
% for is = 1:size(KL,1)
%     for ih = 1:size(KL,2)
%         kl(is,ih) = KL(is,ih).value;
%     end
% end
% 
% kl(isinf(kl)) = nan;
% kl(isnan(kl)) = nanmax(kl(:));
% mKL = mean(kl);
% seKL = std(kl)./sqrt(size(kl,1));
% 
% figName = sprintf('MT_parietal_connection_KL_HCP_%s',mfilename);
% fh(1) = mrvNewGraphWin(figName);
% bar([1 2],mKL,'k');
% hold on
% plot([1 1; 2 2]',[mKL; mKL] + [seKL; -seKL],'r', 'linewidth',4)
% set(gca,'tickdir','out', ...
%         'box','off', ...
%         'ylim',[0 4], ... 
%         'xlim',[0 3], ...
%         'ytick',[0 2 4], ...
%         'xtick',[1 2], ...  
%         'xticklabel',{'Right','Left'}, ...
%         'fontsize',16)
% ylabel(sprintf('Strength of evidence\n(Kullback-Leibler divergence)'),'fontsize',16)
% xlabel('Hemisphere','fontsize',16')
% saveFig(fh(1),fullfile(saveDir,figName),'eps')
% 
% 
% mS  = mean(S,1);
% seS = std(S,[],1) ./ sqrt(size(S,1));
% 
% figName = sprintf('MT_parietal_connection_S_HCP_%s',mfilename);
% fh(1) = mrvNewGraphWin(figName);
% bar([1 2],mS,'k');
% hold on
% plot([1 1; 2 2]',[mS; mS] + [seS; -seS],'r', 'linewidth',4)
% set(gca,'tickdir','out', ...
%         'box','off', ...
%         'ylim',[0 50], ... 
%         'xlim',[0 3], ...
%         'ytick',[0 25 50], ...
%         'xtick',[1 2], ...  
%         'xticklabel',{'Right','Left'}, ...
%         'fontsize',16)
% ylabel('Strength of evidence (S)','fontsize',16)
% xlabel('Hemisphere','fontsize',16')
% saveFig(fh(1),fullfile(saveDir,figName),'eps')
%  
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