function s_fe_plot_directions
%
% Load a file computed by s_ms_directions_96/150/hcp
% Plot the results
%
% Copyright Franco Pestilli Stanford University 2014

saveDir = fullfile('/marcovaldo/frk/Dropbox/pestilli_etal_revision/s_ms_directions_96/average_96_1p5mm/');
load(fullfile(saveDir,'mean_histograms_recomputed.mat'))

fh = figure('name','directions_STN96','color','w');
my = mean(m.nfibers.y,2);
sy = std(m.nfibers.y,[],2)./sqrt(size(m.nfibers.y,2));
y = [my - sy, my + sy]';
x = [m.nfibers.x;m.nfibers.x];
plot(m.nfibers.x,my,'ro-','color',[.45 .45 .95],'markerfacecolor',[.45 .45 .95])
hold on
plot(x,y,'r-','color',[.45 .45 .95])
set(gca,'tickdir','out','box','off','xtick',[8 16 32 64 96],'fontsize',14,'ylim',[5.5*10^4 10*10^4])
ylabel('Number of fascicles','fontsize',14)
xlabel('Number of directions','fontsize',14)
saveFig(fh,fullfile(saveDir,'directions_STN96'),1)

saveDir = fullfile('/marcovaldo/frk/Dropbox/pestilli_etal_revision/s_ms_directions_150/average_150_2mm/');
load(fullfile(saveDir,'mean_histograms.mat'))

fh = figure('name','directions_STN150','color','w');
my = mean(m.nfibers.y,2);
sy = std(m.nfibers.y,[],2)./sqrt(size(m.nfibers.y,2));
y = [my - sy, my + sy]';
x = [m.nfibers.x;m.nfibers.x];
plot(m.nfibers.x,my,'ro-','color',[.45 .45 .95],'markerfacecolor',[.45 .45 .95])
hold on
plot(x,y,'r-','color',[.45 .45 .95])
set(gca,'tickdir','out','box','off','xtick',[8 16 32 64 96],'fontsize',14,'ylim',[5.5*10^4 10*10^4])
ylabel('Number of fascicles','fontsize',14)
xlabel('Number of directions','fontsize',14)
saveFig(fh,fullfile(saveDir,'directions_STN150'),1)
keyboard
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