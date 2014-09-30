function s_pestilli_etal_fig8_major_fascicles(whichSubject)
%
% load an plots the 20 major fascicles generated with AFQ
%
% Copyright Franco Pestilli 2014 Stanford University


% Get the base directory for the data
datapath = '/marcovaldo/frk/2t1/predator/';
subjects = {...
            'KK_96dirs_b2000_1p5iso', ...           
            'MP_96dirs_b2000_1p5iso', ...
            'JW_96dirs_b2000_1p5iso', ...
            'HT_96dirs_b2000_1p5iso', ...
            'KW_96dirs_b2000_1p5iso', ...
            'FP_96dirs_b2000_1p5iso', ...
            };
color = getColors;

for isbj = whichSubject
    fasciclesPath  = fullfile(datapath,subjects{isbj},'afq');
    fasciclesFiles = dir(fullfile(fasciclesPath,'*.mat'));
    t1File         = fullfile(datapath,subjects{isbj},'anatomy','t1.nii.gz');
    savedir        = fullfile(datapath,subjects{isbj},'figures','major_fascicles');
    mkdir(savedir);
    anatomy = niftiRead(t1File);
    
    for iFas = 1:length(fasciclesFiles)
        % Load the fascicles
        load(fullfile(fasciclesPath,fasciclesFiles(iFas).name))
        
        fig_name = [subjects{isbj},fasciclesFiles(iFas).name(1:end-4)];
        viewCoords = [-90,0];
        slice = [1 0 0];
        [fig_h, ~, ~] = plotFascicles(fascicles, color, slice, anatomy, viewCoords, fig_name);
        feSavefig(fig_h,'verbose','yes','figName',[fig_name, 'leftSAG'],'figDir',savedir,'figType','jpg');
        close all; drawnow
    
        viewCoords = [0,90];
        slice = [0 0 6];
        [fig_h, ~, ~] = plotFascicles(fascicles, color, slice, anatomy, viewCoords, fig_name);
        feSavefig(fig_h,'verbose','yes','figName',[fig_name, 'AX'],'figDir',savedir,'figType','jpg');
        close all; drawnow
           
        viewCoords = [90,0];
        slice = [-1 0 0];
        [fig_h, ~, ~] = plotFascicles(fascicles, color, slice, anatomy, viewCoords, fig_name);
        feSavefig(fig_h,'verbose','yes','figName',[fig_name, 'rightSAG'],'figDir',savedir,'figType','jpg');
        close all; drawnow
        
        clear fascicles 
    end
end
            
end % Main function

function [fig_h, light_h, brain_h] = plotFascicles(fascicles, color, slice, anatomy, viewCoords, fig_name)
fig_h = figure('name',fig_name,'color','k');
brain_h = mbaDisplayBrainSlice(anatomy, slice);
hold on 
set(gca,'visible','off','ylim',[-108 69],'xlim',[-75 75],'zlim',[-45 78])
for iFas  = 1:length(fascicles)
    fibers_idx = randsample(1:length(fascicles(iFas).fibers),ceil(length(fascicles(iFas).fibers)*1));
    [~, light_h] = mbaDisplayConnectome(fascicles(iFas).fibers(fibers_idx),fig_h,color{iFas},'single');
    delete(light_h)
end
view(viewCoords(1),viewCoords(2))
light_h = camlight('right');
lighting phong;
%set(fig_h,'Units','normalized', 'Position',[0.5 .2 .4 .8]);
drawnow

end

function colors = getColors

% % Prepare the colors for plotting, Left HM warm, Right HM cool
% numColRes = 12;
% allColors = 1:1:numColRes;
% colormaps = {'spri ng','summer','autumn','winter','bone'};
% for iMap = 1:length(colormaps)
% %figure(iMap)
% for iFas  = 1:length(allColors)
% color{iMap, iFas} = getSmoothColor(allColors(iFas),numColRes,colormaps{iMap});
% %plot(iFas,1,'o','markerfacecolor',color{iMap, iFas},'markersize',16); hold on
% %text(iFas-.1,1,sprintf('%i',iFas),'color','w')
% end
% end
% colors = {color{1,10}, color{1,10}, color{1,5}, color{1,5}, ...
%     color{2,5}, color{2,5},   color{3,1}, color{3,1}, ...
%     color{5,8}, color{5,10},  color{3,6}, color{3,6}, ...
%     color{4,6},  color{4,6},  color{4,9},  color{4,9}, ...
%     color{2,3},  color{2,3},  color{4,3},  color{4,3}};
%       
colors = {[233,150,122]./255, [233,150,122]./255, ... % Salmon
          [255,215,0]./255,   [255,215,0]./255, ... % Gold
          [64,224,208]./255,  [64,224,208]./255, ...% Turquise
          [255,99,71]./255,   [255,99,71]./255,  ...% Tomato
          [220 220 220]./255, [220 220 220]./255, ...% Gainsboro
          [220,20,60]./255,   [220,20,60]./255,   ...
          [221,160,221]./255, [221,160,221]./255, ...
          [199,21,133]./255,  [199,21,133]./255, ...
          [230,230,250]./255, [230,230,250]./255, ...
          [100,149,237]./255, [100,149,237]./255};
                
%figure
%for iFas  = 1:length(colors)
%plot(iFas,1,'o','markerfacecolor',colors{iFas},'markersize',16); hold on
%text(iFas-.1,1,sprintf('%i',iFas),'color','w')
%end
end

function color = getSmoothColor(colorNum,totalColors,colorMap,skipRange)

% default return color
color = [0.5 0.5 0.5];

% check arguments
if ~any(nargin == [1 2 3 4])
  help getSmoothColor
  return
end

% default arguments
if ieNotDefined('totalColors'), totalColors = 256;end
if ieNotDefined('colorMap'), colorMap = 'gray';end
if ~any(strcmp(colorMap,{'hsv','gray','pink','cool','bone','copper','flag'}))
  if ~exist(colorMap,'file')
    disp(sprintf('(getSmoothColor) Unknown colormap function %s',colorMap));
    return
  end
end
if ieNotDefined('skipRange')
  if strcmp(colorMap,'gray')
    skipRange = 0.8;
  else
    skipRange = 1;
  end
end

% get colors to choose from
if skipRange > 0
  colors = eval(sprintf('%s(ceil(totalColors*((1-skipRange)+1)))',colorMap));
else
  colors = eval(sprintf('%s(ceil(totalColors*((1+skipRange)+1)))',colorMap));
  colors = colors(end-totalColors+1:end,:);
end  

% select out the right color
if (colorNum >= 1) & (colorNum <= totalColors)
  color = colors(colorNum,:);
else
  % out of bounds. Warn and return gray
  disp(sprintf('(getSmoothColor) Color %i out of bounds [1 %i]',colorNum,totalColors));
end
end

