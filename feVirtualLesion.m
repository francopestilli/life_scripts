function [streghtOfEvidence, fh, feLesion,  feNoLesion, connectivity, newFascicleIndices, indicesFibersKept, commonCoords] = ...
          feVirtualLesion(feNoLesion,                  fascicleIndices, refitConnectome, displayFibers)
%
% [feLesion, feNoLesion, connectivity, newFascicleIndices, indicesFibersKept, commonCoords] = ...
%  feTestFascicle(feNoLesion,fascicleIndices,[refitConnectome],[displayFibers])
%
% Perform a virtual lesion.
% - Remove a set of fascicles from a whole-brain connectome
% - Find the fascicle's path-neighborhood (the set of fascicles that share the same voxels)
% - Reduce the connectome to only the voxels of the set of fascicels that
%   we are lesioning. These are the voxels that contain the signal and the
%   error useful for the lesion test.
% - Compute a series of statistics on the lesioned connectome.
%
% Copyright Franco Pestilli Stanford University 2014

if notDefined('displayFibers'),   displayFibers = 0;end
if notDefined('refitConnectome'), refitConnectome=0;end

% Handling parallel processing
poolwasopen=1; % if a matlabpool was open already we do not open nor close one
if (matlabpool('size') == 0), matlabpool open; poolwasopen=0; end

tic,fprintf('\n[%s] Performing the lesion (removing fascicles from the connectome)... ',mfilename)
if ~islogical( fascicleIndices )
    nfibers = feGet(feNoLesion,'nfibers');
    fascicles2keep = false(nfibers,1);
    fascicles2keep(fascicleIndices) = true;
else
    fascicles2keep = fascicleIndices;
end
feLesion = feConnectomeReduceFibers(feNoLesion, fascicles2keep );toc

% Extract the fascicle out of the fiber group.
fas = fgExtract(feGet(feNoLesion,'fibers img'), fascicleIndices, 'keep' );

% Get the cordinates of the fascicle we just deleted. These coordinates are
% contained in the full connectome. We want to fin the indices in the M
% matrix to the the voxels of the fascicle and keep only those voxels.
tic, fprintf('\n[%s] Finding the voxels of fascicle... ',mfilename)
fasCoords    = fefgGet(fas,'unique image coords');
allCoords    = feGet(feNoLesion,   'roi coords');
commonCoords = ismember(allCoords, fasCoords,  'rows'); toc

% Now: commonCoords contains the indices of the voxels to keep from feLesion,
% these voxels are part of the connectome feNoLesion. So now we want to delete all
% the rest of the voxels from feLesion and keep only the voxels where
% the fascicle in feFP goes through.
% 
% At the same time we keep all the fibers in the connectome that still pass
% throught the left voxels in the ROI. We will use this subset of voxels
% and fibers to test the quality of fit of the connectome at the location
% where the feFN fascicle went through.
tic,fprintf('\n[%s] Creating a lesioned connectome... ',mfilename)
feLesion = feConnectomeReduceVoxels(feLesion,find(commonCoords));
toc

tic,fprintf('\n[%s] Creating an unleasioned path-neighborhood connectome... ',mfilename)
[feNoLesion, indicesFibersKept]    = feConnectomeReduceVoxels(feNoLesion, find(commonCoords));
toc

tic,fprintf('\n[%s] Finding the pathneighborhood voxels... ',mfilename)
% Here we return the indices of the fascicle in the newly resized
% connectome.
newFascicleIndices = ismember(find(indicesFibersKept),...
                              find(fascicleIndices),'rows');
toc

% Fit the connectomes unlesioned the leftover fibers.
if refitConnectome
    % Refit and install the fit.
    tic, fprintf('\n[%s] Fitting connectome lesioned the fascicle... ',mfilename)
    feLesion = feSet(feLesion,'fit', ...
        feFitModel(feGet(feLesion,'Mfiber'),feGet(feLesion,'dsigdemeaned'), ...
        'bbnnls'));toc
    
    tic, fprintf('\n[%s] Fitting connectome unlesioned the fascicle... ',mfilename)
    feNoLesion    = feSet(feNoLesion,   'fit', ...
        feFitModel(feGet(feNoLesion,'Mfiber'),   feGet(feNoLesion,'dsigdemeaned'), ...
        'bbnnls'));toc
end
keyboard
% Now compute the connectivity measure of the voxel.
% This is the sum of the weights of the fascile divided by the sum of the
% weights of all the fibers in the same volume of white-matter.
%
% Find the weights for the path-neighborhood fibers.
connectivity.wall = feGet(feNoLesion,'fiber weights');
connectivity.wfas = connectivity.wall(newFascicleIndices);
connectivity.wnfas = connectivity.wall(~newFascicleIndices);

% Compute some measures of strength of the connection represented by the
% fascicle, by comparing the fascicle weights unlesioned the weights of all the
% rest of the fascicles going through the same voxels.
connectivity.strength(1) =    sum(connectivity.wfas) /  sum(connectivity.wnfas);
connectivity.strength(2) =   mean(log10(connectivity.wfas)) / mean(log10(connectivity.wnfas));
connectivity.strength(3) = median(log10(connectivity.wfas))/median(log10(connectivity.wnfas));

% Compute S the strength of evidence in favor of the lesioned fascicle
% Make a plot of the R-squared
unlesioned.rmse      = mean(feGetRep(feNoLesion,   'vox  rmse'));
unlesioned.rrmse     = median(feGetRep(feNoLesion,   'vox  rmse ratio'));
unlesioned_rmseall   = (feGetRep(feNoLesion,   'vox  rmse'));

lesioned.rmse    = median(feGetRep(feLesion,'vox  rmse'));
lesioned.rrmse   = median(feGetRep(feLesion,'vox  rmse ratio'));
lesioned_rmseall = (feGetRep(feLesion,'vox  rmse'));

%% The following is the code for the bootstrap test on the MEAN rmse
nboots=10000; nmontecarlo = 10;
sizeunlesioned    = length(unlesioned_rmseall);
nullDistributionW = nan(nboots,nmontecarlo);
nullDistributionWO = nan(nboots,nmontecarlo);
min_x = floor(mean([unlesioned_rmseall]) - mean([unlesioned_rmseall])*.05);
max_x = ceil(mean([lesioned_rmseall]) + mean([lesioned_rmseall])*.05);

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

% Plot the null distribution and the empirical difference
figName = sprintf('virtual_lesion_test_mean_rmse_hist_%s_%s',mfilename,feLesion.name);
fh = mrvNewGraphWin(figName);
patch([xhis,xhis],y_e(:),'b','FaceColor',[.67 .86 .96],'EdgeColor',[.67 .86 .96]); % Distribution as the +/- 2SD
hold on
patch([woxhis,woxhis],ywo_e(:),[.97 .66 .76],'FaceColor',[.97 .66 .76],'EdgeColor',[.97 .66 .76]); % Distribution as the +/- 2SD
set(gca,'tickdir','out', ...
        'box','off', ...
        'ylim',[0 0.45], ... 
        'xlim',[min_x,max_x], ...
        'ytick',[0 0.2 0.4], ...
        'xtick',round(linspace(min_x,max_x,4)), ...
        'fontsize',16)
ylabel('Probability','fontsize',16)
xlabel('rmse','fontsize',16')

% (3) Compute the probability that the empirical difference (1) was
%     observed by chance given th data, by looking at the percentile of the
%     empirical difference in the Nul distribution (2).
streghtOfEvidence = mean(diff([mean(nullDistributionW,1); ...
                          mean(nullDistributionWO,1)])./sqrt(sum([std(nullDistributionW,[],1);std(nullDistributionWO,[],1)].^2,1)));
title(sprintf('Strength of connection evidence %2.3f',(streghtOfEvidence)), ...
    'FontSize',16)

% The following is test code to show where the coordinates of the facicle
% that were removed land inside the connectoem. Also I show in gray hte
% connectome lesioned the fascicle and in red the connectome unlesioned the
% fascicle. Where there is only red in the connectoem that is where the
% fascicle was removed but we are attemtping to explain the variance in the
% data.
if displayFibers
    % Show the coordinates to see if we are in the right spot
    mrvNewGraphWin('Coordinate check');
    plot3(allCoords(commonCoords,1),allCoords(commonCoords,2), ...
          allCoords(commonCoords,3),'ko','MarkerFaceColor','k','MarkerSize',8);
    hold on;
    plot3(fasCoords(:,1),fasCoords(:,2),fasCoords(:,3),'ro', ...
          'MarkerFaceColor','r','MarkerSize',3);
    axis equal
    %view(3,79)
    view(-23,-23);
    
    % Now display the fascicle removed and the connectome lesioned the fascicle
    % in one figure unlesioned different colors.
    %
    feConnectomeDisplay(feSplitLoopFibers(feGet(feLesion,'fibers img')),figure);
    hold on
    %feConnectomeDisplay(feSplitLoopFibers(feGet(feNoLesion,'fibers img')),gcf, [.95 .1 .1])   
    feConnectomeDisplay(feSplitLoopFibers(fas),gcf,[.95 .1 .1]);
    view(-23,-23);
    h= camlight;

    feConnectomeDisplay(feSplitLoopFibers(feGet(feLesion, ...
                        'fibers img')),figure);
    hold on
    plot3(allCoords(commonCoords,1),allCoords(commonCoords,2), ...
        allCoords(commonCoords,3),'go','MarkerFaceColor','g','MarkerSize',12);
    view(-23,-23);
    h= camlight;       
end

if ~poolwasopen, matlabpool close; end

end