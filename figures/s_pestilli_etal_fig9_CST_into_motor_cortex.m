function s_pestilli_etal_fig9_CST_into_motor_cortex()
%
% Show that the CST in the optimized conenctome contains projections into
% lateral portion of the human cortex.
% 
% This is a valuable aspect of the fascicle tha tis difficult to obtain
% with deterministic tractogphy but possible with deterministic
% tractogrpahy.
%
% This example shows that LiFE does not remove this valuable feature.
%
% Copyright Franco Pestilli 2014 Stanford University

% Get the base directory for the data
datapath = '/marcovaldo/frk/2t1/predator/';
subjects = {...
    'KW_96dirs_b2000_1p5iso', ...
    'MP_96dirs_b2000_1p5iso', ...
    'HT_96dirs_b2000_1p5iso', ...
    'FP_96dirs_b2000_1p5iso', ...'JW_96dirs_b2000_1p5iso', ...
    'KK_96dirs_b2000_1p5iso', ...
    };
tractography = {'lmax10','prob'};
for isbj = 1:length(subjects)
    fasciclesPath  = fullfile(datapath,subjects{isbj},'afq','major_fascicles');
    savedir        = fullfile(datapath,subjects{isbj},'afq','CST');
    optimized_connectome  = matchfiles(fullfile(datapath,subjects{isbj},'fibers', ...
                            sprintf('*%s*%s*recomputed-optimized.mat',tractography{1},tractography{2})));
    t1File   = fullfile(datapath,subjects{isbj},'anatomy','t1.nii.gz');
    dtFile   = fullfile(datapath,subjects{isbj},'dtiInit','dt6.mat');
    
    fprintf('[%s] Subject: %s\n',mfilename,subjects{isbj})
    fg_names = matchfiles(fullfile(fasciclesPath,'*Cortico*'));
    for iFas = 1:length(fg_names)
        % Load the Left and Right Cortico Spinal Tract
        fg(iFas) = fgRead(fg_names{iFas});
    end
    
    % Merge the two fiber groups
    fg = fgMerge(fg(1),fg(2));
    
    % Make an ROI out of them.
    %roi = dtiNewRoi('L_R_CST', [], unique(floor(horzcat(fg.fibers{:}))','rows'));
    roi = dtiCreateRoiFromFibers(fg);
    clear fg
    % Show the ROI:
    % plot3(roi.coords(:,1),roi.coords(:,2),roi.coords(:,3),'ro');view(-90,0)
    
    % Clip the ROI to the AC-PC plane
    roi = dtiRoiClip(roi, [], [], [10 90]); 
    % Show it after the trasformation: 
    % hold on; plot3(roi.coords(:,1),roi.coords(:,2),roi.coords(:,3),'k*');

    % Smooth and expand the ROI
    switch subjects{isbj}
        case {'MP_96dirs_b2000_1p5iso','KW_96dirs_b2000_1p5iso'}
            roi    = dtiRoiClean(roi,6);
        otherwise
            roi    = dtiRoiClean(roi,6,['dilate']);
    end
    
    % Show it after the trasformation: 
    % hold on;
    % plot3(roi.coords(:,1),roi.coords(:,2),roi.coords(:,3),'c*');view(0,0)

    % Remove fibers going through the Corpus callosum. These are not the
    % topic of this example
    cc        = dtiNewRoi('CC');
    dt = dtiLoadDt6(dtFile);
    cc.coords = dtiFindCallosum(dt.dt6,dt.b0,dt.xformToAcpc);
    cc        = dtiRoiClean(cc,0,['dilate']);  
    cc        = dtiRoiClean(cc,0,['dilate']); 
    cc        = dtiRoiClean(cc,0,['dilate']);
    clear dt
    % hold on;
    % plot3(cc.coords(:,1),cc.coords(:,2),cc.coords(:,3),'r^');view(90,0)

    % Load an Optimized connectome
    fg = fgRead(optimized_connectome{1});
    
    % Intersect the connectome with the CST ROI
    tic, fprintf('[%s] Removing fibers not going through the CST ROI... ',mfilename)
    fg = feSegmentFascicleFromConnectome(fg, {roi,cc}, {'and','not'}, 'L_R_CST');
    
    % Load the T1 file for display
    t1     = niftiRead(t1File);
    
    % Show te new fiber group
    figName = sprintf('Corona_radiata');
    figureHandle = figure('name',figName);
    h  = mbaDisplayBrainSlice(t1, [0 3 0]);
    hold on
    tractColor = [.2  .4 .95]; %subset = (randsample(1:length(fg.fibers),100));
    [figureHandle, lightHandle] = mbaDisplayConnectome(mbaFiberSplitLoops(fg.fibers),figureHandle,tractColor,'uniform');
    delete(lightHandle);
    view(0,0);
    axis([-65 65 -50 90 -40 75]);
    lightHandle = camlight('right');
    drawnow
    saveFig(figureHandle,fullfile(savedir,figName),'jpeg')
    close all
    
    fprintf('[%s] Saving corona radiata\n',mfilename)
    fgWrite(fg,fullfile(savedir,'Corona_radiata'),'mat');
    clear fg t1 
 end

end % Main function

%-------------------------------%
function saveFig(h,figName,eps)
if ~exist( fileparts(figName), 'dir'), mkdir(fileparts(figName));end
fprintf('[%s] saving figure... \n%s\n',mfilename,figName);

switch eps
    case {0,'jpeg','jpg'}
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