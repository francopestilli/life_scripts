function s_fe_run_afq(bval,tractographyType,optimized,doAFQoutliers)
%
% This function:
%  - Load a series of culled connectomes.
%  - Run AFQ on each one of them.
%  - Make a plots of the fascicles.
%
%  s_fe_run_afq(bval,tractographyType)
% 
% Copyright Franco Pestilli (2014) Stanford University

% Get the base directory for the data
datapath = '/marcovaldo/frk/2t1/predator/';
subjects = {...
            'FP_96dirs_b2000_1p5iso', ...
            'HT_96dirs_b2000_1p5iso', ...
            'KK_96dirs_b2000_1p5iso', ...           
            'MP_96dirs_b2000_1p5iso', ...
            'JW_96dirs_b2000_1p5iso', ...
            'KW_96dirs_b2000_1p5iso', ...
            };
        
if notDefined('optimized'), optimized = true; end  
if notDefined('tractographyType'), tractographyType = 'lmax10'; end
if notDefined('bval'), bval = []; end
if notDefined('doAFQoutliers')
    doAFQoutliers = ''; 
else
    doAFQoutliers = 'noAFQclean'; 
end

% These are default parameters for plotting
for isbj = 1:length(subjects)
    % Directory where to load the whome-brain fiber groups
    fiberPaths = fullfile(datapath,subjects{isbj},'fibers');
    dtFile    = fullfile(datapath,subjects{isbj},'dtiInit','dt6.mat');
    
    if ~isempty(bval)
        if optimized
            fibersFiles       = dir(fullfile(fiberPaths,sprintf('*%s*-optimized*.mat',num2str(bval))));
        else
            fibersFiles       = dir(fullfile(fiberPaths,sprintf('*%s*.pdb',num2str(bval))));
        end
    else
        if optimized
            fibersFiles       = dir(fullfile(fiberPaths,sprintf('*%s*-optimized*.mat',tractographyType)));
        else
            fibersFiles       = dir(fullfile(fiberPaths,sprintf('*%s*.pdb',tractographyType)));
            
        end
    end
            
    % We build one modelper fiber group, whole brain fiber group
    for iFibers = 1:length(fibersFiles)
        % The final connectome and dat astructure will be saved with this name:
        [~,fibersFileName,ext] = fileparts(fibersFiles(iFibers).name);
        
        % Buil a full-file of the fibers and the RESULTS structure to load
        fibersFileName2Load = fullfile(fiberPaths,[fibersFileName,ext]);

        % Initialize the Connectome
        fprintf('[%s] Loading whole-brain connectome: \n%s\n',mfilename,fibersFileName2Load)        
        % Load the connectome
        fg = fgRead(fibersFileName2Load);
        
        % Segment the fibers using AFQ
        [fg_classified,~,classification]= AFQ_SegmentFiberGroups(dtFile, fg);
        
        % Split the fiber groups into individual groups
        fascicles = fg2Array(fg_classified);
        
        if strcmpi(doAFQoutliers,'')
        % Clean the fibers, we apply the same trhesholds to all fiber
        % groups this is the default thrshold used by AFQ. This is done by
        % not passing opts
        [fascicles, classification] = feAfqRemoveFascicleOutliers(fascicles,classification);
        end
        
        % Save the segemented fascicles and the indices into the Mfiber
        fibersInfoToSaveDir = fullfile(datapath,subjects{isbj},'afq');
        mkdir(fibersInfoToSaveDir)
        fibersInfoToSaveName =  fullfile(fibersInfoToSaveDir,[fibersFileName,'-AFQ',doAFQoutliers,'.mat']);
        save(fibersInfoToSaveName,'fg_classified','classification','fascicles')
     end
end

end % End main function

function [fascicles, classification] = feAfqRemoveFascicleOutliers(fascicles,classification,opts)
%
% Removes the outliers in a group of fascicles.
%
% [fascicles classification] = feAfqRemoveFascicleOutliers(fascicles,classification,opts)
%
% INPUTS:
%    fascicles      - a set of fascicles created using feAfqSegment.m
%    classification - a classification strucutre created using
%                     feAfqSegment.m
%    opts           - default options for cleaning the outliers in the
%                     fascicle.
%
% OUTPUTS:
%    fascicles      - a set of fascicles as created by feAfqSegment.m
%                     but withut Outliers.
%    classification - a classification strucutre as created using
%                     feAfqSegment.m, but without outliers
%
% NOTEs: it requres AFQ to be on path.
%
% Copyright Franco Pestilli Stanford University 2014

% Parameters for the outliers removal process.
if notDefined('opts')
  opts.stdCutOff   = 3;   % Standard deviation fo the 3D gaussian distribution used to
  % represent the fascicle when removing outliers. 3.5 z-scores
  opts.maxLen      = 25;  % Max lenght of fibers to be accepted in cm
  opts.maxNumNodes = 100; % This is used only during the computations does not actually change the nodes
end

% Remove fibers outliers, this creates tighter fascicles.
keep = cell(length(fascicles),1);
fprintf('[%s] Removing outliers from fascicles...\n',mfilename)
for in = 1:length(fascicles)
  fprintf('[%s] %i/%i %s.\n',mfilename,in,length(fascicles),fascicles(in).name)
  [fascicles(in), keep{in}] = AFQ_removeFiberOutliers(fascicles(in),opts.stdCutOff,opts.maxLen,opts.maxNumNodes);
  
  % Find the indices to the current fascicle inside the vector of indices
  % for the whol-brain connectome (fg)
   thisFasIndices = find((classification.index == in));

  % Now remove the fibers that we do not want to keep out of the fascicles.
  classification.index(thisFasIndices(~keep{in})) = 0;
end

end