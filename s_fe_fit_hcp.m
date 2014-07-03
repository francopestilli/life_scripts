function fe = s_fe_fit_hcp(bval,hemisphere,doClipping,dataDir,which_subject)
%
% This function:
%  - Initializes a LIFE structure from a candidate connectome
%  - Generates an optimized connectome from a cadidate connectome using 
%  LIFE method
%
%  fe = s_fe_fit_hcp()
% 
% INPUTS:  none
% OUTPUTS: fe structure the optimized life structure
%
% Copyright Franco Pestilli (2014) Vistasoft Stanford University
if notDefined('bval'); bval=2000;end
if notDefined('hemisphere'); hemisphere='both';end
if notDefined('doClipping');
    doClipping=true;
end
if notDefined('dataDir'); dataDir='2t1';end

% Get the base directory for the data
[~,hostname] = system('hostname');
hostname = deblank(hostname);
switch dataDir
    case {'2t2'}
        switch hostname
            case {'marcovaldo'}
                datapath = '/home/frk/2t2/HCP/';
            otherwise
                datapath = '/marcovaldo/frk/2t2/HCP/';
                
        end
        subjects = {...
            '115320', ...
            '117122', ...
            '118730', ...
            };
        
    case {'2t1'}
        switch hostname
            case {'marcovaldo'}
                datapath = '/home/frk/2t1/HCP/';
            otherwise
                datapath = '/marcovaldo/frk/2t1/HCP/';
                
        end
        subjects = {...
            '111312', ...
            '113619', ...
            '105115', ...
            '110411', ...
            };
    otherwise
        keyboard
end

for isbj = which_subject
    % Build the file names for the diffusion data, the anatomical MR, the fiber
    % group containing the connectome and the
    dwiPath  = fullfile(datapath,subjects{isbj},'diffusion_data');
    fprintf('\n Loading data from: %s.',dwiPath)
    if ~isempty(bval)
        dwiFiles = dir(fullfile(dwiPath,sprintf('*%s*.gz',num2str(bval))));
    else
        keyboard
    end
    dwiFile       = fullfile(dwiPath,dwiFiles(2).name);
    dwiFileRepeat = fullfile(dwiPath,dwiFiles(2).name);
    t1File        = fullfile(datapath,subjects{isbj},'anatomy','t1.nii.gz');
    
    % Directory where to save the fe structures
    savedir       = fullfile(datapath,subjects{isbj},'connectomes');
    
    % Now find all the fiber files that we will analyze
    fibersPath    = fullfile(datapath,subjects{isbj},'fibers');
    fgFiles       = dir(fullfile(fibersPath,sprintf('*%s*.pdb',num2str(bval))));
    
    % Define the bounding box to clip the fibers
    % (Posterior portion of the left and right hemisphere)
    switch hemisphere
        case {'left'}
            clip{1} = [0,   80]; % Right - Left deleted
            clip{2} = [-8,   90]; % Anterior - Posterior deleted
            clip{3} = [43,   90]; % Superior - Inferior deleted
        case {'right'}
            clip{1} = [-80,  0]; % Right - Left deleted
            clip{2} = [-8,   90]; % Anterior - Posterior deleted
            clip{3} = [43,  90]; % Superior - Inferior deleted   
        case {'both'}
            clip{1} = [];  % NO Right - Left deleted
            clip{2} = [-8,   90]; % Anterior - Posterior deleted
            clip{3} = [43,  90]; % Superior - Inferior deleted

        otherwise
            disp('computing th whole brain conectome... no clipping applied.')
    end
    
    
    % We build one modelper fiber group, whole brain fiber group
    for iFib = [1]
        % Intialize a local matlab cluster if the parallel toolbox is available.
        feOpenLocalCluster;
        
        % The final connectome and dat astructure will be saved with this name:
        [~,feFileName] = fileparts(fgFiles(iFib).name);
        feFileName = [feFileName,'_recomputed'];
        
        % Build a full-file of the fibers and the FE structure
        fgFileName = fullfile(fibersPath,fgFiles(iFib).name);
        if ~doClipping
            hemisphere = '';
        end
        feFileName = deblank([feFileName,hemisphere]);
        feTempPath = fullfile(savedir,feFileName);
        
        if ~(exist([feTempPath,'.mat'],'file') == 2)
            % Load the fibers and clip them to the correct hemisphere
            fg = fgRead(fgFileName);
            if doClipping
                disp('Clipping the fiber group...')
                fg = dtiClipFiberGroup(fg,clip{1},clip{2},clip{3});
            end
            
            % Initialize the Connectome
            fe = feConnectomeInit(dwiFile,fg,feFileName,savedir,dwiFileRepeat,t1File);
            clear feFileName
            M    = feGet(fe,'mfiber');
            dSig = feGet(fe,'dsigdemeaned');
            
            % Fit the model and cull. This will take some time...
            fit = feFitModel(M,dSig,'bbnnls');
                        
            fe = feSet(fe,'fit',fit);
            clear fit
            
            % Save it
            feConnectomeSave(fe);
            
            % Close the matlab pool, this will limit the amount of memory build up
            matlabpool close force local
        else
            disp('Found FE file not re-computing...')
        end
    end
end

end

