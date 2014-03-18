function fe = s_fe_track_predator(bval)
%
% This function:
%  - Tracks using mrtrix accessign a few subjects at the time
%  - Generates an optimized connectome from a cadidate connectome using 
%  LIFE method
%
%  fe = s_fe_fit()
% 
% INPUTS:  none
% OUTPUTS: fe structure the optimized life structure
%
% Copyright Franco Pestilli (2013) Vistasoft Stanford University.


% Get the base directory for the data
datapath = '/home/frk/2t1/predator/';
subjects = {'FP_96dirs_b2000_1p5iso','JW_96dirs_b2000_1p5iso', ...
    'KK_96dirs_b2000_1p5iso','KW_96dirs_b2000_1p5iso', ...
    'MP_96dirs_b2000_1p5iso','HT_96dirs_b2000_1p5iso'};

for isbj = 1;length(subjects)
    % Build the file names for the diffusion data, the anatomical MR, the fiber
    % group containing the connectome and the
    dwiPath  = fullfile(datapath,subjects{isbj},'diffusion_data');
    dwiFiles = dir(fullfile(dwiPath,sprintf('*%s*.gz',num2str(bval))));
    dwiFile       = fullfile(dwiPath,dwiFiles(1).name);
    dwiFileRepeat = fullfile(dwiPath,dwiFiles(2).name);
    t1File        = fullfile(datapath,subjects{isbj},'anatomy','t1.nii.gz');
    
    % Directory where to save the fe structures
    savedir       = fullfile(datapath,subjects{isbj},'connectomes');
    
    % Now find all the fiber files that we will analyze
    fibersPath    = fullfile(datapath,subjects{isbj},'fibers');
    fgFiles       = dir(fullfile(fibersPath,sprintf('*%s*.pdb',num2str(bval))));
    
    % We build one modelper fiber group, whole brain fiber group
    for iFib = 1:length(fgFiles)
        % The final connectome and dat astructure will be saved with this name:
        [~,feFileName] = fileparts(fgFiles(iFib).name);
        
        % Buil a full-file of the fibers and the FE structure
        fgFileName = fullfile(fibersPath,fgFiles(iFib).name);
        feTempPath = fullfile(savedir,feFileName);
        
        if ~(exist([feTempPath,'.mat'],'file') == 2)
            % Intialize a local matlab cluster if the parallel toolbox is available.
            feOpenLocalCluster;
            
            % Initialize the Connectome
            fe = feConnectomeInit(dwiFile,fgFileName,feFileName,savedir,dwiFileRepeat,t1File);
            feConnectomeSave(fe);
            
            M    = feGet(fe,'mfiber');
            dSig = feGet(fe,'dsigdemeaned');
            clear fe
            
            % Fit the model and cull. This will take some time...
            fit = feFitModel(M,dSig,'bbnnls');
            clear M dSig
            
            % Reload the Fe structure, add the fit and save it back to disk
            disp('Re-loading the Fe structure to save the fit...')
            load([feTempPath,'.mat'])
            
            fe = feSet(fe,'fit',fit);
            clear fit
            
            % Save it
            feConnectomeSave(fe);
            
            % Close the matlab pool, this will limit the amount of memory build up
            matlabpool close
        else
            disp('Found FE file not re-computing...')
        end
    end
end

% Exit matlab
exit

return

