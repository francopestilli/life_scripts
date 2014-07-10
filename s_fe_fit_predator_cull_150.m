function fe = s_fe_fit_predator_cull_150()
%
% This function:
%  - Loads a LIFE structure  and culls it
%
%  fe = s_fe_fit_predator_cull(subjects_to_run)
% 
% INPUTS:  none
% OUTPUTS: fe structure the optimized life structure
%
% Copyright Franco Pestilli (2014) Vistasoft Stanford University

if notDefined('recompute'), recompute = true; end % is true will recumpute the file and overwrite it

% Get the base directory for the data
% Get the base directory for the data
datapath = '/marcovaldo/frk/2t2/predator/';
subjects = {'FP_150dirs_b1000_2000_4000_2iso'};

for isbj = 1
    % Directory where to save the fe structures
    savedir       = fullfile(datapath,subjects{isbj},'connectomes');   
    feFiles       = dir(fullfile(savedir,sprintf('*1000*lmax*.mat')));

    % We build one modelper fiber group, whole brain fiber group
    for iFes = 1:length(feFiles)
        % The final connectome and dat astructure will be saved with this name:
        [~,feFileName] = fileparts(feFiles(iFes).name);
        feTempPath     = fullfile(savedir,feFileName);
        
        if (exist([feTempPath,'.mat'],'file') == 2) 
            if ~(exist([feTempPath,'_culled.mat'],'file') == 2) || recompute == true
            % Intialize a local matlab cluster if the parallel toolbox is available.
            feOpenLocalCluster;
            
            fprintf('[%s] Loading the FE structure...\n %s \n', mfilename,[feTempPath,'.mat'])
            load([feTempPath,'.mat'])
            
            % Initialize the Connectome
            [fe,cull] = feConnectomeCull(fe);
            
            fprintf('[%s] Saving the FE structure and the culling info...\n %s \n', mfilename,[feTempPath,'_culled.mat'])
            save([feTempPath,'_culled.mat'],'fe','cull','-v7.3');
            clear fe
            
            % Close the matlab pool, this will limit the amount of memory build up
            matlabpool close force local
            end
        else
            fprintf('[%s] FE file not found, skipping...\n %s \n', mfilename,[feTempPath,'.mat'])
        end
    end
end

% Exit matlab
exit

return

