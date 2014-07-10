function fe = s_fe_fit_predator_cull(subjects_to_run)
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
datapath = '/marcovaldo/frk/2t1/predator/';
subjects = {...
    'FP_96dirs_b2000_1p5iso', ...
    'KW_96dirs_b2000_1p5iso', ...
    'KK_96dirs_b2000_1p5iso', ...
    'MP_96dirs_b2000_1p5iso', ...
    'HT_96dirs_b2000_1p5iso', ...
    'JW_96dirs_b2000_1p5iso', ...
    };


for isbj = subjects_to_run
    % Directory where to save the fe structures
    savedir       = fullfile(datapath,subjects{isbj},'connectomes');   
    feFiles       = dir(fullfile(savedir,sprintf('*_recomputed.mat')));

    % We build one modelper fiber group, whole brain fiber group
    for iFes = length(feFiles):-1:1
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

