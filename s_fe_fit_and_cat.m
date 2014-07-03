function fe = s_fe_fit_and_cat(subjects2run,bval)
%
% This function:
%  - Initializes a LIFE structure from a candidate connectome
%  - Generates an optimized connectome from a cadidate connectome using
%  LIFE method
%  - Eliminates the 0-weight fibers
%  - Concatenates several connectomes together
%
%  fe = s_fe_fit_and_cat()
%
% INPUTS:  none
% OUTPUTS: fe structure the optimized life structure
%
% Copyright Franco Pestilli (2013) Vistasoft Stanford University

% Get the base directory for the data
if strcmpi(deblank(evalc('!hostname')),'marcovaldo') || ...
        strcmpi(deblank(evalc('!hostname')),'black')
    datapath = '/marcovaldo/frk/2t1/predator';
elseif   strcmpi(deblank(evalc('!hostname')),'celadon')
    datapath = '/marcovaldo/frk/2t1/predator';
else    datapath = '/hsgs/projects/wandell/frk/life';
end
fprintf('\n datapath: %s \n',datapath)
subjects = {...
    'JW_96dirs_b2000_1p5iso', ...
    'KK_96dirs_b2000_1p5iso', ...
    'KW_96dirs_b2000_1p5iso', ...
    'MP_96dirs_b2000_1p5iso', ...
    'HT_96dirs_b2000_1p5iso', ...
    'FP_96dirs_b2000_1p5iso', ...
    };

if notDefined('subjects2run');
    subjects2run = 1:length(subjects);
end

for isbj = subjects2run
    fprintf('\n\n\n subject: %s \n\n\n\n',subjects{isbj})
    
    % Build the file names for the diffusion data, the anatomical MR, the fiber
    % group containing the connectome and the
    dwiPath  = fullfile(datapath,subjects{isbj},'diffusion_data');
    if ~isempty(bval)
        dwiFiles = dir(fullfile(dwiPath,sprintf('*%s*.gz',num2str(bval))));
    else
        dwiFiles = dir(fullfile(dwiPath,sprintf('run*.gz')));
    end
    dwiFile       = fullfile(dwiPath,dwiFiles(1).name);
    dwiFileRepeat = fullfile(dwiPath,dwiFiles(2).name);
    t1File        = fullfile(datapath,subjects{isbj},'anatomy','t1.nii.gz');
    
    % Directory where to save the fe structures
    savedir       = fullfile(datapath,subjects{isbj},'connectomes');
    
    % Now find all the fiber files that we will analyze
    fibersPath    = fullfile(datapath,subjects{isbj},'fibers');
    if ~isempty(bval)
        fgFiles       = dir(fullfile(fibersPath,sprintf('*%s*.pdb',num2str(bval))));
    else
        fgFiles       = dir(fullfile(fibersPath,sprintf('*.pdb')));
    end
    keyboard
    % We build one modelper fiber group, whole brain fiber group
    for iFib = [5 1 3 4 2] %length(fgFiles):-1:1
        
        % The final connectome and dat astructure will be saved with this name:
        [~,feFileName] = fileparts(fgFiles(iFib).name);
        
        diary_file = fullfile(datapath,sprintf('diary-%s-%s-%s-%s.txt', ...
            mfilename,subjects{isbj},feFileName,date));
        
        if ~(exist([diary_file],'file') == 2)
            diary(diary_file)
            
            % Build a full-file of the fibers and the FE structure
            fgFileName = fullfile(fibersPath,fgFiles(iFib).name);
            feTempPath = fullfile(savedir,feFileName);
            
            if ~(exist([feTempPath,'.mat'],'file') == 2)
                % Initialize the Connectome
                fe = feConnectomeInit(dwiFile,fgFileName,feFileName,savedir,dwiFileRepeat,t1File);
                feConnectomeSave(fe);
                
                M    = feGet(fe,'mfiber');
                dSig = feGet(fe,'dsigdemeaned');
                clear fe
                
                % Fit the model and cull. This will take some time...
                fit = feFitModel(M,dSig,'bbnnls');
                clear M dSig
                
                fe = feSet(fe,'fit',fit);
                clear fit
                
                % Save it
                feConnectomeSave(fe);
                
            else
                fprintf('\n[%s] Found FE file not re-computing...\n\n%s\n\n',mfilename,  feTempPath)
            end
        else
            fprintf('\n\nFOUND DIARY FILE: \n%s \n\n\n',diary_file)
        end
        fprintf('\n\n\n\n\nDONE processing: %s ...\n\n\n\n',diary_file)
        diary off
    end
    
    fprintf('\n\n\n DONE subject: %s \n\n\n\n',subjects{isbj})
    
end
        
% Exit matlab
exit