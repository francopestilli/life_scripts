function s_fe_save_culled_fg_hcp(tractographyType,bval)
%
% This function:
%  - Load a series of precomputed connectomes (fe strucutres)
%  - Extracts the optimized fiber group.
%  - Saves it to the corresponding folder containing the candidate connectome.
%
% The output of this function is used by: s_fe_weights_hcp.m
%
%  fe = s_fe_save_culled_fg()
% 
% Copyright Franco Pestilli (2014) Vistasoft Stanford University

% Get the base directory for the data
% Get the base directory for the data
datapath = '/marcovaldo/frk/2t2/HCP/';
subjects = {'115320','117122','118730'};
if notDefined('tractographyType'), tractographyType = 'lmax10'; end
if notDefined('bval'), bval = []; end

for isbj = 1:length(subjects)
    % Directory where to save the fibers and the results
    fibersSaveDir       = fullfile(datapath,subjects{isbj},'fibers');
    resultsSaveDir       = fullfile(datapath,subjects{isbj},'results');
  
    % Now find all the fiber files that we will analyze
    fePath    = fullfile(datapath,subjects{isbj},'connectomes');
    
    if ~isempty(bval)
       feFiles       = dir(fullfile(fePath,sprintf('*%s*.mat',num2str(tractographyType))));
    else
       feFiles       = dir(fullfile(fePath,sprintf('*.mat')));
    end
            
    % We build one modelper fiber group, whole brain fiber group
    for iFe = 1:length(feFiles)
        % The final connectome and dat astructure will be saved with this name:
        [~,feFileName] = fileparts(feFiles(iFe).name);
        
        % Buil a full-file of the fibers and the FE structure
        feFileName2Load = fullfile(fePath,feFiles(iFe).name);
        fgGoodFileName = fullfile(fibersSaveDir,[feFileName,'-optimized.mat']); 
        fgBadFileName = fullfile(fibersSaveDir,[feFileName,'-rejected.mat']);
        resultsFileName2Save = fullfile(resultsSaveDir,[feFileName,'-fiberStatsResults.mat']);

        % Initialize the Connectome
        fprintf('[%s] Loading a FE: \n%s\n',mfilename,feFileName2Load)
        load(feFileName2Load);

        % Get the weights
        fprintf('[%s] Extracting the weights\n',mfilename)
        %xformimg2acpc = feGet(fe,'xformimg2acpc');
        %mapsize       = feGet(fe,'mapsize');
        w    = feGet(fe,'fiber weights');
        goodFibers = w > 0;
        badFibers  = w == 0;
        results.weights = w;
        
        fprintf('[%s] Extracting the FG...\n',mfilename)
        fg = feGet(fe,'fibers acpc');
        clear fe
    
        fprintf('[%s] Extracting fiber density and length of the candidate FG\n',mfilename)
        %results.candidate.density = dtiComputeFiberDensityNoGUI(fg,xformimg2acpc,mapsize);
        results.candidate.length  = cellfun(@length,fg.fibers); 
        results.candidate.n = length(w);
          
        fprintf('[%s] Extracting the optimized FG\n',mfilename)
        fgB = fgExtract(fg,badFibers,'keep');   
       
        fprintf('[%s] Extracting fiber density and length of the rejected fibers\n',mfilename)
        %results.rejected.density = dtiComputeFiberDensityNoGUI(fgB,xformimg2acpc,mapsize);
        results.rejected.length  = cellfun(@length,fgB.fibers); 
        results.rejected.n = sum(badFibers);
        fprintf('[%s] Saving a Rejected FG: \n%s\n',mfilename,fgBadFileName)
        fgWrite(fgB,fgBadFileName);
        clear fgB badFibers

        fprintf('[%s] Extracting the optimized FG\n',mfilename)
        fgG = fgExtract(fg,goodFibers,'keep');        
        clear fg
        fprintf('[%s] Extracting fiber density and length of the candidate FG\n',mfilename)
        %results.optimized.densityw = dtiComputeFiberDensityNoGUI(fgG,xformimg2acpc,mapsize,[],[],[],[],w);
        %results.optimized.density  = dtiComputeFiberDensityNoGUI(fgG,xformimg2acpc,mapsize);
        results.optimized.length   = cellfun(@length,fgG.fibers); 
        results.optimized.n = sum(goodFibers);
        clear goodFibers
        
        fprintf('[%s] Saving a Optimized FG: \n%s\n',mfilename,fgGoodFileName)
        fgWrite(fgG,fgGoodFileName);
        clear fgG
        
        fprintf('[%s] Saving a Results: \n%s\n',mfilename,resultsFileName2Save)
        mkdir(resultsSaveDir)
        save(resultsFileName2Save,'results')
        clear results w
        
    end
end


return

