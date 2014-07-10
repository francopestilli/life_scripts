function s_fe_compute_length_hist(bval)
%
% This function:
%  - Load FG strucutres
%  - Compute fiber statistics (length, etc) 
%
%  fe = s_fe_compute_length_hist()
% 
% Copyright Franco Pestilli (2014) Vistasoft Stanford University
keyboard
% DO NOT USE THIS SCRIPT
% Get the base directory for the data
if ~isempty(bval)
    datapath = '/marcovaldo/frk/2t2/predator/';
    subjects = {...
            'FP_150dirs_b1000_2000_4000_2iso', ...
            };
else
    datapath = '/marcovaldo/frk/2t1/predator/';
    subjects = {...
            'KK_96dirs_b2000_1p5iso', ...           
            'MP_96dirs_b2000_1p5iso', ...
            'JW_96dirs_b2000_1p5iso', ...
            'HT_96dirs_b2000_1p5iso', ...
            'KW_96dirs_b2000_1p5iso', ...
            'FP_96dirs_b2000_1p5iso', ...
            };
end

for isbj = 1:length(subjects)
    % Directory where to save the fibers and the results
    fibersDir       = fullfile(datapath,subjects{isbj},'fibers');
          resultsSaveDir       = fullfile(datapath,subjects{isbj},'results');

    if ~isempty(bval)
       fgFiles       = dir(fullfile(fibersDir,sprintf('*%s*.mat',num2str(bval))));
    else
       fgFiles       = dir(fullfile(fibersDir,sprintf('*lmax10*.pdb')));
    end
  keyboard          
    % We build one modelper fiber group, whole brain fiber group
    for iFe = 1:length(fgFiles)
        % The final connectome and dat astructure will be saved with this name:
        [~,feFileName] = fileparts(fgFiles(iFe).name);
        resultsToLoad = fullfile(resultsSaveDir,[feFileName,'-fiberStatsResults.mat']);

        % Initialize the Connectome
        fprintf('[%s] Loading RESULTS: \n%s\n',mfilename,resultsToLoad)
        load(resultsToLoad);
        
        fprintf('[%s] Extracting fiber density and length of the candidate FG\n',mfilename)
        xbins = [1,2,4,8,16,32,64,128,256,512];
        [c.y(isbj,:),c.x(isbj,:)] = hist(results.candidate.length,xbins);
        [o.y(isbj,:),o.x(isbj,:)] = hist(results.optimized.length,xbins);
          
        
        fprintf('[%s] Extracting the optimized FG\n',mfilename)
        fgB = fgExtract(fg,find(badFibers),'keep');   
       
        fprintf('[%s] Extracting fiber density and length of the rejected fibers\n',mfilename)
        %results.rejected.density = dtiComputeFiberDensityNoGUI(fgB,xformimg2acpc,mapsize);
        results.rejected.length  = cellfun(@length,fgB.fibers); 
        results.rejected.n = sum(badFibers);
        fprintf('[%s] Saving a Rejected FG: \n%s\n',mfilename,fgBadFileName)
        fgWrite(fgB,fgBadFileName);
        clear fgB badFibers

        fprintf('[%s] Extracting the optimized FG\n',mfilename)
        fgG = fgExtract(fg,find(goodFibers),'keep');        
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
        
        fprintf('[%s] Saving a Results: \n%s\n',mfilename,resultsToLoad)
        mkdir(resultsSaveDir)
        save(resultsToLoad,'results')
        clear results w
        
    end
end


return

