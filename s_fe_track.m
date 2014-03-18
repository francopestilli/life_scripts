function s_fe_track()
%
% This function:
%  - Tracks using mrtrix accessign a few subjects at the time
%  - Generates an optimized connectome from a cadidate connectome using 
%  LIFE method
%
%  fe = s_fe_track()
% 
% INPUTS:  none
% OUTPUTS: fe structure the optimized life structure
%
% Copyright Franco Pestilli (2013) Vistasoft Stanford University.


% Get the base directory for the data
datapath = '/marcovaldo/frk/2t1/predator/';
subjects = {'MP_96dirs_b2000_1p5iso','HT_96dirs_b2000_1p5iso'};
%           'FP_96dirs_b2000_1p5iso','JW_96dirs_b2000_1p5iso'};%, ...
%            'KK_96dirs_b2000_1p5iso','KW_96dirs_b2000_1p5iso', ...
%            'MP_96dirs_b2000_1p5iso','HT_96dirs_b2000_1p5iso'};
nSeeds            = 500000;
trackingAlgorithm = {'tensor','prob','stream'};
lmax               = [2 8 10];

for isbj = 1:length(subjects)
    % Build the file names for the diffusion data, the anatomical MR, the fiber
    % group containing the connectome and the
    baseDir                     = fullfile(datapath,subjects{isbj});
    dtFile.files.alignedDwRaw   = fullfile(baseDir,'diffusion_data','run01_fliprot_aligned_trilin.nii.gz');
    dtFile.files.alignedDwBvecs = fullfile(baseDir,'diffusion_data','run01_fliprot_aligned_trilin.bvecs');
    dtFile.files.alignedDwBvals = fullfile(baseDir,'diffusion_data','run01_fliprot_aligned_trilin.bvals');
    dtFile.files.brainMask      = fullfile('diffusion_data','brainMask.nii.gz');
    dtFile.files.wmMask         = fullfile('anatomy','wm_mask.nii.gz');
    
    % Directory where to save the fe structures
    fibersFolder       = fullfile(baseDir,'fibers');
    
    % Track
    feTrack(trackingAlgorithm, dtFile,fibersFolder,lmax,nSeeds,fullfile(baseDir,dtFile.files.wmMask))
end

end

