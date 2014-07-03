function s_ms_hcp_track
%
% This function creates whole-brain white-matter fiber groups for the HCP data in:
%   /home/frk/2t2/HCP
%
% Copyright Franco Pestilli (2013) Vistasoft Stanford University.
datapath = {'/marcovaldo/frk/2t1/HCP','/marcovaldo/frk/2t2/HCP'};
subjects = {{'105115','110411','111312','113619'}, {'115320','117122','118730'}};
data_dir = 'diffusion_data';
anatomypath = getenv('SUBJECTS_DIR');

bval = [2000,1000,3000];
    
% Tracking parameters
trackingAlgorithm = {'tensor','stream','prob'};
lmax         = [10,2,8,12,6];
nSeeds       = 500000;
for ib = 1:numel(bval)
    for it = 1:numel(trackingAlgorithm)
        for id = 1:numel(datapath)
            for is = 1:numel(subjects{id})
                dt_info.files.alignedDwRaw   = fullfile(datapath{id},subjects{id}{is},data_dir,sprintf('dwi_data_b%i_aligned_trilin.nii.gz',bval(ib)));
                dt_info.files.alignedDwBvecs = fullfile(datapath{id},subjects{id}{is},data_dir,sprintf('dwi_data_b%i_aligned_trilin.bvecs', bval(ib)));
                dt_info.files.alignedDwBvals = fullfile(datapath{id},subjects{id}{is},data_dir,sprintf('dwi_data_b%i_aligned_trilin.bvals', bval(ib)));
                dt_info.files.brainMask      = fullfile(sprintf('dt6_b%itrilin', bval(ib)),'bin','brainMask.nii.gz');
                dt_info.files.wmMask         = fullfile(anatomypath,subjects{id}{is},'wm_mask.nii.gz');
                wmMask                       = fullfile(anatomypath,subjects{id}{is},'wm_mask.nii.gz');
                outputFibersFolder           = fullfile(datapath{id},subjects{id}{is},'fibers');
                
                switch trackingAlgorithm{it}
                    case 'tensor'
                        feTrack({trackingAlgorithm{it}}, dt_info ,outputFibersFolder,2,nSeeds,wmMask);
                        
                    otherwise
                        for il = 1:numel(lmax)
                            feTrack({trackingAlgorithm{it}}, dt_info ,outputFibersFolder,lmax(il),nSeeds,wmMask);
                        end
                end
            end
        end
    end
end

return
