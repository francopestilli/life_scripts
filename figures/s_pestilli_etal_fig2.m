function s_pestilli_etal_fig2

% Load the data from the Human Connectome Project
datapath = '/home/frk/2t2/HCP/';
subject = '115320';
diffusionData = fullfile(datapath,subject,'diffusion_data','dwi_data_b2000.nii.gz');
fprintf('[%s] loading diffusion data...',mfilename)
dwi = dwiLoad(diffusionData);

% This is the coordinate we will use for the example
coord = [67 70 62]; % CC
coord = [83 124 62]; % frontal
coord = [81 119 62]; % frontal
fibers{1} = [0.15 0.15 0.1; 0.75 0.5 0.95; ]';

% Load the bvecs and plot them
bvecs = dwiGet(dwi,'diffusion bvecs');
dwiPlot(dwi,'bvecs');

% Load the diffusion distance, this is the actual data
dDist = dwiGet(dwi,'diffusiondistanceimage',coord,'um');

% Generate the model tensor from the fibers in this voxel
dParms(1) = 10; 
dParms(2) = 1; 
dParms(3) = 1;

nFibers = length(fibers); % The number of Fibers.
Q       = cell(1,nFibers); % Memory for the tensors of each fiber.
D       = diag(dParms);    % The diagonal form of the Tensors' model parameters.

for ii = 1:nFibers
 % Compute the diffusion gradient at each node of the fiber.
 fiberGradient = gradient(fibers{ii});
 
 % Number of nodes fro this fiber
 numNodes = size(fibers{ii},2);
 
 % preallocated memory for the vector representation of tensors.
 T = zeros(numNodes,3,3);
 
 for jj = 1:numNodes
  % Rotate the tensor toward the gradient of the fiber.
  %
  % Calculate a rotation matrix for the tensor so that points in the fiberGradient
  % direction and has two perpendicular directions (U)
  % Leaving the 3 outputs for this function is the fastest use of it.
  [Rot,~, ~] = svd(fiberGradient(:,jj)); % Compute the eigen vectors of the gradient.
  
  % Create the quadratic form of the tensor.
  %
  % The principal eigenvector is in the same direction of the
  % fiberGradient. The direction of the other two are scaled by dParms.
  % Human friendly version fo the code:
  % tensor = Rot*D*Rot'; % tensor for the current node, 3x3 matrix.
  % T(jj,:) = reshape(tensor,1,9); % reshaped as a vector 1,9 vector
  T(jj,:,:) = Rot*D*Rot';
 end
 % T is a matrix; each row is a 1,9 version of the tensor.
 tensor{ii} = T;
end

% This should be a function like dtiPlotDist(tensor)
%
[X,Y,Z] = sphere(50);
[r,c] = size(X);

v = [X(:),Y(:),Z(:)];

iFib = 1;
for iNode = 1:size(tensor{1},1)
    adcPredicted = diag(v*squeeze(tensor{1}(iNode,:,:))*v');
    
    % The diffusion distance is the length of the vector, v, such
    % that v' tensor v = 1.  We know that for unit length vectors, the
    % adc = u' tensor u.  So, v = u / sqrt(adc).
    v = diag(1./sqrt(adcPredicted))*v;
    x = reshape(v(:,1),r,c);
    y = reshape(v(:,2),r,c);
    z = reshape(v(:,3),r,c);
    
    if iFib == 1
        % Compute and plot vectors of measured distances
        dDistV = diag(dDist)*bvecs;
        plot3(dDistV(:,1),dDistV(:,2),dDistV(:,3),'k.','MarkerSize',30)
        grid off
        box on
        axis equal,
        hold on
        cmap = autumn(255);
    end
    
    surf(fibers{iFib}(1,iNode)+x,fibers{iFib}(2,iNode)+y,fibers{iFib}(3,iNode)+z,repmat(256,r,c),'EdgeAlpha',0.1);
    colormap([cmap; .25 .25 .25]), alpha(0.25)
    camlight; lighting phong; material shiny;
    set(gca, 'Projection', 'perspective');
end
