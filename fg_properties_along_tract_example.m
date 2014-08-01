


direction = 'AP'; % Specify direction in which paths should be ordered
%direction = 'LR';

subDir = fullfile('somedirectory','dt6.mat');
dt     = dtiLoadDt6(fullfile(subDir));

fiberDir   = fullfile('somedirectorytothefibers', 'myfibers.mat/pdb');
fg_clean   = fgRead(fiberGroup);

% Align fibers across subjects
fg_clean = dtiAlignFiberDirection(fg_clean,direction);

% extract fiber properties
[fa,md,rd,ad,cl,SuperFibersGroup]=dtiComputeDiffusionPropertiesAlongFG(fg_clean, dt);

