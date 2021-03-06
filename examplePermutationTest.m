% Run the permutationTest function with different parameters and subsets of
% data:

addpath('./SDICA/');
addpath('./validation/');

% Load labels and data
    
xfile = './data/X.mat';
mfile = matfile(xfile);
labels = mfile.labels;

qfile = './data/Q.mat';
zfile = './data/Z.mat';
load(qfile);
load(zfile);

% Number of elements in Z
N=500;

Z = Z(1:N, :);
Q = Q(1:N, :);

% Set parameters:
kappa = 0.9;
n = 5;
nOfPerm = 2500;

tic 
tmaps = permutationTest(labels, Z, Q, kappa, n, nOfPerm);
toc

save('tmaps_K09_n5_Perm2500','tmaps');

