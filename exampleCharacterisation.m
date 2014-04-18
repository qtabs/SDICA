addpath('./SDICA/');
addpath('./import');


% -- Data (matrices generated using the import libs) --

xfile = './data/X.mat';
qfile = './data/Q.mat';
zfile = './data/Z.mat';
normsfile = './data/Xnorm.mat';
resultsFold = './results/';
mfile = matfile(xfile);
labels = mfile.labels;
% Need to convert the problem into a two-classes problem?
% labels = 2 * (not(mfile.labels == 0) - 0.5);

% You'll need to supply one of the original NIFITs of the data for reference
% to save the obtained maps in a NIFIT format
referenceNIFTI = './data/NIFTIs/1585708.nii';


% -- Parameters --

% Number of components to extract
n = 2;
% Number of elements in Z
N = 500;
% Value of the weighting coefficient kappa
kappa = 0.75;
		
load(qfile);
load(zfile);
load(normsfile);

Z = Z(1:N, :);
Q = Q(1:N, :);


% -- Analysis --
solutions = SDICA(labels, Z, n, Q, kappa);


% -- Writting down the solutions --
for i = 1:length(solutions)
	map = solutions{i}.w * Z;
	targetPath = [resultsFold, 'map', sprintf('%02d', i), '.nii'];
	writefMRI(map, referenceNIFTI, targetPath);
end