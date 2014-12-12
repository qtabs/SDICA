% In this example we import an fMRI dataset to the SDICA data structure
% In our case we used the preprocessed peking subsample of the ADHD200:
% http://www.nitrc.org/plugins/mwiki/index.php/neurobureau:AthenaPipeline

% The importation of the NIFTI volumes requires SPM:
% http://www.fil.ion.ucl.ac.uk/spm/software/spm8/


addpath('./import')
% You will need some SPM thingies to import the NIFITs. If you haven't already
% done it, add the SPM path to MATLAB:
addpath(genpath('~/Apps/MatlabApps/spm8/'));


% -- Definitions --

% Folder with the NIFTI files
NIFTIfold = './data/NIFTIs/';
datafold = './data/';

% Amount of samples in the Z-set
N = 1000;

% Limits on the amount of volumes kept for each subject in the PCA:
minNumOfComp = 10;
maxNumOfComp = 100;

% Load the list of subject IDs and path from a txt file 
filename = 'path_subjects_ADHD_ExtDisk_subsetPekin.txt';
subjects = loadFromFile(filename);


% Labels of the data
DX = [0, 0, 1, 0, 3, 1, 0, 3, 0, 0, 0, 0, 3, 0, 0, 3, 0, 0, 0, 0, 3, ...
      0, 0, 0, 1, 3, 0, 0, 0, 1, 0, 0, 0, 0, 3, 0, 0, 0, 3, 3, 0, 0, ...
      0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 3, 0, ...
      3, 1, 0, 3, 0, 0, 0, 3, 3, 0, 0, 0, 0, 3, 0, 0, 0, 3, 0, 0, 0, ...
      0, 0, 0, 0, 0, 3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, ...
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 1, 3, 1, 1, 3, ...
      3, 3, 3, 3, 1, 1, 1, 3, 1, 1, 3, 3, 1, 3, 3, 3, 1, 3, 3, 3, 1, ...
      1, 1, 3, 3, 1, 0, 0, 1, 1, 0, 0, 0, 3, 0, 3, 0, 0, 1, 0, 3, 0, ...
      3, 1, 3, 3, 3, 0, 0, 3, 0, 0, 0, 3, 3, 3, 1, 0, 3, 0, 0, 1, 0, ...
      0, 1, 0, 0, 0];

% Subset of the samples considered in the analysis (0 to exclude)
% (Useful to produce balanced subdatasets)
subs = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ...
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ...
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ...
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ...
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ...
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ...
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ...
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ...
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ...
        1, 1, 1, 1, 1];
    
% Transform into a 2-classes problem (all ADHD kinds belong to the same class)
classes = ((DX == 0) * -2) + 1;



% --- Importing the fMRI volumes to a Matlab matrix ---
mkdir([datafold 'mats/']);
subjectVols = [];

for i = 1:length(subjects)
    if subs(i) == 1
        % 1. Transform the NIFTIs into standard MATLAB vectors
        fprintf('\nTransforming subject %d/%d....\n', i, length(subjects)); 
        f = subjects{i};

        fprintf('-- Loading NIFTI file...')
        x = loadfMRI(f);
        label = classes(i);

        fprintf(' merging volumes...')
        X = [];
        for j = 1:length(x)
            X = [X; x{j}];
        end

        [M, D] = size(X);
        fprintf(' Samples: %d, Dimension: %d \n', M, D);

        % 2. Subject specific preprocessing: PCA
        fprintf('-- Preprocessig...\n')
        % In this case we change the maxNumOfComp = lastEig:
        if maxNumOfComp > size (X, 1)       
            [E, D] = pcamat(X, minNumOfComp, size (X, 1), 'off', 'on'); 
        else
            [E, D] = pcamat(X, minNumOfComp, maxNumOfComp, 'off', 'on');         
        end
         X = (inv(sqrt(D)) * E') * X;

        [M, D] = size(X);
        fprintf('......... Samples: %d, Dimension: %d\n', M, D);

        % 3. Label each volume as an individual sample
        labels = repmat(classes(i), [1 M]);
        subjectVols = [subjectVols; repmat(i, [M, 1])];

        % 4. Save the preprocessed subject as a MATLAB matrix
        fprintf('-- Saving file.... \n\n');
        [pth, name_subj, ext] = fileparts(subjects{i});
        filename = [datafold 'mats/' name_subj '.mat'];
        save(filename, 'X', 'labels');        
    else
        fprintf('..not in subset!\n');
    end
end

% Useful information to keep for performing subject-like analysis:
% save('volsToSubjectData.mat', 'subjects', 'subjectVols', 'DX');


% 5. Merge all the subject-specific volumes into a single matrix
fprintf('\nGenerating dataset file...\n')

X = [];
labels = [];
Xfilename = [datafold 'X.mat'];
save(Xfilename, 'X', 'labels', '-V7.3');

mfile = matfile(Xfilename, 'Writable', true);
st = 1;

for i = 1:length(subjects)
    fprintf('-- Loading subject %d/%d....', i, length(subjects));
    if subs(i) == 1
        [pth, name_subj, ext] = fileparts(subjects{i});
        filename = [datafold 'mats/' name_subj '.mat'];
        load(filename);
        [M, D] = size(X);
        nd = st + M - 1;
        mfile.X(st:nd, 1:D) = X;
        mfile.labels(1, st:nd) = repmat(classes(i), [1 M]);
        fprintf('..done!\n');
        st = nd + 1;
    else
        fprintf('..not in subset!\n');
    end
end
fprintf('-- Process finished! <3\n\n')



% -- Computing extra matrices for the SDICA analysis --

% Most of the import functions work by dividing the set in several 
% parts to avoid memory problems. If you feel like it is taking too long or
% you still have memory issues computing your matrices, adjust the number of 
% subdivisions in the import scripts.

% 1. Compute the Z-set 
% (eigenvalues l can be used to further reduce the size of Z)
[Z, l] = computeZ(Xfilename, N);
% You can reduce the size of Z to reduce the searc-space of the algorithm
Z = Z(1:N, :); 
save([datafold 'Z.mat'], 'Z');

% 2. Compute the norm of the samples
Xnorm = computeNorms(Xfilename);
save([datafold 'Xnorm.mat'], 'Xnorm');

% 3. Compute the Q-matrix (this can take a while...)
Q = computeQ(Xfilename, Z);
save([datafold, 'Q.mat'], 'Q');


% That's it!