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

% List of subject IDs 
% (Semimanually extracted from the phenotopic data, not very elegant, but...)
subjects = {'1056121', '1113498', '1133221', '1139030', '1186237', ...
             '1240299', '1258069', '1282248', '1302449', '1391181', ...
             '1408093', '1469171', '1561488', '1686092', '1689948', ...
             '1791543', '1805037', '1875711', '1879542', '1912810', ...
             '1947991', '2081754', '2106109', '2123983', '2174595', ...
             '2196753', '2240562', '2249443', '2266806', '2367157', ...
             '2408774', '2427408', '2535087', '2538839', '2697768', ...
             '2703336', '2714224', '2833684', '2897046', '2910270', ...
             '3004580', '3086074', '3212536', '3233028', '3239413', ...
             '3262042', '3269608', '3306863', '3390312', '3554582', ...
             '3587000', '3593327', '3672854', '3707771', '3732101', ...
             '3739175', '3809753', '3889095', '3967265', '3976121', ...
             '3983607', '4028266', '4053836', '4091983', '4095748', ...
             '4256491', '4334113', '4383707', '4475709', '4921428', ...
             '5150328', '5193577', '5600820', '6187322', '7093319', ...
             '7135128', '7390867', '8328877', '8838009', '9093997', ...
             '9210521', '9221927', '9783279', '9887336', '9890726', ...
             '3655623', '3248920', '2659769', '9640133', '4055710', ...
             '1050975', '1809715', '1916266', '1117299', '2559537', ...
             '2296326', '7407032', '3157406', '3994098', '3561920', ...
             '2140063', '3610134', '3494778', '2033178', '1562298', ...
             '1093743', '1177160', '1860323', '1494102', '2498847', ...
             '1875013', '1068505', '3993793', '2310449', '3308331', ...
             '4265987', '3562883', '9578631', '2377207', '7011503', ...
             '8278680', '4075719', '2737106', '1159908', '4053388', ...
             '7689953', '2884672', '2601519', '2529026', '6500128', ...
             '1628610', '3691107', '4221029', '3856956', '1643780', ...
             '3910672', '3194757', '2141250', '9002207', '3446674', ...
             '2031422', '4073815', '1341865', '7253183', '5993008', ...
             '4225073', '1094669', '2950754', '2207418', '2919220', ...
             '3827352', '3205761', '1050345', '1132854', '1356553', ...
             '1399863', '1404738', '1411536', '1662160', '1771270', ...
             '1794770', '1843546', '2107404', '2208591', '2228148', ...
             '2268253', '2276801', '2493190', '2524687', '2780647', ...
             '2907951', '2940712', '2984158', '3224401', '3277313', ...
             '3291029', '3385520', '3473830', '3624598', '3672300', ...
             '3712305', '3803759', '3870624', '3930512', '4006710', ...
             '4048810', '4136226', '4241194', '5575344', '5669389', ...
             '6383713', '6477085', '7994085', '8191384'};

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

subjectVols = [];

for i = 1:length(subjects)
    
    % 1. Transform the NIFTIs into standard MATLAB vectors
    fprintf('\nTransforming subject %d/%d....\n', i, length(subjects)); 
    f = [NIFTIfold subjects{i} '.nii'];

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
    filename = [datafold 'mats/' subjects{i} '.mat'];
    save(filename, 'X', 'labels');

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
        filename = [datafold 'mats/' subjects{i} '.mat'];
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