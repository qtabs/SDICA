function tmaps = permutationTest(labels, Z, Q, kappa, n, nOfPerm)
    % -----------------------------------------------------------------------
    % maps = permutationTest
    % 
    % Returns n t-maps representing the components extracted with SDICA. The 
    % t-maps are obtained in a non-parametric way using a random permutation 
    % test
    %
    % Inputs: 
    %            labels:  vector with the labels of the samples
    %                 Z:  whitened and possibly reduced sample matrix 
    %                 n:  number of basis elements to find
    %                 Q:  Q-matrix
    %             kappa:  Relative importance of negentropy and SDICA [0-1]
    %           nOfPerm:  Number of permutations used in the test
    %
    % Output:
    %             tmaps: an ordered array of n structs containing the 
    %                 tmaps{i}.component: the original component
    %                      tmaps{i}.tmap: associated t-map
    %                    tmaps{i}.ppvals: map of p-values for the t > 0
    %                    tmaps{i}.npvals: map of p-values for the t < 0
    %                       tmaps{i}.neg: negentropy of the solution
    %                      tmaps{i}.sfld: SFLD of the solution
    %                         tmaps{i}.J: joint objective function value
    % -----------------------------------------------------------------------

    if nargin < 6, nOfPerm = 2000; end;
    if nargin < 5, n = 20; end; 
    if nargin < 4, kappa = 1; end; 

    % First we extract the maps using SDICA
    solutions = SDICA(labels, Z, n, Q, kappa);

    % Then we estimate the baseline map distribution using permTest
    for i = 1:nOfPerm
        permutedLabels = labels(randperm(length(labels))); 
        s = SDICA(permutedLabels, Z, 1, Q, kappa); 
        voxDist(i, :) = s{1}.w * Z;
    end
    
    voxSigma = std(voxDist);

    % Now we can normalise the voxels in the maps to a t-measure
    for i = 1:length(solutions)
     
        tmaps{i}.component = solutions{i}.w * Z;
        tmaps{i}.tmap = (tmaps{i}.component - mode(tmaps{i}.component)) ./ voxSigma;
       
        tmaps{i}.ppvals = 1 - tcdf(tmaps{i}.tmap, nOfPerm - 1);
        
        tmaps{i}.npvals = tcdf(tmaps{i}.tmap, nOfPerm - 1);
        
        tmaps{i}.negentropy = solutions{i}.neg;
        tmaps{i}.sfld = solutions{i}.sfld;
        tmaps{i}.J = solutions{i}.J;
    
    end

end