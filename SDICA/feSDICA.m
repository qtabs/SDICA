function Xhat = feSDICA(filename, labels, n, Z, Q, Xnorm, kappa, subset)   
    % -----------------------------------------------------------------------
    % Xhat = feSDICA(filename, labels, n, Z, Q, Xnorm, kappa, subset)   
    % 
    % Returns a peak in the SDICA function
    % Inputs: 
    %         filename:   path to the file with the matrix (the matrix
    %                     should be contained in the file with the name
    %                     "X". The dimensionality reduction will be 
    %                     performed over the rows i.e. the final number
    %                     of cols will be always the same as the 
    %                     initial)
    %            labels:  vector with the labels of the samples
    %                 n:  number of basis elements to find
    %                 Z:  whitened and possibly reduced sample matrix 
    %                 Q:  Q-matrix
    %             Xnorm:  vector with the norms of the samples in X
    %             kappa:  relative importance of the SDICA components 
    %            subset:  subset of instances to use for the learning 
    %                     process (i.e. when finding the components).
    %                     Useful for cross-validation or statistical tests.
    %                     Note that this library will return the whole 
    %                     transformed dataset in any case. (Use -1 or left 
    %                     empty for using all instances)
    %
    % Output:
    %              Xhat:  matrix with the transformed samples
    % -----------------------------------------------------------------------

    if nargin < 7
        kappa = -1;
    end

    if nargin < 8 || ismember(-1, subset)
        subset = 1:M;
    end

    [N, M] = size(Q);

    solutions = SDICA(labels(subset), Z, n, Q(:, subset), kappa);

    
    fprintf('Projecting samples...');

    for i = 1:n
        [u, v] = transformSamples(Q, Z, Xnorm, solutions{i}.w);
        Xhat(:, 2 * i - 1) = u;
        Xhat(:, 2 * i) = v;
        fprintf('.');
    end

    fprintf('done! \n\n');



function [u, v] = transformSamples(Q, Z, Xnorm, w)
    % -----------------------------------------------------------------------
    % [u, v] = transformSamples(X, Z, w)
    % 
    % Returns the BD-transformed samples. Inputs:
    %              Q:   Q-matrix
    %              Z:   reduced sample matrix (rows samples, cols variables)
    %          Xnorm:   norm of the elements in the data matrix X
    %              w:   considered basis projector
    %
    % Outputs:
    %              u:   direct projection of the samples over the Ind. Comp.
    %              v:   the square of the u values
    % -----------------------------------------------------------------------

    u = w * Q;
    u = u .* Xnorm;
    v = u .^ 2;
