function solutions = SDICA(labels, Z, n, Q, kappa)
    % -----------------------------------------------------------------------
    % solutions = SDICA(labels, Z, n, Q)
    % 
    % Returns a set of n solutions of SDICA for the problem {Z, Q, labels}
    %
    % Inputs: 
    %            labels:  vector with the labels of the samples
    %                 Z:  whitened and possibly reduced sample matrix 
    %                 n:  number of basis elements to find
    %                 Q:  Q-matrix
    %             kappa:  Relative importance of negentropy and SDICA [0-1]
    %
    % Output:
    %         solutions: an ordered array of n structs containing the 
    %                            solution.w:  a maxima of the SDICA function
    %                                         (as a col vector)
    %                          solution.neg:  negentropy at that point
    %                         solution.sfld:  SFLD score at that point
    %                            solution.J:  joint objective function value
    %                        solution.index:  the place in which the solution 
    %                                         was found by the algorithm
    % -----------------------------------------------------------------------

    % -- HARDCODED PARAMETERS --
    % Relative importance of negentropy and SDICA (0 < kappa < 1; def. 0.7)
    kappaDefault = 0.7;   
    
    % Learning rate (eta > 0; def. 0.15)
    hcPars.eta = 0.15;
    % Importance of the inertia term (mass >= 0; def. 0.1)
    hcPars.mass = 0.1;

    % Tolerance for the stopping criterion (def. 1E-5)
    hcPars.tol = 1E-5;
    % Max. number of iterations (maxIt > 0; def. 99)
    hcPars.maxIt = 99;
    % Min. number of iterations (minIt >= 0; def. 5)
    hcPars.minIt = 2;

    % Decay of the learning rate (0 < dec <= 1; def. 1)
    hcPars.dec = 0.99;
    % Decay of kappa (0 < dec <= 1; def. 1)
    hcPars.kappaDecay = 1;
    % Decay of the inertia multiplier (0 < dec <= 1; def. 1)
    hcPars.massDecay = 0.99;

    % Parameters for the intialisation steps
    initPars.kappa = kappa; % 0.4;
    initPars.eta = 3 * hcPars.eta;
    initPars.mass = 0.5 * hcPars.mass;
    % Number of different random itialisations to test (nOfTries >= 1, def. 1)
    initPars.nOfTries = 1;
    % Number of iteration in the initialisation (its >= 0, def. 5)
    initPars.its = 1;

    if nargin == 4 | kappa < 0
        fprintf('Using Default kappa....\n')
        kappa = kappaDefault;
    end


    % Initialising search variables...
    pastWs = {};

    % Looking for the n discriminant networks
    for i = 1:n
        fprintf('Looking for component %2d/%2d...\n', i, n);
        % initPars.nOfTries = ceil(initPars.nOfTries / i);
        [w0, delta0] = initialisation(labels, Z, Q, initPars, pastWs);
        s{i} = hillClimbing(labels, Z, Q, w0, delta0, kappa, hcPars, pastWs);
        pastWs{i} = s{i}.w;
        J(i) = s{i}.J;
    end

    % Ordering solutions
    % [drop, indices] = sort(J, 'descend');
    indices = 1:n;
    for i = 1:n
        solutions{i} = s{indices(i)};
        solutions{i}.index = i;
    end

    save(['./K', int2str(kappa * 10), 'sols.mat'], 'solutions');



function solution = hillClimbing(labels, Z, Q, w0, delta0, kappa, p, pastWs)
    % -----------------------------------------------------------------------
    % solution = hillClimbing(X, labels, Z, Q, w0, delta0, kappa, p, pastWs)
    % 
    % Returns a peak in the SDICA function
    % Inputs: 
    %            labels:  vector with the labels of the samples
    %                 Z:  whitened and possibly reduced sample matrix 
    %                 Q:  Q-matrix 
    %                w0:  a proposed initialisation for hill climbing
    %            delta0:  a initialisation for the original velocity (e.g. 0)
    %             kappa:  relative importance of negentropy and SFLD in the
    %                     search: J = (1-kappa) * negentropy + kappa * SFLD
    %                 p:  a struct containing the search parameters:
    %                                 p.eta: learning rate
    %                                p.mass: inertial mass (0 for no inertia)
    %                                 p.tol: difference between ws for stop
    %                               p.maxIt: maximum number of iterations
    %                               p.minIt: minimum number of iterations
    %                                 p.dec: decay for eta (1 for no decay)
    %                          p.kappaDecay: decay for kappa (1 for no decay)
    %                           p.massDecay: decay for mass (1 for no decay)
    %            pastWs:  an array containing all the already extracted ws
    %
    % Output:
    %          solution: a struct containing the solution:
    %                            solution.w:  a maxima of the SDICA function
    %                                         (as a col vector)
    %                          solution.neg:  negentropy at that point
    %                         solution.SFLD:  SFLD at that point
    %                            solution.J:  joint objective function value
    % -----------------------------------------------------------------------

    fprintf('\nLooking for maxima......\n');

    [N, dim] = size(Z);

    % Initialising search variables...
    t = cputime;
    eta = p.eta;
    mass = p.mass;
    theta = pi;
    it = 0;

    % Hill climbing cycle
    while (theta > p.tol) & (it < p.maxIt) % | (it < p.minIt)

        % Tracking iterations...
        it = it + 1;
        fprintf('-- It %3d:', it);

        % Measuring contributions to the gradient...
        if not(kappa == 1)
            [Neg, deltaInd] = deltaICA(Z, w0);
        else
            Neg = 0;
            deltaInd = 0;
        end

        if not(kappa == 0)
            [SFLD, deltaDisc] = deltaSFLD(Q, labels, w0);
        else
            SFLD = 0; 
            deltaDisc = 0;
        end

        % Computing delta...
        eta = eta * p.dec;
        mass = mass * p.massDecay;
        kappa = kappa * p.kappaDecay;
        delta = (1 - kappa) * deltaInd + kappa * deltaDisc;
        delta = delta + mass * delta0;

        % Computing new projector vector...
        w = w0 + eta * delta;
        w = orthoNormForcing(w, pastWs);
        w = w / norm(w);

        % Printing information about the cycle...
        theta = 1 - abs(w0 * w');
        J = (1 - kappa) * Neg + kappa * SFLD;
        fprintf(' theta = %.5f', theta/pi * 180);
        fprintf(', Neg = %.3f, SFLD = %.3f; J = %.3f \n', Neg, SFLD, J);

        % Initialising initial values for next iteration...
        w0 = w;
        delta0 = delta;

    end

    % Measuring and printing results...
    t = cputime - t;

    fprintf('  ...done! ---> #Its = %d, Time = %.0fs', it, t);

    [neg,  drop] = deltaICA(Z, w);
    [disc, drop] = deltaSFLD(Q, labels, w);
    J = (kappa - 1) * neg + kappa * disc;

    fprintf(', J = %.3f \n\n', J);

    solution.w = w;
    solution.J = J;
    solution.neg = neg;
    solution.sfld = disc;



function [w0, delta0] = initialisation(labels, Z, Q, p, pastWs)
    % -----------------------------------------------------------------------
    % initialisation(labels, Z, p, pastWs)
    % 
    % Returns a intialisation seed
    % Inputs: 
    %            labels:  vector with the labels of the samples
    %                 Z:  whitened and possibly reduced sample matrix 
    %                 Q:  Q-matrix 
    %                 p:  a struct containing the search parameters:
    %                                 p.eta: learning rate
    %                                p.mass: inertial mass (0 for no inertia)
    %                            p.nOfTries: number of random initialisations
    %                                        tried
    %                                 p.its: number of iterations for each 
    %                                        tested random w 
    %            pastWs:  an array containing all the already extracted ws
    %
    % Output:
    %                w0:  a proposed initialisation for hill climbing
    %            delta0:  the last delta in the rough search
    % -----------------------------------------------------------------------

    fprintf('\n\nInitialising search vector...\n')

    epsilon = 0;         % To avoid rejecting good initial position
    N = size(Z, 1);

    % Initialising variables
    t = cputime;
    maxJ = -inf;
    alreadyTriedWs = {};

    % Running initialisations...
    for i = 1:p.nOfTries
    
        wsToAvoid = {pastWs{:}, alreadyTriedWs{:}};

        w0 = rand([1 N]);
        w0 = orthoNormForcing(w0, wsToAvoid);

        fprintf('-- Trying %1d of %1d...\n', i, p.nOfTries);

        [w, delta0, J] = roughSearch(labels, Z, Q, w0, p, wsToAvoid);

        % Storing greater results...
        if J > (maxJ + epsilon)
            maxJ = J;
            bestW = w;
            bestDelta = delta0;
        end
    
        % Tracking previous attempts to improve initialisation performance...
        alreadyTriedWs{i} = w0;

    end

    % Measuring and printing results...
    t = cputime - t;

    if p.nOfTries > 0
        w0 = bestW;
        delta0 = bestDelta;
        fprintf('...done! Best J: %.3f, Time = %.0fs\n', maxJ, t);
    else 
        w0 = rand([1 N]);
        w0 = orthoNormForcing(w0, wsToAvoid);
        delta0 = 0;
        fprintf(' Random initialisation done!\n');
    end



function [w0, delta0, J] = roughSearch(labels, Z, Q, w0, p, pastWs)
    % -----------------------------------------------------------------------
    % [w0, delta0, SFLD] = roughSearch(X, labels, Z, w0, p, pastWs)
    % 
    % Runs hill climbing using only the SFLD. Aimed for the initialisation
    % of the SDICA hill climbing algorithm
    %
    % Inputs: 
    %            labels:  vector with the labels of the samples
    %                 Z:  whitened and possibly reduced sample matrix 
    %                 Q:  Q-matrix 
    %                w0:  a proposed initialisation
    %                 p:  a struct containing the search parameters:
    %                                 p.eta: learning rate
    %                                p.mass: inertial mass (0 for no inertia)
    %                                 p.its: a fixed number of iterations
    %            pastWs:  an array containing all the already extracted ws
    %
    % Output:
    %                w0:  a proposed initialisation for hill climbing
    %            delta0:  the last delta in the rough search
    %                 J:  value of the J at the proposed initialisation
    % -----------------------------------------------------------------------

    delta0 = 0;

    for i = 1:p.its
        fprintf('---- It %1d of %1d:', i, p.its);

        % Measuring contributions to the gradient...
        [SFLD, deltaDisc] = deltaSFLD(Q, labels, w0);
        [Neg, deltaInd] = deltaICA(Z, w0);

        delta = (1 - p.kappa) * deltaInd + p.kappa * deltaDisc;
        delta = delta + p.mass * delta0;

        % Computing new projector vector...
        w = w0 + p.eta * delta;
        % w = orthoNormForcing(w, pastWs);
        w = w / norm(w);

        % Verbose information...
        theta = 1 - abs(w0 * w');
        J = (1 - p.kappa) * Neg + p.kappa * SFLD;
        fprintf(' theta = %.5f', theta/pi * 180);
        fprintf(', Neg = %.3f, SFLD = %.3f; J = %.3f \n', Neg, SFLD, J);


        % Initialising initial values for next iteration...
        w0 = w;
        delta0 = delta;
    end



function w = orthoNormForcing(w, pastWs)
    % -----------------------------------------------------------------------
    % w = orthoNormForcing(w, pastWs)
    % 
    % Removes the contributions of the already found basis projection vectors
    % in pastWs of the given projection vector and performs normalisation. As
    % a result, te set {pastWs, w} is an orthonormal set of vectors.
    %
    % Inputs: 
    %                 w:  given projection vector w
    %            pastWs:  an array containing a set of orthonormal ws
    %
    % Output:
    %                 w:  a representation of the input w such that the set
    %                     {w, pastWs} is now an orthonormal set of vectors
    % -----------------------------------------------------------------------

    w = w / norm(w);

    for i = 1:length(pastWs)
        w = w - (w * pastWs{i}') * pastWs{i};
    end

    w = w / norm(w);

