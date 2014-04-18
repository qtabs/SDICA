function [SFLD, grad] = deltaSFLD(Q, labels, w)
    % -----------------------------------------------------------------------
    % [SFLD, grad] = deltaSFLD(Z, labels, w)
    % 
    % Returns the gradient of the Spatial-FLD and the measure 
    % itself for a given w.
    % Inputs:
    %               Q:   Q-matrix
    %          labels:   vector with the labels of the samples
    %               w:   considered basis projector (col vector)
    %
    % The output is:
    %            SFLD:   The basis-decomposition FLD for the given w
    %            grad:   vector of partial derivatives of the SFLD 
    %                    respect with the components of w (col vector)
    % -----------------------------------------------------------------------

    [xi, eta] = projectSamples(Q, w);

    [xid, etad] = projectionDerivatives(Q, w, xi);

    [SW, SB] = scatterMatrices(xi, eta, labels);

    [SWd, SBd] = scatterDerivatives(xi, eta, labels, xid, etad);

    SFLD = logFLD(SW, SB);
    grad = logSFLDd(SB, SW, SBd, SWd);



function lBCd = logSFLDd(SB, SW, SBd, SWd)
    % -----------------------------------------------------------------------
    % lBDd = logSFLDd(SB, SW, SBd, SW)
    % 
    % Returns the partial derivatives of the log-basis-decomposition FLD
    % Inputs:
    %               SB:    Diagonal elements of the Between scatter matrix 
    %                      in vector form [11, 22]
    %               SW:    Diagonal elements of the Within scatter matrix 
    %                      in vector form [11, 22]
    %              SBd:    Partial Derivatives of the diagonal elements of 
    %                      the Between scatter matrix
    %              SWd:    Partial Derivatives of the diagonal elements of 
    %                      the Within scatter matrix
    % -----------------------------------------------------------------------

    trSB = SB(1) + SB(2);
    trSW = SW(1) + SW(2);

    trSB = max(trSB, eps);
    trSW = max(trSW, eps);

    trSBd = SBd(:, 1) + SBd(:, 2);
    trSWd = SWd(:, 1) + SWd(:, 2);

    lBCd = trSBd' / trSB - trSWd' / trSW;



function fld = logFLD(SW, SB)
    % -----------------------------------------------------------------------
    % fld = logFLD(SW, SB)
    %
    % Performs log-Fisher-Linear-Discriminant given the transformed
    % scatter matrices
    % Inputs:
    %               SB:    Diagonal elements of the Between scatter matrix 
    %                      in vector form [11, 22]
    %               SW:    Diagonal elements of the Within scatter matrix 
    %                      in vector form [11, 22]
    % -----------------------------------------------------------------------

    trSB = SB(1) + SB(2);
    trSW = SW(1) + SW(2);

    trSB = max(trSB, eps);
    trSW = max(trSW, eps);

    fld = log(trSB) - log(trSW);



function [xi, eta] = projectSamples(Q, w)
    % -----------------------------------------------------------------------
    % [xi, eta] = projectSamples(Z, w)
    % 
    % Returns the xi-eta representation of the samples. Inputs:
    %              Q:   Q-matrix
    %              w:   considered basis projector
    %
    % Outputs:
    %             xi:   direct projection of the samples over the Ind. Comp.
    %            eta:   the square of the xi values
    % -----------------------------------------------------------------------

    xi = w * Q;
    eta = xi .^ 2;



function [xid, etad] = projectionDerivatives(Q, w, xi)
    % -----------------------------------------------------------------------
    % [xid, etad] = projectionDerivatives(Z, w, xi)
    % 
    % Returns the partial derivatives of the projection variables xi and eta.
    % Inputs: 
    %              Q:    Q-matrix
    %              w:    considered basis projector
    %             xi:    direct projection of the samples over w
    %
    % Outputs:
    %            xid:    which stands for xi dot, a matrix where xid(i, b)
    %                    is the partial derivative of xi(i) respect to b
    %           etad:    which stands for eta dot, where etad(i, b) is the
    %                    partial derivative of eta(i) respect to b  
    % -----------------------------------------------------------------------

    N = length(w);

    xid = Q';
    etad = 2 * repmat(xi, [N 1])' .* xid;


function [SW, SB] = scatterMatrices(xi, eta, labels)
    % -----------------------------------------------------------------------
    % [SW, SB] = scatterMatrices(xi, eta, muP, muN, nuP, nuN)
    % 
    % Returns the within and between matrices SW, SB. Inputs: 
    %               xi:    direct projection of the samples over w
    %              eta:    the partial derivatives of eta
    %           labels     vector with the labels of the samples
    %
    % Output:
    %              SW:     the within matrix in a vector representation: 
    %                      SW = [SW(1,1), SW(2,2), SW(1,2)]
    %              SB:     the between matrix in a vector representation: 
    %                      SB = [SB(1,1), SB(2,2), SB(1,2)]
    % -----------------------------------------------------------------------

    % Preliminary Calculations
    xiN = xi(labels == -1);
    xiP = xi(labels == 1);
    etaN = eta(labels == -1);
    etaP = eta(labels == 1);

    N = length(xi);
    nN = length(xiN);
    pN = length(xiP);

    mu = mean(xi);
    nu = mean(eta);
    muN = mean(xiN);
    muP = mean(xiP);
    nuN = mean(etaN);
    nuP = mean(etaP);

    % Within Matrix
    nSW = 0;
    pSW = 0;

    for i = 1:nN
        v = [(xiN(i) - muN), (etaN(i) - nuN)];
        nSW = nSW + v' * v;
    end

    for i = 1:pN
        v = [(xiP(i) - muP), (etaP(i) - nuP)];
        pSW = pSW + v' * v;
    end

    W = 1/N * (nSW + pSW);

    SW = [W(1, 1), W(2, 2)];

    % Between Matrix
    nSB = nN * [muN - mu, nuN - nu]' * [muN - mu, nuN - nu];
    pSB = pN * [muP - mu, nuP - nu]' * [muP - mu, nuP - nu];

    B = 1/N * (nSB + pSB);

    SB = [B(1, 1), B(2, 2)];



function [SWd, SBd] = scatterDerivatives(xi, eta, labels, xid, etad)
    % -----------------------------------------------------------------------
    % [SWd, SBd] = catterDerivatives(xi, eta, labels, xid, etad)
    % 
    % Returns the partial derivatives of the within and between matrices 
    % Inputs: 
    %               xi:   direct projection of the samples over w
    %              eta:   the partial derivatives of eta
    %           labels:   vector with the labels of the samples
    %              xid:   partial derivatives of xi
    %             etad:   partial derivatives of eta
    % Outputs:
    %              SBd:   the partial derivatives of the three elements of 
    %                     the withing matrix, were SWd(i, beta) is the 
    %                     partial derivative respect to beta of the ith 
    %                     element of the vector [SB(1,1), SB(2,2), SB(1,2)]
    %              SWd:   the partial derivatives of the three elements of 
    %                     the withing matrix, were SWd(i, beta) is the 
    %                     partial derivative respect to beta of the ith 
    %                     element of the vector [SW(1,1), SW(2,2), SW(1,2)]
    % -----------------------------------------------------------------------
    
    % Preliminaries
    xiN = xi(labels == -1);
    xiP = xi(labels == 1);
    etaN = eta(labels == -1);
    etaP = eta(labels == 1);

    xiNd = xid(labels == -1, :);
    xiPd = xid(labels == 1, :);
    etaNd = etad(labels == -1, :);
    etaPd = etad(labels == 1, :);

    N = length(xi);
    nN = length(xiN);
    pN = length(xiP);

    mu = mean(xi);
    muN = mean(xiN);
    muP = mean(xiP);
    nu = mean(eta);
    nuN = mean(etaN);
    nuP = mean(etaP);

    mud = rowmean(xid);
    muNd = rowmean(xiNd);
    muPd = rowmean(xiPd);
    nud = rowmean(etad);
    nuNd = rowmean(etaNd);
    nuPd = rowmean(etaPd);

    % Within Matrix...
    pSW11d = 0;     
    pSW22d = 0;     

    for i = 1:pN
        pSW11d = pSW11d + 2 * (xiP(i) - muP) * (xiPd(i, :) - muPd);
        pSW22d = pSW22d + 2 * (etaP(i) - nuP) * (etaPd(i, :) - nuPd);
    end

    nSW11d = 0;
    nSW22d = 0;

    for i = 1:nN
        nSW11d = nSW11d + 2 * (xiN(i) - muN) * (xiNd(i, :) - muNd);
        nSW22d = nSW22d + 2 * (etaN(i) - nuN) * (etaNd(i, :) - nuNd);
    end


    SW11d = 1/N * (nSW11d + pSW11d);
    SW22d = 1/N * (nSW22d + pSW22d);


    SWd = [SW11d', SW22d'];


    % Between Matrix...
    pSB11d = 2 * pN * (muP - mu) * (muPd - mud);
    nSB11d = 2 * nN * (muN - mu) * (muNd - mud);

    pSB22d = 2 * pN * (nuP - nu) * (nuPd - nud);
    nSB22d = 2 * nN * (nuN - nu) * (nuNd - nud);

    SB11d = 1/N * (nSB11d + pSB11d);
    SB22d = 1/N * (nSB22d + pSB22d);

    SBd = [SB11d', SB22d'];



function m = rowmean(X)
    % -----------------------------------------------------------------------
    % m = rowmean(X)
    % 
    % Returns the mean of the rows of X
    % -----------------------------------------------------------------------
    
    if size(X, 1) == 1
        m = X;
    else
        m = mean(X);
    end

