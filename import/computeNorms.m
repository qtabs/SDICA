function Xnorm = computeNorms(filename)
    % -----------------------------------------------------------------------
    % Xnorm = computeNormOfX(filename)
    % 
    % Returns the norm of the imstances in X
    % Input:
    %         filename:   path to the file with the matrix
    % -----------------------------------------------------------------------

    fprintf('Computing norms..')
    
    mfile = matfile(filename);
    [M, D] = size(mfile, 'X');

    k = 7;
    d = ceil(M/k);
    e = 0;
    i = 0;

    while e < M
        i = i + 1;
        f = min(e + d, M);
        m{i} = (e + 1):f;
        e = f;
    end

    for i = 1:length(m)
        X = mfile.X(m{i}, 1:D);
        Xnorm(m{i}) = sqrt(sum((X').^2));
        fprintf('.');
    end

    fprintf('done! \n')