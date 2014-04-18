function [negentropy, delta] = deltaICA(Z, w0)
    % -----------------------------------------------------------------------
    % [negentropy, delta] = deltaICA(Z, w)
    % 
    % Returns the delta for a gradien ascend technique corresponding to a 
    % negentropy-based ICA using: 
    %             Z:      sample matrix (rows samples, cols variables)
    %            w0:      considered basis projector from previous iteration
    %
    % Output:
    %    negentropy:       value for negentropy of provided w0
    %         delta:       delta vector for the iteration: w <- w0 + delta
    % -----------------------------------------------------------------------

    [N, dim] = size(Z);

    % Proyecting samples
    y = w0 * Z;

    % Computing the g(·) needed for the expectation values
    G = -exp(-y.^2 / 2);
    g = -y .* G;
    
    % Approximation multiplicative constant
    k = 24/(16 * sqrt(3) - 27);

    % Computing gamma
    gamma = sum(G) / dim + 1/sqrt(2);

    for i = 1:N
        delta(i) = 2 * k * (gamma / dim) * (g * Z(i, :)');
    end

    % Computing negentropy
    negentropy = k * (sum(G) / dim + 1/sqrt(2))^2;