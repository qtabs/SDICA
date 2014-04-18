function Q = computeQ(Xfilename, Z)
    % -----------------------------------------------------------------------
    % Q = computeQ(Xfilename, Z)
    % 
    % Computes PCA for huge files that do not fit into the computers memory 
    % 
    % Inputs:
    %              filename:     path to the file with the matrix (the matrix
    %						     should be contained in the file with the name
    %						     "X". The dimensionality reduction will be 
    %							 performed over the rows i.e. the final number
    %							 of cols will be always the same as the 
    %						     initial)
    %					  Z: 	 The PCA-ed matrix
    %  Output:			  Q: 	 The Q-matrix
    % -----------------------------------------------------------------------

    % Harcoded parameter: number of blocks to deal with the file
    k = 7;

	mfile = matfile(Xfilename);
	[M, D] = size(mfile, 'X');
	N = size(Z, 1);

	fprintf('Normalising Z...\n')

	for i = 1:N
		Zn(i, :) = Z(i, 1:D) / norm(Z(i, 1:D));
	end

	clear Z;

	fprintf('Computing Q-matrix...\n')

	d = ceil(M/k);
	e = 0;
	i = 0;

	while e < M
		i = i + 1;
		f = min(e + d, M);
		m{i} = (e + 1):f;
		e = f;
	end

	clear e f i d;

	for i = 1:length(m)
		tic
		fprintf('--- Block %d of %d...', i, length(m));
		k = 1;
		for j = m{i}
			X(k, :) = mfile.X(j, 1:D) / norm(mfile.X(j, 1:D));
			k = k + 1;
		end
		Q(:, m{i}) =  Zn * X';
		clear X;
		t = toc; 
		fprintf(' time: %.1fm\n', t/60);
	end
	
	fprintf('...done!\n')	
