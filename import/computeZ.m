function [Z, l] = computeZ_svd(filename, N)
    % -----------------------------------------------------------------------
    % computeZ(filename, N)
    % 
    % Computes PCA for huge files
    % 
    % Inputs:
    %              filename:     path to the file with the matrix (the matrix
    %						     should be contained in the file with the name
    %						     "X". The dimensionality reduction will be 
    %							 performed over the rows i.e. the final number
    %							 of cols will be always the same as the 
    %						     initial)
    %					  N: 	 max number of desired components
    % Output: 
    %					  Z: 	 The PCA-ed matrix
    %					  l: 	 Vector with the eigenvalues
    % -----------------------------------------------------------------------

    % Harcoded parameter: number of blocks to deal with the file 
    % (adjust manually depending on the size of the dataset and the memory 
    % of the computer)
    k = 25;

	mfile = matfile(filename);
	[M, D] = size(mfile, 'X');

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

	fprintf('Removing the mean... \n')
	X = [];
	save('tmpX.mat', 'X', '-V7.3');
	tmpXfile = matfile('tmpX.mat', 'Writable', true);

	for i = 1:length(m)
		fprintf('--- Block %d of %d...', i, length(m));
		tic;
		Xi = mfile.X(m{i}, 1:D);
		Xi = Xi - repmat(mean(Xi')', [1, D]);
		tmpXfile.X(m{i}, 1:D) = Xi;
		t = toc; 
		fprintf(' time: %.1fmin', t/60);
		fprintf('\n');
	end	

	clear tmpXfile mfile
	mfile = matfile('tmpX.mat'); 
 
	fprintf('Computing X-square...\n')

	for i = 1:length(m)
		fprintf('--- Block %d of %d...', i, length(m));
		tic;
		Xi = mfile.X(m{i}, 1:D);
		for j = 1:length(m)
			X2(m{i}, m{j}) = Xi * mfile.X(m{j}, 1:D)'; 
		end
		t = toc; 
		fprintf(' time: %.1fmin', t/60);
		fprintf('\n');
	end
	
	fprintf('...done!\n')	

    save('X2.mat', 'X2','-V7.3');
	fprintf('Computing eigenvalues...')
	tic;
    [U, S, V] = svd(X2);
    U=U(:, 1:N);
    S=S(1:N, 1:N);
    V=V(:, 1:N);    
    
	L = sqrt(S);
	t = toc;
	fprintf(' time: %.1fmin', t/60);
    
    for i = 1:size(L)
		l(i) = L(i, i);
	end

	save('./tmp.mat')
	clear U X2 L S; 
	fprintf(' ...done!\n')

	fprintf('Projecting the data matrix over the principal components...\n')
	d = ceil(D/k);
	e = 0;
	i = 0;

	while e < D
		i = i + 1;
		f = min(e + d, D);
		m{i} = (e + 1):f;
		e = f;
	end

	for i = 1:length(m)
	 	fprintf('--- Block %d of %d...', i, length(m));
	 	tic
	 	Xi = mfile.X(1:M, m{i});
 		Z(1:N, m{i}) = V' * Xi; 
 		t = toc; 
		fprintf(' time: %.1fmin', t/60);
	 	fprintf('\n');
	end

	for i = 1:N 
		Z(i, :) = Z(i, :) / sqrt(var(Z(i, :)));
	end

	fprintf('...done! \n\n')
