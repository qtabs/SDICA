addpath('./SDICA/');

% -- Data (matrices generated using the import libs) --
xfile = './data/X.mat';
qfile = './data/Q.mat';
zfile = './data/Z.mat';
normsfile = './data/Xnorm.mat';
resultsFold = './results/';
mfile = matfile(xfile);
labels = mfile.labels;
% Need to convert the problem into a two-classes problem?
% labels = 2 * (not(mfile.labels == 0) - 0.5);


% -- Parameters --
[M, dim] = size(mfile, 'X');
% Number of components to extract
n = [1, 2, 4, 6, 10, 15, 20];
% Number of elements in Z
N = [10, 20, 50, 100, 200, 500, 1000, 2000];
% K for the K-folds crossvalidation
K = 10;
% Value of the weighting coefficient kappa
Kappa = [0, .2, .4, .6, .8, 1];


% -- Experiments --
for kappa = Kappa

	for i = 1:length(N)

		load(qfile);
		load(zfile);
		load(normsfile);

		Z = Z(1:N(i), :);
		Q = Q(1:N(i), :);

		cv = cvpartition(M, 'Kfold', K);
		
		for l = 1:K
			fprintf('\n\nLearning projections for part. %d', l);
			fprintf(' (N = %d, k = %.1f)\n', N(i), kappa);
			Xh{i}{l} = feSDICA(xfile, labels, max(n), Z, Q, Xnorm, ...
									                   kappa, cv.training(l));
		end

		% Saves the new representation of the data for later use
		saveXhatfile = [resultsFold, './K', int2str(kappa * 10), 'Xhat.mat'];
		save(saveXhatfile, 'Xh', 'n', 'N');
	end

	clear Z Q


	for i = 1:length(N)

		for j = 1:length(n)

			fprintf('Testing for (N, n) = (%d, %d)\n', N(i), n(j));

			for l = 1:K

				X = Xh{i}{l}(:, 1:2*n(j));

				Xl = X(cv.training(l), :);
				Yl = labels(cv.training(l));
				Xt = X(cv.test(l), :);
				Yt = labels(cv.test(l));

				fprintf('-- Training SVM %d/%d', l, K);
				svm = svmtrain(Xl, Yl, 'kernel_function','linear', ...
							   'method','LS');
				YpredT = svmclassify(svm, Xt);
				perfL(l) = sum(YpredT == Yt') / length(Yt);
				YpredL = svmclassify(svm, Xl);
				perfT(l) = sum(YpredL == Yl') / length(Yl);
				fprintf('... done! Perf: (%.2f, %.2f)\n', perfT(l), perfL(l));
			
			end
			
			perfTest(i, j) = mean(perfL);
			devTest(i, j) = sqrt(var(perfL));		
			perfTrain(i, j) = mean(perfT);
			devTrain(i, j) = sqrt(var(perfT));

		end
	end

	clear Xh 

	perfFile = [resultsFold, 'K', int2str(kappa * 10), 'perf.mat'];
	save(perfFile, 'perfTest', 'devTest', 'perfTrain', 'devTrain', 'N', 'n');

	fprintf('\n\n --------------------- RESULTS --------------------- \n')
	fprintf(' |   N   |   n   |  Training Set  |  Testing Set   |\n');
	fprintf(' ---------------------------------------------------\n');

	for i = 1:length(N)
		for j = 1:length(n)
			fprintf(' | %4d  |  %2d   |', N(i), n(j));
			fprintf('  %.2f +- %.2f  |', perfTrain(i, j), devTrain(i, j));
			fprintf('  %.2f +- %.2f  |\n', perfTest(i, j), devTest(i, j));
		end
	end

	fprintf(' ---------------------------------------------------\n');

end

