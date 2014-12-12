function subjects = loadFromFile(filename)
    % -----------------------------------------------------------------------
    % subjects = loadFromFile(filename)
    % 
    % Returns the list of subject paths and IDs in an array of cells. 
    % Input:
    %         filename: path to the file with the subject paths and IDs
    %         information
    % -----------------------------------------------------------------------

fid = fopen(filename);
nsubjects = 0;

disp('Getting number of subjects...')
%--------------------------------------------------------------------------
% Compute the number of subject in the file auto_control_subjects.txt
while 1
    tline = fgets(fid);
    nsubjects = nsubjects+1;
    if ~ischar(tline), break, end
    disp(tline)
end
fclose(fid);

nsubjects = nsubjects-1;
fprintf('Number of subjects: %u\n',nsubjects);

%--------------------------------------------------------------------------
disp('Getting the subject names: path name and extension ...')
subjects=[];

fid = fopen(filename);
for i = 1:nsubjects
    P = fgets(fid);
    subjects{i} = P;
end
fclose(fid);

end
