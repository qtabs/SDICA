function x = loadfMRI(filepath)
    % -----------------------------------------------------------------------
    % x = loadfMRI(filepath)
    %
    % Hey! Don't forget I need SPM to work! :)
    %  
    % Loads an fMRI volume (using spm_vol) in a vector representation
    % Input:
    %            filepath:    a string indicating the location of the volume
    %                         (NIFTI files)
    % Output:
    %                  x:     a vector of size [1, #ofVoxels]. If volumes 
    %                         are also defined along time (i.e. 4D volumes)
    %                         x is actually an array where x{i} is the ith
    %                         3D volume of the time series
    % -----------------------------------------------------------------------

    headers = spm_vol(filepath);
    vols = spm_read_vols(headers);

    for i = 1:size(vols, 4)
        x{i} = reshape(vols(:, :, :, i), 1, numel(vols(:, :, :, i)));
    end