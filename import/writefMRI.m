function writefMRI(vector, originalVolPath, targetPath)
    % -----------------------------------------------------------------------
    % writefMRI(vector, originalVolPath, targetPath)
    %  
    % Hey! Don't forget I need SPM to work! :)
    
    % Writes an fMRI volume represented as a vector in an specified path
    % Inputs:
    %               vector:   vector representation of the volume as given 
    %                         by loadfMRI
    %      originalVolPath:   The original volume imported by SPM, necessary
    %                         to reuse the header of the file
    %           targetPath:   Path to write the volume as a NIFTI file
    % -----------------------------------------------------------------------

    headers = spm_vol(originalVolPath);
    headers = headers(1);
    headers.fname = targetPath;
    spmVol = reshape(vector, headers.dim(1), headers.dim(2), headers.dim(3));
    spm_write_vol(headers, spmVol);
