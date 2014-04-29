SDICA - Spatial Discriminant ICA (MATLAB implementation)
========================================================

Spatial Discriminant-ICA MATLAB implementation, with some scripts for 
importing fMRI data using SPM.  This version: April2014

Libs written by Alejandro Tabas <alextabas-AT-gmail-DOT-com>. Write me if you 
want to report a bug, ask for help or have a beer! ;)

This implementation is part of my Master's Thesis: "Basis Decomposition 
Discriminant ICA". If you are curious, you can download it from here:
https://upcommons.upc.edu/pfc/bitstream/2099.1/20739/1/Tabas.pdf

Much shorter and readable (but also less detailed and less geekie) is its
corresponding conference publication at the PRNI2014 (still in press). Please,
cite this last publication if you are using our method/libs for your paper :)

------------------------------------------------------------------------------

These libs work under MATLAB in Linux systems. SDICA libs should work softly 
and nicely also with windows and OCTAVE. The examples should require some
adaptation of the paths to work under windows (you know, the classical 
slash / and \ thingies) . The importation scripts could requiere some extra 
adaptation to work under octave because of the external libraries I use, 
namely: pcamat, for which I have included an Octave version provided by the 
authors (I haven't test it tho) and the SPM thingies for importing NIFTI 
files. Anyway, let me know if you need some help to adapt it, we can include
the new version in the repository afterwards :).

------------------------------------------------------------------------------

Main libs for the algorithm are in the SDICA folder. The main analysis is 
launch through SDICA.m. This function performs a characterisation-based SDICA
and return an object containing the relevant results. 
