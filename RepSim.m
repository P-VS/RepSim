function RepSim

maskfile = '/Volumes/LaCie/UZ_Brussel/ME_fMRI_GE/data/sub-01/ses-001/SPMMAT-SE-EmoFaces_default/mask.nii';
outdir = '/Volumes/LaCie/UZ_Brussel/ME_fMRI_GE/Groeps_Analyses/EmoFaces/';

outname = 'repSim_se-fmri_default';

pthr = 0.001; %p threshold of the ressults map
wwidth = 0.001; %weighting width in the weighting function: '1 if p=<pthr'; 'exp(-(1/2)*((p-pthr)/wwidth)^2) if p>pthr'

nsub = 10; %number of simulated subjects
niter = 1000; %number of simulated experiiments

fwhm = rp_Smoothest_gui; %smoothnes in the image. Default [4,4,4] or estimated from a statisticall map with rp_Smoothest_gui

repSim_simulation(maskfile,fwhm,pthr,wwidth,nsub,niter,outdir,outname);
