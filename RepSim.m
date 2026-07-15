function RepSim

maskfile = '/Volumes/LaCie/UZ_Brussel/ASLBOLD_Manon/Groepsanalyse/Test31PP_PREcog Ses-001/DUNE-BOLD_HRF/mask_bold.nii';
outdir = '/Volumes/LaCie/UZ_Brussel/ASLBOLD_Manon/Groepsanalyse/Test31PP_PREcog Ses-001/';

outname = 'repSim_ClusSize_Go-NoGo_p0.001';

ind_type = 'FALFF'; %Type of underlaying individual maps 'spmT'
thr_type = 'clustersize'; %Type of MC correction 'overlap' or 'clustersize'

mean_falff = 0.42;
sd_falff = 0.2;

pthr_ind = 0.001; %p threshold on the individual result maps
wwidth = 0.001; %weighting width in the weighting function: '1 if p=<pthr_ind'; 'exp(-(1/2)*((p-pthr_ind)/wwidth)^2) if p>pthr_ind'
pthr_group = 0.001; %p threshold on the individual result maps

nsub = 31; %number of simulated subjects
niter = 1000; %number of simulated experiiments

fwhm = [6 6 6];%rp_Smoothest_gui; %smoothnes in the image. Default [4,4,4] or estimated from a statisticall map with rp_Smoothest_gui

spm_path = '/Users/petervanschuerbeek/Library/Mobile Documents/com~apple~CloudDocs/Matlab/spm25';
addpath(genpath(spm_path));

repSim_simulation(maskfile,fwhm,ind_type,pthr_ind,thr_type,wwidth,mean_falff,sd_falff,pthr_group,nsub,niter,outdir,outname);
