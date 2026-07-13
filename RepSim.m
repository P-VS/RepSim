function RepSim

maskfile = '/Volumes/LaCie/UZ_Brussel/ASLBOLD_Manon/Groepsanalyse/Test31PP_PREcog Ses-001/DUNE-BOLD_HRF/mask_bold.nii';
outdir = '/Volumes/LaCie/UZ_Brussel/ASLBOLD_Manon/Groepsanalyse/Test31PP_PREcog Ses-001/';

outname = 'repSim_FALFF_Go-NoGo_p0.005';

ind_type = 'FALFF'; %Type of underlaying individual maps 'FALFF' or 'spmT'
mean_falff = 0.42;
sd_falff = 0.2;

pthr_ind = 0.55; %p threshold on the individual result maps
wwidth = 0.001; %weighting width in the weighting function: '1 if p=<pthr_ind'; 'exp(-(1/2)*((p-pthr_ind)/wwidth)^2) if p>pthr_ind'
pthr_group = 0.005; %p threshold on the individual result maps

nsub = 31; %number of simulated subjects
niter = 1000; %number of simulated experiiments

fwhm = [4 4 4]; %rp_Smoothest_gui; %smoothnes in the image. Default [4,4,4] or estimated from a statisticall map with rp_Smoothest_gui

repSim_simulation(maskfile,fwhm,ind_type,pthr_ind,wwidth,mean_falff,sd_falff,pthr_group,nsub,niter,outdir,outname);
