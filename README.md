The RepSim toolbox is developed to perform Monte Carlo simulations of fMRI experiments. The resulting text file gives an idea of how likely it is, to have a false positive result repeated in n subjects out of a sample of N subjects.

The RepSim toolbox is created with the app designer in Matlab R2020b. The code to do the simulations is based on the AlphaSim tool as implemented in RESTplus v1.24 (http://restfmri.net/forum/). In order to work, SPM12 (https://www.fil.ion.ucl.ac.uk/spm/) should be installed.

The script performs a S simulatiosn of statistical maps for N subjects. For each stimulation an activation count map is created and a group analysis (t-test) is performed.
 
The activation count maps can be created using an exact threshold (no weighted threshold) or a weighted threshold. In case of the weighted threshold, the activation count is determined as
- 1 if p<=pthres
- exp(-(1/2)*((p-pthres)/width)^2) if p>pthres

The results are given in a tekst file:
- Number of subjects (n): a significant result is found in at least n out of N subjects
- Percentage of subjects: n given as percentage of N
- Frequency >=n: Number of times (in voxels) an overlap of a false positive result is seen in at least n subjects (summed over all simulations)
- Probability of >=n: Probability of having an overlap of a false positive result in at least n subjects
- Freq significant group: Number of voxels showing an overlap of a false positive result is seen in at least n subjects and a false positive group effect
- Probability significant group: Probability of having an overlap of a false positive result in at least n subjects and a false positive group effect
- Found in x simulations: Number of simulation in which an overlap of a false positive result in at least n subjects and a false positive group effect is seen in at least 1 voxel
- Mean freq per simulation: Mean number of voxels per simulation showing an overlap of a false positive result in at least n subjects and a false positive group effect
- SD freq per simulation: Standard deviation of the number of voxels per simulation showing an overlap of a false positive result in at least n subjects and a false positive group effect
- Max freq per simulation: Maximum number of voxels per simulation showing an overlap of a false positive result in at least n subjects and a false positive group effect


The 'Result_analysis' folder contains a SPM toolbox to create activation count maps (ACM) or weighted overlap maps (WOM) from real fMRI results. The whole 'Result_analysis' folder should be saved as a subfolder in the 'SPM12/toolbox' folder.