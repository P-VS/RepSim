function Result_analysis = tbx_cfg_Result_analysis

if ~isdeployed, addpath(fullfile(spm('dir'),'toolbox','Result_analysis')); end

% Individual result overlap
%--------------------------------------------------------
indTmaps         = cfg_files;
indTmaps.tag     = 'indTmaps';
indTmaps.name    = 'Individual T-maps';
indTmaps.help    = {'Select the individual thresholded T-maps'};
indTmaps.filter  = 'image';
indTmaps.ufilter = 'spmT.*';
indTmaps.num     = [1 Inf];

outdir         = cfg_files;
outdir.tag     = 'outdir';
outdir.name    = 'Output directory';
outdir.help    = {'Select the directory to save the overlap map in.'};
outdir.filter  = 'dir';
outdir.ufilter = '.*';
outdir.num     = [1 1];

omtype         = cfg_menu;
omtype.tag     = 'omtype';
omtype.name    = 'Type overlap map?';
omtype.help    = {['Selec the type of overlap map: '...
                    'ACM: Activation Count Map'...
                    'WOM: Weighted Overlap Map']};
omtype.values  = {1 2};
omtype.val     = {1};
omtype.labels  = {'ACM' 'WOM'};

thresok         = cfg_menu;
thresok.tag     = 'thresok';
thresok.name    = 'Apply p(unc) threshold?';
thresok.help    = {['Should a p(unc) threshold be applied?'...
                    'Always done for WOM!']};
thresok.values  = {1 2};
thresok.val     = {1};
thresok.labels  = {'Yes' 'No'};

pthres         = cfg_entry;
pthres.tag     = 'pthres';
pthres.name    = 'p(unc) threshold';
pthres.help    = {'Give the p(unc) threshold'};
pthres.val     = {0.001};
pthres.strtype = 'r';
pthres.num     = [1 1];

wwidth         = cfg_entry;
wwidth.tag     = 'wwidth';
wwidth.name    = 'Weighting width?';
wwidth.help    = {['Give the weighting width (w).'...
                   'Weightinig function: '...
                   '1 if p<=pthres'...
                   'exp(-(1/2)*((p-pthres)/w)^2) if p<pthres']};
wwidth.strtype = 'r';
wwidth.val     = {0.001};
wwidth.num     = [1 1];

percout         = cfg_menu;
percout.tag     = 'percout';
percout.name    = 'Output as percentage?';
percout.help    = {'Should the output be expressed as percentage or absolute number.'};
percout.values  = {1 2};
percout.val     = {1};
percout.labels  = {'Yes' 'No'};

olthres         = cfg_entry;
olthres.tag     = 'olthres';
olthres.name    = 'Threshold overlap in percentage?';
olthres.help    = {'Give the threshold for the minimal overlap in percentage.'};
olthres.strtype = 'r';
olthres.val     = {0.0};
olthres.num     = [1 1];

% Effect size
%--------------------------------------------------------
spmmat         = cfg_files;
spmmat.tag     = 'spmmat';
spmmat.name    = 'Select SPM.mat';
spmmat.help    = {'Select the SPM.mat file. The effect size is calculated for all contrasts as spmT/sqrt(df) with df the residual degrees of freedom.'};
spmmat.filter  = 'mat';
spmmat.ufilter = 'SPM.mat';
spmmat.num     = [1 Inf];

% ROI extraction
%--------------------------------------------------------
ROImaps         = cfg_files;
ROImaps.tag     = 'ROImaps';
ROImaps.name    = 'Index ROIs image';
ROImaps.help    = {'Select the ROI index file'};
ROImaps.filter  = 'image';
ROImaps.ufilter = '.*';
ROImaps.num     = [1 1];

statmaps         = cfg_files;
statmaps.tag     = 'statmaps';
statmaps.name    = 'Statistical maps';
statmaps.help    = {'Select the statistical image files'};
statmaps.filter  = 'image';
statmaps.ufilter = '.*';
statmaps.num     = [1 Inf];

excel          = cfg_entry;
excel.tag      = 'excel';
excel.name     = 'Name result file';
excel.help     = {'Give the name for the result table file'};
excel.strtype  = 's';
excel.val      = {'ROIexctraction'};
excel.num      = [1 Inf];

selvox         = cfg_menu;
selvox.tag     = 'selvox';
selvox.name    = 'Which voxels should be selected?';
selvox.help    = {['All: All voxels within the ROI are selected. ' ...
                   'Only positive: only voxels with a poistive value are selected']};
selvox.values  = {1 2};
selvox.val     = {1};
selvox.labels  = {'All' 'Only positive'};

% Main structure
%--------------------------------------------------------

ROIanalysis      = cfg_exbranch;
ROIanalysis.tag  = 'ROIanalysis';
ROIanalysis.name = 'ROI extraction';
ROIanalysis.val  = {ROImaps statmaps outdir excel selvox};
ROIanalysis.help = {'This function extractes the mean and disppersion metrix within ROIs from statistical maps'};
ROIanalysis.prog = @(job)vout_Result_analysis('ROI',job);
ROIanalysis.vout = @(job)vout_Result_analysis('vout',job);

effectsize      = cfg_exbranch;
effectsize.tag  = 'effectsize';
effectsize.name = 'Effect size';
effectsize.val  = {spmmat};
effectsize.help = {'This function determines the effect size of a contrast'};
effectsize.prog = @(job)vout_Result_analysis('effectsize',job);
effectsize.vout = @(job)vout_Result_analysis('vout',job);

overlap      = cfg_exbranch;
overlap.tag  = 'overlap';
overlap.name = 'Activation overlap';
overlap.val  = {indTmaps outdir omtype thresok pthres wwidth percout olthres};
overlap.help = {'This function determines the overlap between individual fMRI results'};
overlap.prog = @(job)vout_Result_analysis('overlap',job);
overlap.vout = @(job)vout_Result_analysis('vout',job);

Result_analysis         = cfg_choice;
Result_analysis.tag     = 'ResultAnalysis';
Result_analysis.name    = 'Result analysis';
Result_analysis.help    = {'This toolbox is meant to analyse fMRI results'};
Result_analysis.values  = {overlap effectsize ROIanalysis};

function out=vout_Result_analysis(cmd,job)

switch lower(cmd)
    case 'overlap'
        [out.files]=create_overlap_maps(job);
    case 'effectsize'
        [out.files]=calc_effect_size(job);
    case 'roi'
        [out.files]=ra_ROI_extraction(job);
    case 'vout'
        out(1)           =cfg_dep;
        out(1).sname     =sprintf('Overlap maps');
        out(1).src_output=substruct('.','files');
        out(1).tgt_spec  =cfg_findspec({{'filter','image','strtype','e'}});
end