function [out]=calc_effect_size(job)

[spth,snam,cext] = spm_fileparts(job.spmmat{1});

load(job.spmmat{1})

df = SPM.xX.erdf;
xCon = SPM.xCon;

ncontrasts = numel(xCon);

for ic=1:ncontrasts
    spmT = fullfile(spth,xCon(ic).Vspm.fname);
    
    split_tnam = split(xCon(ic).Vspm.fname,'_');
    outname = ['spmES_' split_tnam{2}];
    
    tim = spm_vol(spmT);
    timage = spm_read_vols(tim);
    
    esimage = timage./sqrt(df);
    
    ESM.dim = size(esimage);
    ESM.mat = tim.mat;
    ESM.fname = fullfile(spth,outname);
    ESM = spm_create_vol(ESM);
    ESM = spm_write_vol(ESM,esimage);
end

out = '';