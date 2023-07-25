function [out]=calc_effect_size(job)

[spth,snam,cext] = spm_fileparts(job.spmmat{1});

load(job.spmmat{1})

df = SPM.xX.erdf;
xCon = SPM.xCon;
betas = SPM.Vbeta;
cbeta = fullfile(spth,[betas(numel(betas)).fname]);

CB = spm_vol(cbeta);
cbetaim = spm_read_vols(CB);
pcsmask = find(abs(cbetaim)>0);

xX = SPM.xX.X;

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
    ESM.descrip = 'Effect size (Cohen D-score)';
    ESM.pinfo = [1,0,0];
    ESM.dt = [spm_type('float32'),spm_platform('bigend')];
    ESM.n = [1 1];

    ESM = resanat_write_vol_4d(ESM,esimage);

    spmC = fullfile(spth,xCon(ic).Vcon.fname);
    cim = spm_vol(spmC);
    cimage = spm_read_vols(cim);

    pcsimage = zeros(size(cbetaim));
    pcsimage(pcsmask) = (cimage(pcsmask)*max(xX,[],'all')*100)./cbetaim(pcsmask);

    outname = ['spmPCS_' split_tnam{2}];

    PCS.dim = size(pcsimage);
    PCS.mat = cim.mat;
    PCS.fname = fullfile(spth,outname);
    PCS.descrip = 'Percent Signal Change (PCS)';
    PCS.pinfo = [1,0,0];
    PCS.dt = [spm_type('float32'),spm_platform('bigend')];
    PCS.n = [1 1];

    PCS = resanat_write_vol_4d(PCS,pcsimage);
end

out = '';