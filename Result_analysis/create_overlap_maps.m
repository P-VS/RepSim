function [out]=create_overlap_maps(job)

for i=1:numel(job.indTmaps)
    V=spm_vol(job.indTmaps{i});
    mat=V.mat;
    
    Tim=spm_read_vols(V);
    
    pim=zeros(size(Tim));
    
    mask=find(abs(Tim)>0);
    
    [tpth,tnam,text] = spm_fileparts(job.indTmaps{i});
    
    spmmat=fullfile(tpth,'SPM.mat');
    
    if exist(spmmat)>0
        load(spmmat);
        if numel(mask)>0
            pim(mask)=(1-spm_Tcdf(Tim(mask),SPM.xX.erdf));
        else
            pim(:)=1;
        end
        
        stnam=split(tnam,'_');
    else
        stnam=split(tnam,'T');
        
        ipf=fullfile(tpth,['lP' stnam{2} '.nii']);
        
        VIP=spm_vol(ipf);
        ipim=spm_read_vols(VIP);
        
        ipmask=find(abs(ipim)>0);
        
        pim(ipmask)=10.^(-ipim(ipmask));
    end
    
    PM.dim=size(pim);
    PM.mat=mat;
    PM.fname=fullfile(tpth,['spmP_' stnam{2} '.nii']);
    PM=spm_create_vol(PM);
    PM=spm_write_vol(PM,pim);

    if job.omtype==1
        if job.thresok==1
            Tmask=find((pim<=job.pthres) & (pim>0));
        else
            Tmask=find(abs(Tim)>0);
        end
    
        if i==1
            dim=size(Tim);
            overlapmap=zeros(dim);
        
            overlapmap(Tmask)=1;
        else
            overlapmap(Tmask)=overlapmap(Tmask)+1;
        end     
    else
        Tmask=find(abs(Tim)>0);
        
        tempmap=zeros(size(Tim));
        
        p1mask=find(pim(Tmask)<=job.pthres);
        
        tempmap(Tmask(p1mask))=1;
        
        p2mask=find(pim(Tmask)>job.pthres);
        
        tempmap(Tmask(p2mask))=exp(-1/2*((pim(Tmask(p2mask))-job.pthres)/job.wwidth).^2);
        
        if i==1
            overlapmap=tempmap;
        else
            overlapmap=overlapmap+tempmap;
        end
    end
end

if job.percout==1
    overlapmap=overlapmap/numel(job.indTmaps);
    
    OLThres=job.olthres;
else
    OLThres=job.olthres*numel(job.indTmaps);
end

overlapmask = overlapmap;
THmask=find(overlapmap<OLThres);
overlapmask(THmask)=0;
overlapmask(overlapmap>=OLThres)=1;

[tpth,tnam,text] = spm_fileparts(job.indTmaps{1});

stnam=split(tnam,'_');

if size(stnam)<2, stnam{2}=''; end

thresname=strrep(num2str(OLThres),'.','');

if job.omtype==1
    ofname=fullfile(job.outdir{:},['ACM_spmT_' stnam{2} '_thres' thresname '.nii']);
    mfname=fullfile(job.outdir{:},['mask_ACM_spmT_' stnam{2} '_thres' thresname '.nii']);
else
    wname=strrep(num2str(job.wwidth),'.','');
    ofname=fullfile(job.outdir{:},['WOM_spmT_' stnam{2} '_thres' thresname '_w' wname '.nii']);
    mfname=fullfile(job.outdir{:},['mask_WOM_spmT_' stnam{2} '_thres' thresname '_w' wname '.nii']);
end

OM.dim=size(overlapmap);
OM.mat=mat;
OM.fname=ofname;
OM=spm_create_vol(OM);
OM=spm_write_vol(OM,overlapmap);

MM.dim=size(overlapmask);
MM.mat=mat;
MM.fname=mfname;
MM=spm_create_vol(MM);
MM=spm_write_vol(MM,overlapmask);

out=[OM.fname];