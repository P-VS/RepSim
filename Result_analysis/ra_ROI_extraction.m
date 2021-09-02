function [out]=ra_ROI_extraction(job)

fid=fopen(fullfile(job.outdir{:},[job.excel '.txt']),'w');

fprintf(fid,'%s\t','ROI index');
fprintf(fid,'%s\t','Stat file');
fprintf(fid,'%s\t\t','ROI size');

fprintf(fid,'%s\t','Mean');
fprintf(fid,'%s\t\t','SD');

fprintf(fid,'%s\t','Minimum');
fprintf(fid,'%s\t','Q1 (25% quantile)');
fprintf(fid,'%s\t','Q2 (Median)');
fprintf(fid,'%s\t','Q3 (75% quantile)');
fprintf(fid,'%s\n','Maximum');

VROI=spm_vol(job.ROImaps{1});
ROIim=spm_read_vols(VROI);

ROIim=round(ROIim);

amask=find(ROIim>0);

for nr=1:max(ROIim(:))
    mask=find(ROIim(amask)==nr);
    for ni=1:numel(job.statmaps)
        fprintf(fid,'%s\t',num2str(nr));
        fprintf(fid,'%s\t',job.statmaps{ni});
        
        VS=spm_vol(job.statmaps{ni});
        Sim=spm_read_vols(VS);

        if job.selvox==1
            Svals=Sim(amask(mask));
        else
            mask2 =find(Sim(amask(mask))>0);
            Svals=Sim(amask(mask(mask2)));
        end
        
        fprintf(fid,'%s\t\t',num2str(numel(Svals)));
        
        MN=mean(Svals(:));
        SD=std(Svals(:));
        
        fprintf(fid,'%8.2f\t',MN);
        fprintf(fid,'%8.2f\t\t',SD);
        
        minR=min(Svals(:));
        maxR=max(Svals(:));
        q1R=quantile(Svals(:),0.25);
        q2R=quantile(Svals(:),0.50);
        q3R=quantile(Svals(:),0.75);
        
        fprintf(fid,'%8.2f\t',minR);
        fprintf(fid,'%8.2f\t',q1R);
        fprintf(fid,'%8.2f\t',q2R);
        fprintf(fid,'%8.2f\t',q3R);
        fprintf(fid,'%8.2f\n',maxR);
        
    end
end

fclose(fid);
    
out=[''];