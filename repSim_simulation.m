%This code is largely copied from the AphaSim tool as implemeted in the REST
%toolbox

%Author: dr. Peter Van Schuerbeek (UZ Brussel - VUB)

function repSim_simulation(maskfile,fwhm,pthr,wwidth,nsub,iter,outdir,outname)

[maskpath, maskname, masketc]=fileparts(maskfile);
[mask,voxdim,header]=rp_readfile(maskfile);
[nx,ny,nz] = size(mask);
mask = logical(mask);
nxyz = numel(find(mask));
dvoxel=max(voxdim);

V=spm_vol(maskfile);

outfile=fullfile(outdir,[outname,'.txt']);

if fwhm(1) == 4
    fwhm(1) = 4.55;  %afni 4 : spm 4.55 
end
if fwhm(2) == 4
    fwhm(2) = 4.55;  %afni 4 : spm 4.55 
end
if fwhm(3) == 4
    fwhm(3) = 4.55;  %afni 4 : spm 4.55 
end

count=nx*ny*nz;
counttabel = zeros(1,nsub+1);

count_sim = zeros(1,nsub+1);
max_count = zeros(1,nsub+1);
iter_grcount = zeros(iter,nsub+1);
iter_count = zeros(iter,nsub+1);

gr_count = zeros(1,nsub+1);

for nt = 1:iter
    fprintf(['Performing iteration ' num2str(nt) ' of ' num2str(iter) '\n'])
    
    countim = zeros(nx,ny,nz);
    
    fim4d = zeros(nx,ny,nz,nsub);
    
    for ns = 1:nsub
        fim=randn(nx,ny,nz);
        
        if fwhm(1)*fwhm(2)*fwhm(3) ~= 0
            fim = gauss_filter(fwhm,fim,voxdim); 
        end
        
        fim4d(:,:,:,ns)=fim;
        
        fim = 1-normcdf(fim,mean(fim,'all'),std(fim,0,'all'));
        
        fim2=fim;
        if wwidth==0
            fim(fim2<=pthr)=1;      
            fim(fim2>pthr)=0;
        else
            fim(fim2<=pthr)=1;      
            fim(fim2>pthr)=exp(-(1/2)*((fim2(fim2>pthr)-pthr)/wwidth).^2);
        end
        
        fim = fim.*mask;
        
        countim = countim+fim;
    end
    
    meangr = mean(fim4d,4);
    stdgr = std(fim4d,0,4);
    ttestgr = meangr./(stdgr/sqrt(nsub));
    ttestgr(~isfinite(ttestgr))=0;
    
    ttestgr = 1-tcdf(ttestgr,nsub-1,'upper');

    ttestgr2=ttestgr;
    ttestgr(ttestgr2<=pthr)=1;
    ttestgr(ttestgr2>pthr)=0;

    ttestgr = ttestgr.*mask;
    
    for ci = 0:nsub
        
        indmask = find(and(countim(mask)>=ci-0.5,countim(mask)<ci+0.5));
        
        counttabel(ci+1) = counttabel(ci+1)+numel(indmask);
        iter_count(nt,ci+1) = numel(indmask);

        grmask=find(and(and(countim>=ci-0.5,countim<ci+0.5),ttestgr>0));
        
        gr_count(ci+1) = gr_count(ci+1)+numel(grmask);

        max_count(ci+1) = max(max_count(ci+1),numel(grmask));
        iter_grcount(nt,ci+1) = numel(grmask);

        if numel(grmask)>0; count_sim(ci+1)=count_sim(ci+1)+1; end  
    end     
end

mean_count = round(mean(iter_grcount,1))';
sd_count = round(std(iter_grcount,0,1))';

prob_table = counttabel/(nxyz*iter);
gr_prob = gr_count/(nxyz*iter);

repsimulation.resultfile = outfile;
repsimulation.maskfile = maskfile;
repsimulation.fwhm = fwhm;
repsimulation.pthr = pthr;
repsimulation.wwidth = wwidth;
repsimulation.nsub = nsub;
repsimulation.iter = iter;
repsimulation.nxyz = nxyz;

repsimulation.counttabel = counttabel;
repsimulation.prob_table = prob_table;
repsimulation.gr_count = gr_count;
repsimulation.gr_prob = gr_prob;
repsimulation.count_sim = count_sim;
repsimulation.mean_count = mean_count;
repsimulation.sd_count = sd_count;
repsimulation.max_count = max_count;

repsimulation.iter_count = iter_count;
repsimulation.iter_grcount = iter_grcount;

matfile = fullfile(outdir,[outname,'.mat']);

save(matfile,'repsimulation');

fid=fopen(sprintf('%s',outfile),'w');

fprintf(fid,'RepSim: Monte Carlo simulations to determine the probability of repeated false positive results');
fprintf(fid,'\nThis tool is bassed on AlphaSim as implemented in REST');
fprintf(fid,'\nAuthor: dr. Peter Van Schuerbeek (UZ Brussel - VUB)\n');

fprintf(fid,'\nMask filename = %s\n',maskfile);
fprintf(fid,'Voxels in mask = %d\n',nxyz);
fprintf(fid,'Gaussian filter width (FWHMx, in mm) = %.3f\n',fwhm(1));
fprintf(fid,'Gaussian filter width (FWHMy, in mm) = %.3f\n',fwhm(2));
fprintf(fid,'Gaussian filter width (FWHMz, in mm) = %.3f\n',fwhm(3));
fprintf(fid,'Individual voxel threshold probability = %.3f\n',pthr);
if wwidth>0
    fprintf(fid,'Weidthed filter width = %.3f\n',wwidth);
    fprintf(fid,'Applied weighting function:\n');
    fprintf(fid,'1 if p<=pthres\n');
    fprintf(fid,'exp(-(1/2)*((p-pthres)/width)^2) if p>pthres\n');
else
    fprintf(fid,'No wighting filter applied\n');
end
    
fprintf(fid,'Number of subjects = %d\n',nsub);
fprintf(fid,'Number of Monte Carlo simulations = %d\n',iter);

fprintf(fid,['\nNumber of subjects (n)', ...
             '\tPercentage of subjects',...
             '\tFrequency =n',...
             '\tProbability of =n', ...
             '\tFreq significant group',...
             '\tProbability significant group', ...
             '\tFound in x simulations',...
             '\tMean freq per simulation',...
             '\tSD freq per simulation',...
             '\tMax freq per simulation',...
             '\tp corrected']);

for i=0:nsub
    p_corrected = 1-(1-gr_prob(i+1))^nxyz;
    fprintf(fid,'\n\t%d\t%d\t%d\t%.3e\t%d\t%.3e\t%d\t%d\t%d\t%d\t%.3e',i,(i/nsub)*100,counttabel(i+1),prob_table(i+1),gr_count(i+1),gr_prob(i+1),count_sim(i+1),mean_count(i+1),sd_count(i+1),max_count(i+1),p_corrected);
end

fclose(fid);

fprintf('Done\n')

end


function Q=gauss_filter(s,P,VOX)
if length(s) == 1; s = [s s s];end 
s  = s./VOX;					% voxel anisotropy
s  = max(s,ones(size(s)));			% lower bound on FWHM
s  = s/sqrt(8*log(2));				% FWHM -> Gaussian parameter

x  = round(6*s(1)); x = [-x:x];
y  = round(6*s(2)); y = [-y:y];
z  = round(6*s(3)); z = [-z:z];
x  = exp(-(x).^2/(2*(s(1)).^2)); 
y  = exp(-(y).^2/(2*(s(2)).^2)); 
z  = exp(-(z).^2/(2*(s(3)).^2));
x  = x/sum(x);
y  = y/sum(y);
z  = z/sum(z);

i  = (length(x) - 1)/2;
j  = (length(y) - 1)/2;
k  = (length(z) - 1)/2;
Q=P;
spm_conv_vol(P,Q,x,y,z,-[i,j,k]);
end