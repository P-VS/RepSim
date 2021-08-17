%This code is largely copied from the AphaSim tool as implemeted in the REST
%toolbox

%Author: dr. Peter Van Schuerbeek (UZ Brussel - VUB)

function repSim_simulation(maskfile,fwhm,pthr,nsub,iter,outdir,outname)

[maskpath, maskname, masketc]=fileparts(maskfile);
[mask,voxdim,header]=rp_readfile(maskfile);
[nx,ny,nz] = size(mask);
mask = logical(mask);
nxyz = numel(find(mask));
dvoxel=max(voxdim);

outfilename=outname;
outfile=fullfile(outdir,[outfilename,'.txt']);

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
counttabel = zeros(1,nsub);

count_sim = zeros(1,nsub);
max_count = zeros(1,nsub);

for nt = 1:iter
    fprintf(['Performing iteration ' num2str(nt) ' of ' num2str(iter) '\n'])
    
    countim = zeros(nx,ny,nz);
    
    for ns = 1:nsub
        fim=randn(nx,ny,nz);
        
        if fwhm(1)*fwhm(2)*fwhm(3) ~= 0
            fim = gauss_filter(fwhm,fim,voxdim); 
        end
        
        fimca=reshape(fim,1,[]);
        suma=sum(fimca); 
        sumsq=sum(fimca.*fimca);
        mean=suma/count;
        sd = sqrt((sumsq - (suma * suma)/count) / (count-1));
        
        zthr=-sqrt(2) * erfcinv((1-pthr)*2);
        
        xthr=sd*zthr+mean;
        
        fim2=fim;
        fim2(fim<=xthr)=0;      
        fim2(fim>xthr)=1;
        fim=fim2;
        
        fim = fim.*mask;
        
        countim = countim+fim;
    end
    
    for ci = 1:nsub
        counttabel(ci) = counttabel(ci)+numel(find(countim==ci));
        
        max_count(ci) = max(max_count(ci),numel(find(countim==ci)));
        
        if numel(find(countim==ci))>0; count_sim(ci)=count_sim(ci)+1; end
    end
        
end

prob_table = counttabel/(nxyz*iter);

p_bonf = 0.05/nxyz;

cum_prob = zeros(1,nsub);
for i=1:nsub
    cum_prob(i)=sum(prob_table(i:nsub));
end

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
fprintf(fid,'Number of subjects = %d\n',nsub);
fprintf(fid,'Number of Monte Carlo simulations = %d\n',iter);

fprintf(fid,'\nBonferroni corrected p(unc) = %.2e\n',p_bonf);

fprintf(fid,'\nNumber of subjects\tFrequency\tProbability of =n\tProbability of >=n\tFound in x simulations\tIn max voxels per simulation');

for i=1:nsub
    fprintf(fid,'\n%d\t%d\t%.3e\t%.3e\t%d\t%d',i,counttabel(i),prob_table(i),cum_prob(i),count_sim(i),max_count(i));
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