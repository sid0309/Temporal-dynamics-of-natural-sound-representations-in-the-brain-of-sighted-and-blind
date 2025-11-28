function [src,pval] = fn_signperm_src(mat1,mat2,permnum,thr)

% 3D Sign permutation test for only source searchlight
% Written by Siddharth Talwar
% Last edited on 26-11-2025
% The script performs sign permutation testing and clusterwise correction.
% Inputs
% 1) mat1, mat2 = data should be source x subjects. Time should be
%                 averaged before.
% 2) permnum    = number of permutations
% 3) thr        = threshold (0.05)
%
% Output
% 1) src  = fieldtrip format src data structure that can be
%           plotted further. See fields pow and mask.
% 2) pval = uncorrected at every sourcepoint
%

src = [];
src.pow = 0;
if isempty(mat2) == 1
    meanvec  = squeeze(mean(mat1,2)); % mean
    nulldist = zeros(permnum,1);% null dist
    nd2pmaps = zeros(length(meanvec),size(nulldist,1)); %nulldist to pmaps
    
    tic
    for srcpt = 1:size(meanvec,1)
        num = randi([1,size(mat1,2)],1,permnum); % choose a number
        
        for perm = 1:permnum
            
            r       = randperm(size(mat1,2),num(perm)); % choose indices
            tmp     = mat1(srcpt,:); % vector of all observed values
            tmp(r)  = tmp(r)*-1; %flip sign the chosen indices
            nulldist(perm) = mean(tmp); % mean the result
            
        end
        
        % pval is fraction of values higher in n.d. than observed mean and
        % total number of permutations
        pval(srcpt) = length(find(nulldist >= meanvec(srcpt))) / permnum;
        
        % convert the null distribution to p-maps
        % find unique values, saves time removes redundancy
        [uniq, ~,ic] = unique(nulldist);
        
        tmp = zeros(size(uniq,1),1); %find pvalues of unique values in n.d.
        for ii = 1:size(uniq,1)
            tmp(ii) = sum((nulldist - uniq(ii))>0) / permnum;
        end
        nd2pmaps(srcpt,:) = tmp(ic); % get all p-values
    end
    
else
    meanvec  = squeeze(mean(mat1,2)) - squeeze(mean(mat2,2)); % mean
    nulldist = zeros(permnum,1);% null dist
    nd2pmaps = zeros(length(meanvec),size(nulldist,1)); %nulldist to pmaps
    
    tic
    for srcpt = 1:size(meanvec,1)
        num = randi([1,size(mat1,2)],1,permnum); % choose a number
        
        for perm = 1:permnum
            
            r1       = randperm(size(mat1,2),num(perm)); % choose indices
            tmp1      = mat1(srcpt,:); % vector of all observed values
            tmp1(r1)  = tmp1(r1)*-1; %flip sign the chosen indices
            
            r2       = randperm(size(mat2,2),num(perm)); % choose indices
            tmp2      = mat2(srcpt,:); % vector of all observed values
            tmp2(r2)  = tmp2(r2)*-1; %flip sign the chosen indices
            
            nulldist(perm) = mean(tmp1) - mean(tmp2); % mean the result
            
        end
        
        % pval is fraction of values higher in n.d. than observed mean and
        % total number of permutations
        pval(srcpt) = length(find(nulldist >= meanvec(srcpt))) / permnum;
        
        % convert the null distribution to p-maps
        % find unique values, saves time removes redundancy
        [uniq, ~,ic] = unique(nulldist);
        
        tmp = zeros(size(uniq,1),1); %find pvalues of unique values in n.d.
        for ii = 1:size(uniq,1)
            tmp(ii) = sum((nulldist - uniq(ii))>0) / permnum;
        end
        nd2pmaps(srcpt,:) = tmp(ic); % get all p-values
        

    end
    
end

%%
load('/media/siddharth/DATA/CPP/Projects/Aud_Cat/templategrid_gray_10mm.mat');
pos = template_grid.pos(template_grid.inside==1,:);

% find max cluster size of each nd
for ii = 1:permnum
    x1 = pos(find(nd2pmaps(:,ii)<0.05),:);
    T  = clusterdata(x1,0.1);
    maxcl(ii) = max(histc(T,unique(T)));
end
cluster_thr = prctile(maxcl,95); % cluster theshold

pos_emp   = pos(find(pval<thr),:); % find points which are significant
if isempty(pos_emp) == 0
    T_emp     = clusterdata(pos_emp,0.1); % cluster the significant points
    tmp1      = histc(T_emp,unique(T_emp)); %find size of all clusters
    tmp       = find(tmp1>cluster_thr); %find the clusters bigger than cluster threshold
    
    sigpos    = pos_emp(ismember(T_emp,tmp),:); %sigificant positions
    
    % create power and mask
    template_grid.pow = zeros(length(template_grid.inside),1);
    ind = find(template_grid.inside==1);
    
    jj=1;
    for ii = 1:length(ind)
        template_grid.pow(ind(ii)) = meanvec(jj,:);
        jj= jj +1;
    end
    
    template_grid.mask = zeros(length(template_grid.inside),1);
    template_grid.mask(find(ismember(template_grid.pos,sigpos,'rows'))) = 1;
    template_grid.pow(template_grid.mask==0) = 0;
      
%     src = template_grid;

else
    disp('No significant clusters found');
end