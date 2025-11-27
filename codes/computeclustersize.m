function [maxcl,clustersize,clusterind] = computeclustersize(pval,critthr)
% Written by Siddharth Talwar
% Last edited on 26-11-2025
% computes clustersizes

clustersize = [];
clusterind = [];
maxcl=0;
ind = find(pval<critthr); %find pval less than threshold

dif = diff(ind); %find difference to get consecutive points
dif_ind = [dif ind(1:end-1)];

ind2 = find(dif == 1); % get indices of consecutive points

if isempty(ind2) == 1
    return
else
    
    ind2_ind = [ind2 dif_ind(ind2,2)];
    
    clustersize =[];clusterind={};
    
    if isempty(ind2)==0
        
        seed = ind2(1);
        
        ii=1;
        
        while seed < ind2(end)+1 %while starting point is less than end
            
            z=1;sz=1;
            
            while z == 1 && seed < ind2(end)+1
                
                if dif(seed) > 1
                    z = 2;    % diff ~=1, so exit loop through z
                else
                    seed = seed + 1; %found a consective point
                    sz  = sz + 1; % size of cluster increases by 1
                end
            end
            
            clusterind{ii} = []; %new cluster size found
            if sz > 1
                clustersize(ii) = sz; % note cluster size
                clusterind{ii} = [clusterind{ii}...
                    ind(seed)-sz+1:ind(seed-1)+1]; %note cluster indices
                ii = ii+1;
            else
                seed = seed + 1;
            end
            
        end
        maxcl = max(clustersize);
    end
end
