function [pval,nulldist,sigind] = fn_signperm_chan(vec1,vec2,sample,permnum,thr,lay)

% 2D Sign permutation test for only sensor searchlight
% Written by Siddharth Talwar
% Last edited on 26-11-2025
% The script performs sign permutation testing and clusterwise correction.
% Not compatible for other 2D tests. See fn_signperm2D.
% Inputs
% 1) vec1, vec2 = data should be channels x subjects. Time should be
%                 averaged before.
% 2) sample     = which tails? 1 or 2
% 3) permnum    = number of permutations
% 4) thr        = threshold (0.05)
% 5) lay        = layout (see fieldtrip doc)
%
% Output
% 1) pval     = uncorrected pvalues at all timepoints
% 2) nulldist = null distribution
% 3) sigind   = significant timepoints after clusterwise correction

%%

% 1 sample
if sample == 1
    meanvec  = squeeze(mean(vec1,2)); % mean
    nulldist = zeros(permnum,1);% null dist
    nd2pmaps = zeros(length(meanvec),size(nulldist,1)); %nulldist to pmaps
    
    tic
    for chan = 1:size(meanvec,1)
        num = randi([1,size(vec1,2)],1,permnum); % choose a number
        
        for perm = 1:permnum
            
            r       = randperm(size(vec1,2),num(perm)); % choose indices
            tmp     = vec1(chan,:); % vector of all observed values
            tmp(r)  = tmp(r)*-1; %flip sign the chosen indices
            nulldist(perm) = mean(tmp); % mean the result
            
        end
        
        % pval is fraction of values higher in n.d. than observed mean and
        % total number of permutations
        pval(chan) = length(find(nulldist >= meanvec(chan))) / permnum;
        
        % convert the null distribution to p-maps
        % find unique values, saves time removes redundancy
        [uniq, ~,ic] = unique(nulldist);
        
        tmp = zeros(size(uniq,1),1); %find pvalues of unique values in n.d.
        for ii = 1:size(uniq,1)
            tmp(ii) = sum((nulldist - uniq(ii))>0) / permnum;
        end
        nd2pmaps(chan,:) = tmp(ic); % get all p-values
    end
    toc
elseif sample ==2
    
    meanvec  = squeeze(mean(vec1,2)) - squeeze(mean(vec2,2)); % mean
    nulldist = zeros(permnum,1);% null dist
    nd2pmaps = zeros(length(meanvec),size(nulldist,1)); %nulldist to pmaps
    
    tic
    for chan = 1:size(meanvec,1)
        num = randi([1,size(vec1,2)],1,permnum); % choose a number
        
        for perm = 1:permnum
            
            r1       = randperm(size(vec1,2),num(perm)); % choose indices
            tmp1     = vec1(chan,:); % vector of all observed values
            tmp1(r1) = tmp1(r1)*-1; %flip sign the chosen indices
            
            r2       = randperm(size(vec2,2),num(perm)); % choose indices
            tmp2     = vec2(chan,:); % vector of all observed values
            tmp2(r2) = tmp2(r2)*-1; %flip sign the chosen indices
            
            nulldist(perm) = mean(tmp1) - mean(tmp2); % mean the result
            
        end
        
        % pval is fraction of values higher in n.d. than observed mean and
        % total number of permutations
        pval(chan) = length(find(nulldist >= meanvec(chan))) / permnum;
        
        % convert the null distribution to p-maps
        % find unique values, saves time removes redundancy
        [uniq, ~,ic] = unique(nulldist);
        
        tmp = zeros(size(uniq,1),1); %find pvalues of unique values in n.d.
        for ii = 1:size(uniq,1)
            tmp(ii) = sum((nulldist - uniq(ii))>0) / permnum;
        end
        nd2pmaps(chan,:) = tmp(ic); % get all p-values
    end
    toc
    
end

%% correction
% Layout to image
Cartesian = (lay.pos + abs(min(min(lay.pos))) + 0.1) * 100;
imSize=120;
binaryImage= false(imSize);
Cartesian(:,2) = imSize-(Cartesian(:,2)-1); % Comment if you want zero at the top-left corner
for idx=1:size(Cartesian,1)
    binaryImage(floor(Cartesian(idx,2)),floor(Cartesian(idx,1))) = true; % Note that X-Y index are inverted in indexing
end

%% Max Clustersize in n.d.
tic

disp('calculating max cluster size in the null distribution');
maxcl_perm = zeros(1,perm);
for perm = 1:(permnum + 1) % for each perm of null dist (+1 is the empirical)
    fprintf(num2str(perm));
    if perm == permnum+1
        seedchan = find(pval<thr);
    else
        seedchan = find(nd2pmaps(:,perm)<thr); % find seed channels
        seedchan = seedchan';
    end
    
    if isempty(seedchan)==0
        null_clustersize_chan = zeros(1,128);
        
        maxcl_chan = zeros(1,128);
        neigh=cell(1,128);
        x={};
        
        for idx = seedchan % finding maxcl taking each channel as seed
            
            previous_maxclsize = -1;
            d=[];
            searchsphere = 1;
            chancount = -1;
            
            curr_idx = [floor(Cartesian(idx,2)),floor(Cartesian(idx,1))]; % find the channel in image
            
            for idx2 = 1:128 % find distance of seed channel with others
                d(idx2) = pdist2(curr_idx,...
                    [floor(Cartesian(idx2,2)),floor(Cartesian(idx2,1))]);
            end
            
             % sort distance of all channels from seed channel
            [d,srt] = sort(d);
            
            for search = 1:1000 % increasing search sphere by iteration
                
                % search for channels in sphere
                neigh{idx} = srt(find(d<=searchsphere));
                
                % if the new bigger sphere size includes same num of
                % channels, break
                if chancount == numel(neigh{idx})
                    searchsphere = searchsphere+0.01;
                    
                    continue
                else
                    chancount = numel(neigh{idx});
                    
                end
                
                % find if the closest channels are also significant
                for ii = 1:length(neigh{idx})
                    x{idx}(ii) = ismember(neigh{idx}(ii),seedchan);
                end
                
                % update max cluster size of seed channel
                if any(x{idx})==1
                    null_clustersize_chan(idx) = numel(find(x{idx}==1));
                end
                
                % break loop if the null_clustersize_chan does not update in successive search
                % iterations, else keep going on (cluster is bigger)
                if null_clustersize_chan(idx) == previous_maxclsize
                    % found the biggest cluster, ciao
                    maxcl_chan(idx) = null_clustersize_chan(idx);
                    break
                    
                else
                    % make the sphere bigger and search again
                    previous_maxclsize = null_clustersize_chan(idx);
                    searchsphere = searchsphere+0.01;
                
                end
            end 
        end
        
        maxcl_perm(perm) = max(maxcl_chan);
    end
    fprintf(repmat('\b', 1, length(mat2str(perm))));
end
maxclsize = prctile(maxcl_perm(1:permnum),95);
toc

%% find if empirical clusters are bigger than maxlsize

if maxclsize<=1
    maxclsize = 1; % min neighbours reqd
end
% maxclsize = 1;
sig_seedchan = find(null_clustersize_chan>=maxclsize);

if isempty(sig_seedchan) == 0
    disp('finding significant clusters');
    sigind = [];
    for sigcl = 1:length(sig_seedchan)
        sigind = [sigind sig_seedchan neigh{sig_seedchan(sigcl)}(x{sig_seedchan(sigcl)}==1)];
    end
    sigind = unique(sigind);
    
else
    sigind = [];
    disp('no significant clusters found');
end

