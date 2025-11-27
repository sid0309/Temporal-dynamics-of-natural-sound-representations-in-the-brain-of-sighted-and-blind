function [pseudo_tr,group,n_tr] = make_pseudobycat(preproc,ps_k,ps_n,ncat)
%%
if (isempty(ps_k)==1 && isempty(ps_n)==1) || (isempty(ps_k)==0 && isempty(ps_n)==0)
    error('enter number of trials to avg OR number of pseudotrials reqd.');
end

% epoch by condition
trinfo = preproc.trialinfo;

for trig = 1:24
    epochs{trig} = zeros(size(preproc.trial{1},1),size(preproc.trial{1},2),1); % create data structure
end

for tr = 1:length(preproc.trial)
    epochs{trinfo(tr)}(:,:,size(epochs{trinfo(tr)},3)+1) = preproc.trial{tr}; % ep_ds_bl_tr{trigno}(ch x tp x trials)
end

if ncat == 2
    
    sounds_in_cat{1} = [1:6 13:18]; % sounds arranged like this
    sounds_in_cat{2} = [7:12 19:24];
    
elseif ismember(ncat,[4 8 24]) == 1
    
    for ii = 1:ncat
        sounds_in_cat{ii} = (1:length(epochs) /...
            ncat)...
            + ((length(epochs) /ncat)*(ii-1));
    end
    
else
    error('param.n_category can be 2,4,8 or 24');
end

dataset = {}; 
for ii = 1 : ncat
    dataset{ii} = cat(3,epochs{sounds_in_cat{ii}});
    
end
epochs = [];
epochs = dataset;

for trig = 1:length(epochs)
    epochs{trig}(:,:,1) = []; % remove the 1st empty dummy trial (because of initializing and +1 in the above step)
    mintr(trig)         = size(epochs{trig},3);
end

%%
% condition if no avg is required
if ps_k==1
    group     = [];
    pseudo_tr = epochs;
    n_tr      = [];
    return
end

% condition if num of trials to avg are given
if isempty(ps_k) == 0
    
    mintr = min(mintr);
    temp  = mod(mintr,ps_k);
   
    if temp == 0
        n_tr = mintr/ps_k;
%         disp(['Making ',num2str(mintr/ps_k),' pseudotrials']);
    else
        n_tr = (mintr-temp)/ps_k;
%         disp(['Making ',num2str((mintr-temp)/ps_k),' pseudotrials']);
    end
end

% condition if num of pseudotrials reqd are given
if isempty(ps_n) == 0
    
    mintr = min(mintr);
    temp  = mod(mintr,ps_n);
%     disp(['Making ',num2str(ps_n),' pseudotrials']);
    ps_k  = (mintr-temp)/ps_n;
    temp  = mod(mintr,ps_k);
    n_tr  = ps_n;
end

% choose trials to avg further
for trig = 1:length(epochs)
    r = randperm(size(epochs{trig},3),(mintr-mod(mintr,ps_n)));
    if isempty(temp) == 0
        epochs{trig} = epochs{trig}(:,:,r);
    end
end

% pseudotrials
pseudo_tr=[];
for trig = 1:length(epochs)
    r     = randperm(size(epochs{1},3),size(epochs{1},3)); % create a unique order
    total_gr = (size(epochs{1},3)/ps_k);
    group{trig} = zeros(total_gr,ps_k); % create groups of n
    
    temp = mod(size(epochs{1},3),ps_k);
    
    if mod(temp,ps_k)~=0
        error('Enter n to divide all trials equally');
    else
        
        ps_tr=1; %reset group counter
        for num_group = 1:total_gr
            
            temp  = randperm(length(r),ps_k); %choose random n indices
            group{trig}(num_group,:) = r(temp); %make a group of values of those indices
            r(temp) = []; %remove the used values
            
            pseudo_tr{trig}(:,:,ps_tr) = mean(epochs{trig}(:,:,group{trig}(num_group,:)),3); %mean of the grouped trials
            ps_tr = ps_tr+1; %next group
        end
    end
end

end
