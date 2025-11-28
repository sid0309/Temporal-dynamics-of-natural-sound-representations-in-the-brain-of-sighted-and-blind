function fn_svm_decode_libsvm(param)

% Binary decoding using SVM.
% Written by Siddharth Talwar
% Last edited on 26-11-2025

% Inputs -
% 1) param.n_category   : number of categories
% 2) param.timewindow   : [begin end] time in data, [] for all timepoints
% 3) param.fs           : sampling frequency of the dataset
% 4) param.ds           : downsample time by a factor (2,4,8 etc)
% 5) param.n_iterations : number of iterations (from creating pseudotrials
%                         to decoding.
% 6) param.pseudo_k     : number of trials to be averaged to make pseudotr
% 7) param.pseudo_n     : number of pseudotrials reqd for decoding
%                         Put EITHER pseudo_k OR pseudo_n.
% 8) param.popn         : population name / folder name
% 9) param.system       : change paths for different system: monster or pc

% Outputs-
% 1) acc_comp_iter   : computes accuracy of classification between each
%                          comparison. comparison x time x iter
% 2) auc_comp_iter   : computes AUC of classification between each
%                          comparison. comparison x time x iter
% 3) comparison:     : indexed list of all binary comparisons
% 4) param           : input parameters 

% NOTES
% - Uses fieldtrip toolbox to perform few functions before deocoding
% - Uses LIBSVM library for decoding
% - All outputs are saved in the derivatives folder, in the corresponding
%   analysis folder and sub-folder for each subject
% - The function uses parallel toolbox. Make sure there are enough workers.
% - Pseudotrials are created using an inbuit function.
% - leave one trial out cross validation implemented.

%% Paths (change)

% change
parent  = 'xxx/Temporal-dynamics-of-natural-sound-representations-in-the-brain-of-sighted-and-blind'; %change
addpath 'xxx/fieldtrip-20191024' %fieldtrip change
ft_defaults
addpath('xxx/libsvm-master/matlab') %libsvm change

% don't change
cd(parent);
preproc = fullfile(parent,'preprocess');
deriv   = fullfile(parent,'derivatives');
dir_preproc = dir(fullfile(preproc,param.popn));
names       = {dir_preproc([dir_preproc.isdir]).name};
dir_preproc(ismember(names,{'.','..'})) = [];

%% Binary decoding and weights

for sub = param.sub
    
    clc
    % input: load input data
    data = load(fullfile(preproc,param.popn,...
        ['sub',sprintf('%03d',sub)],'preproc.mat'));
    
    % downsample
    if isempty(param.ds) ~= 1
        cfg=[];
        cfg.resamplefs = param.fs/param.ds;
        data.preprocdata = ft_resampledata(cfg,data.preprocdata);
    end
    
    % extract timewindow of interest
    if isempty(param.timewindow) ~= 1
        cfg              = [];
        cfg.toilim       = param.timewindow;
        data.preprocdata = ft_redefinetrial(cfg,data.preprocdata);
    end
    
    % num of pseudotrials to name a folder
    if isempty(param.pseudo_n) == 1
        [~,~,n_tr] = make_pseudo(data.preprocdata,param.pseudo_k,...
            param.pseudo_n);
    else
        n_tr = param.pseudo_n;
    end
    
    % make output folders
    add  = fullfile(deriv,param.popn,['sub',sprintf('%03d',sub)],...
        ['binarydecoding_svm_',num2str(param.n_category),'cat']);
    mkdir(add);

    dest = fullfile(add,['iter',num2str(param.n_iterations),...
        '_pt',num2str(n_tr)]);
    mkdir(dest);
    
    % compute number of comparisons
    comparison = flip(combnk(1:param.n_category,2));
    
    if comparison(1,1)~=1 % make sure the combinations start with the 1st
        comparison = flip(comparison);
    end
    
    % prepare datasets
    acc_comp_iter = zeros(size(comparison,1),...
        length(1:size(data.preprocdata.trial{1},2)),...
        param.n_iterations);
    
    auc_comp_iter = acc_comp_iter;
    
    
    clc
    disp(['subject num ',num2str(sub)]);
    
    % iteration parfor loop
    for iter = 1 : param.n_iterations
        tic
        disp(['iteration num ' ,num2str(iter)]);
        
        % Make pseudotrials
        [pseudo_tr,~] = make_pseudo(data.preprocdata,param.pseudo_k,...
            param.pseudo_n);
        
        % divide data as per predefined categories
        sounds_in_cat={}; % indices of sounds in each category
        
        if param.n_category == 2
            
            sounds_in_cat{1} = [1:6 13:18]; % sounds arranged like this
            sounds_in_cat{2} = [7:12 19:24];
            
        elseif ismember(param.n_category,[4 8 24]) == 1
            
            for ii = 1:param.n_category
                sounds_in_cat{ii} = (1:length(pseudo_tr) /...
                    param.n_category)...
                    + ((length(pseudo_tr) / param.n_category)*(ii-1));
            end
            
        else
            error('param.n_category can be 2,4,8 or 24');
        end
        
        % dataset_re = trials x (for timepoint 1 - 128 chan values,
        % then 128 val for next timepoint...)
        dataset = {}; dataset_re={};
        for ii = 1 : param.n_category
            dataset{ii} = cat(3,pseudo_tr{sounds_in_cat{ii}});
            
            % reshape dataset
            dataset_re{ii} = reshape(dataset{ii},...
                [size(dataset{ii},1)*size(dataset{ii},2),...
                size(dataset{ii},3)] )';
        end
        dataset = [];
        
        % prepare datasets
        acc_comp = zeros(size(comparison,1),...
            length(1:size(pseudo_tr{1},2)));
        
        auc_comp = acc_comp;
        
        % main loop for each comparison
        for comp = 1 : size(comparison,1)
            
            % labels of each category
            target_label=[ones(size(dataset_re{comparison(comp)},1),1)*-1;...
                ones(size(dataset_re{comparison(comp)},1),1)];
            
            % prepare dataset for C.V. as per current comparison
            ds          = cat(1,dataset_re{comparison(comp,1)},...
                dataset_re{comparison(comp,2)});
            indx        = 1:length(target_label);
            sampleindx  = 1 : size(data.preprocdata.trial{1},1) :...
                size(dataset_re{comparison(comp)},2);
            targets = [-1 1];
            
            % check if number of trials are equal across 2 categories
            if length(target_label==1) ~= length(target_label==2)
                error(['Unequal number of trials in category ',...
                    num2str(comparison(comp,1)),...
                    ' and ',num2str(comparison(comp,2))]);
            end
            
            % make datastructures for plt
            acc_time = zeros(1, 1 : length(size(pseudo_tr{1},2)));
            
            auc_time = acc_time;
            
            % compute for each time point
            jj=1;
            for time = 1 : size(pseudo_tr{1},2)
                
                % make datastructures for plt
                acc_fold = zeros(1,length(target_label)/2);
                auc_fold = acc_fold;
                
                % compute for each fold
                for fold = 1: length(target_label)/2
                    
                    % prepare training and testing data
                    testindx     = [fold length(target_label)/2+fold];
                    trainindx    = (~ismember(indx,testindx)==1);
                    target_train = target_label(trainindx==1);
                    
                    samples_train = ds(trainindx,...
                        sampleindx(time) :...
                        sampleindx(time)+(size(data.preprocdata.trial{1},1)-1));
                    
                    samples_test  = ds(testindx,...
                        sampleindx(time):...
                        sampleindx(time)+(size(data.preprocdata.trial{1},1)-1));
                    
                    model = svmtrain(target_train, samples_train,'-q');
                    [label, auc_fold(fold), ~] = do_binary_predict([-1;1], samples_test, model);
                    
                    % compute accuracy
                    count = 0;
                    for ii = 1:length(label)
                        if targets(ii) == label(ii)
                            count = count+1;
                        end
                    end
                    acc_fold(fold) = count/2;
                    
                end
                
                % get accuracy and weights for each time point
                acc_time(jj)   = mean(acc_fold); % mean(folds)
                auc_time(jj)   = mean(auc_fold); % mean(folds)
                jj=jj+1;
            end
            
            % get accuracy and weights for each comparison
            acc_comp(comp,:)   = acc_time;
            auc_comp(comp,:)   = auc_time;
            
        end
        
        % put each iteration in its place
        acc_comp_iter(:,:,iter) = acc_comp;
        auc_comp_iter(:,:,iter) = auc_comp;
        
        disp(['Iteration ',num2str(iter),' finished']);
        toc
        
    end

    save(fullfile(dest,'acc_libsvm.mat'),'acc_comp_iter','auc_comp_iter','comparison','param');
    
    clear acc_fold wt_fold acc_comp wt_comp acc_comp_iter...
        wt_comp_iter pseudo_tr ds svmmdl label
end
