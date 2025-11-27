function [acc_mean, auc_mean] = fn_svm_decode_libsvm_src(param,roi_timeseries,comparison)

% Binary decoding using SVM within an roi / sphere at source level
% Written by Siddharth Talwar
% Last edited on 26-11-2025
% Inputs/param initialized to call fn_computesourcetime_search
% 1) param.n_category   : number of categories
% 2) param.n_iterations : number of iterations
% 3) param.pseudo_k     : number of trials to be averaged to make pseudotr
% 4) param.pseudo_n     : number of pseudotrials reqd for decoding
%                         Put EITHER pseudo_k OR pseudo_n.
%
% roi_timeseries        : source-time seriess from
%                         fn_computesourcetime_search
%
% comparison            : computed in fn_computesourcetime_search

% Outputs-
% 1) acc_mean           : computes accuracy of classification between each
%                         comparison. Averages all iterations. 
%                         structure = comparison x time
% 1) auc_mean           : computes AUC of classification between each
%                         comparison. Averages all iterations. 
%                         structure = comparison x time
% NOTES
% -This function is nested within another function to compute source time
% series (fn_computesourcetime_search)
% - Pseudotrials are created using an inbuilt function.
% - Leave one trial out cross validation implemented.
% - Uses LIBSVM library for decoding

%% Source time binary decoding

acc_mean = zeros(length(roi_timeseries),length(1 :...
    size(roi_timeseries.trial{1},2)));
auc_mean = acc_mean;

% prepare datasets
acc_comp_iter = zeros(size(comparison,1),...
    length(1:...
    size(roi_timeseries.trial{1},2)),...
    param.n_iterations);
auc_comp_iter = acc_comp_iter;

% iteration parfor loop
for iter = 1 : param.n_iterations
    
    % Make pseudotrials
    [pseudo_tr,~] = make_pseudo(roi_timeseries,param.pseudo_k,...
        param.pseudo_n);

    % divide data as per predefined categories
    sounds_in_cat=cell(1,param.n_category); % indices of sounds in each category
    
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
    dataset = cell(1,param.n_category); dataset_re=cell(1,param.n_category);
    for ii = 1 : param.n_category
        dataset{ii} = cat(3,pseudo_tr{sounds_in_cat{ii}});
        
        % reshape dataset
        dataset_re{ii} = reshape(dataset{ii},...
            [size(dataset{ii},1)*size(dataset{ii},2),...
            size(dataset{ii},3)] )';
    end
    
    % prepare datasets
    acc_comp = zeros(size(comparison,1),...
        length(1:...
        size(roi_timeseries.trial{1},2)));
    auc_comp = acc_comp;
   
    % main loop for each comparison
    for comp = 1:size(comparison,1)
      
        % labels of each category
        sz = size(dataset_re{comparison(comp)},1);
        target_label=[ones(sz,1)*-1;...
            ones(sz,1)];
        
        % prepare dataset for C.V. as per current comparison  
        ds          = cat(1,dataset_re{comparison(comp,1)},...
            dataset_re{comparison(comp,2)});
        indx        = 1:length(target_label);
        sampleindx  = 1 : roi_timeseries.num_s :...
            size(dataset_re{comparison(comp)},2);
        targets = [-1 1];
        
        % check if number of trials are equal across 2 categories
        if length(target_label==1) ~= length(target_label==2)
            error(['Unequal number of trials in category ',...
                num2str(comparison(comp,1)),...
                ' and ',num2str(comparison(comp,2))]);
        end
        
        % make datastructures for plt
        acc_time = zeros(1,length(1 :...
            size(roi_timeseries.trial{1},2)));
        auc_time = acc_time;
        
        
        % compute for each time point
        jj=1;

        for time = 1 : size(roi_timeseries.trial{1},2)
            
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
                    sampleindx(time)+(roi_timeseries.num_s-1));
                
                samples_test  = ds(testindx,...
                    sampleindx(time) :...
                    sampleindx(time)+(roi_timeseries.num_s-1));
               
                model = svmtrain(target_train, samples_train,'-q');
                [label, auc_fold(fold), ~] = do_binary_predict...
                    ([-1;1], samples_test, model);
                 
                % compute accuracy
                count = 0;
                for ii = 1:length(label)
                    if targets(ii) == label(ii)
                        count = count+1;
                    end
                end
                acc_fold(fold) = count/2;
             
            end
             
            % get accuracy for each time point
            acc_time(jj)   = mean(acc_fold); % mean(folds)
            auc_time(jj)   = mean(auc_fold); % mean(folds)
            
            jj=jj+1;
        end
        
        % get accuracy and weights for each comparison
        
        acc_comp(comp,:)   = acc_time;
        auc_comp(comp,:)   = auc_time;
         
    end
   
    acc_comp_iter(:,:,iter) = acc_comp;
    auc_comp_iter(:,:,iter) = auc_comp;
 
end
acc_mean = mean(acc_comp_iter,3);
auc_mean = mean(auc_comp_iter,3);

