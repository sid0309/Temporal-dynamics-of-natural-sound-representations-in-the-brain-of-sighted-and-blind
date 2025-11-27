% Make eeg dsm from precomputed decoding AUC for sensor and source levels
% Written by Siddharth Talwar
% Last edited on 26-11-2025

clear
clc
close all

parent  = '/media/siddharth/DATA/CPP/Projects/Aud_Cat/codes/Temporal-dynamics-of-natural-sound-representations-in-the-brain-of-sighted-and-blind';

cd(parent);
deriv   = fullfile(parent,'derivatives');

param              = [];
param.factor       = 0; % to average neighbouring bins
param.n_category   = 24; % categorical 24,8,4
param.pseudo_n     = 3; % number of pseudotrials used
param.iter         = 100; % iterations

%% EEG sensor level DSM

c=[]; ii=1;
for popn = {'Blind','Control'}
    dsm = [];
    avg_dsm = [];
    popn = char(popn);
    
    dir_deriv = dir(fullfile(deriv,popn));
    names       = {dir_deriv([dir_deriv.isdir]).name};
    dir_deriv(ismember(names,{'.','..'}))=[];
    
    eegdsm{ii} = [];
    avg_eegdsm{ii} = [];
    
    % compute mean accuracies across all iterations
    for sub = 1:length(dir_deriv)
        
        load(fullfile(deriv,popn,['sub',sprintf('%03d',sub)],...
            ['binarydecoding_svm_',num2str(param.n_category),'cat'],...
            ['iter',num2str(param.iter),...
            '_pt',num2str(param.pseudo_n)],'acc_libsvm.mat'),'auc_comp_iter');
        
        dsm = [];avg_dsm = [];
        for iter = 1:param.iter
            
            [dsm(:,:,:,iter), avg_dsm(:,:,:,iter), comparison] =...
                make_dsm(auc_comp_iter(:,:,iter),param.n_category,param.factor);
            % dsm = ncat x ncat x timepoints x iterations
        end
        
        eegdsm{ii}(:,:,:,sub) = mean(dsm,4); % average iterations for sub
        avg_eegdsm{ii}(:,:,:,sub) = mean(avg_dsm,4); % average iterations
        
    end
    ii=ii+1;
end

time = -0.5 + (2.5/640):2.5/640:2;
timesmooth = time(1+param.factor:length(time)-param.factor);

save(fullfile(parent,'others','dsm','eeg_dsm_libsvm.mat'),'eegdsm','avg_eegdsm',...
    'time','timesmooth','param');

%% Source level DSM

c=[];
ii=1;
for popn = {'Blind','Control'}
    
    popn = char(popn);
    
    dir_deriv = dir(fullfile(deriv,popn));
    names       = {dir_deriv([dir_deriv.isdir]).name};
    dir_deriv(ismember(names,{'.','..'}))=[];
    
    eegdsm{ii} = [];
    avg_eegdsm{ii} = [];
    
    % compute mean accuracies across all iterastions
    for sub = 1:length(dir_deriv)

        if ismember(sub,[1 8:15]) && ii==1
            load(fullfile(deriv,popn,['sub',sprintf('%03d',sub)],...
                ['binarydecoding_svm_',num2str(param.n_category),'cat'],...
                ['iter',num2str(param.iter)],'sourcetime',...
                'auc_sourcetime_svd_tmplatefwdmodel_glasser_-0.1s.mat'),'auc_mean');
            
        elseif ismember(sub,[2:7 16:18]) && ii==1
            
            load(fullfile(deriv,popn,['sub',sprintf('%03d',sub)],...
                ['binarydecoding_svm_',num2str(param.n_category),'cat'],...
                ['iter',num2str(param.iter)],'sourcetime',...
                'auc_sourcetime_svd_indivfwdmodel_glasser_-0.1s.mat'),'auc_mean');
            
        elseif ismember(sub,[1 5 7:10 12:14 16:18]) && ii==2
            load(fullfile(deriv,popn,['sub',sprintf('%03d',sub)],...
                ['binarydecoding_svm_',num2str(param.n_category),'cat'],...
                ['iter',num2str(param.iter)],'sourcetime',...
                'auc_sourcetime_svd_tmplatefwdmodel_glasser_-0.1s.mat'),'auc_mean');
            
        elseif ismember(sub,[2:4 6 11 15]) && ii==2
            load(fullfile(deriv,popn,['sub',sprintf('%03d',sub)],...
                ['binarydecoding_svm_',num2str(param.n_category),'cat'],...
                ['iter',num2str(param.iter)],'sourcetime',...
                'auc_sourcetime_svd_indivfwdmodel_glasser_-0.1s.mat'),'auc_mean');
        end
        
        for roi = 1:size(auc_mean,1)
            
            [srcdsm{ii}(roi,:,:,:,sub), avg_srcdsm{ii}(roi,:,:,:,sub), comparison] =...
                make_dsm(squeeze(auc_mean(roi,:,:)),param.n_category,param.factor);
            % eegdsm = ncat x ncat x timepoints x subjects x iterations
        end
                
    end
    ii=ii+1;
end

time = -0.5 + (2.5/320):2.5/320:2;
timesmooth = time(1+param.factor:length(time)-param.factor);

save(fullfile(parent,'others','dsm','srcdsm.mat'),'srcdsm','avg_srcdsm',...
    'time','timesmooth','param');

%% EEG sensor searchlight DSMs

c=[]; ii=1;
for popn = {'Blind','Control'}
    dsm = [];
    avg_dsm = [];
    popn = char(popn);
    
    dir_deriv = dir(fullfile(deriv,popn));
    names       = {dir_deriv([dir_deriv.isdir]).name};
    dir_deriv(ismember(names,{'.','..'}))=[];
    
    eegdsm{ii} = [];
    avg_eegdsm{ii} = [];
    
    % compute mean accuracies across all iterations
    for sub = 1:length(dir_deriv)
        
        load(fullfile(deriv,popn,['sub',sprintf('%03d',sub)],...
            ['binarydecoding_svm_',num2str(param.n_category),'cat'],...
            ['iter',num2str(param.iter),...
            '_pt',num2str(param.pseudo_n),...
            '_tr',num2str(param.trial_radius)],'acc_chan_libsvm.mat'),'auc_comp_iter_chan');
        
        dsm = [];avg_dsm = [];
        for chan = 1:128
            
            [dsm(:,:,:,chan), ~, comparison] =...
                make_dsm(squeeze(auc_chan_comp_tp(chan,:,:)),param.n_category,param.factor);
            % dsm = ncat x ncat x timepoints x chan x sub
        end
        eegdsm{ii}(:,:,:,:,sub) = dsm; % average iterations for sub
        
        sub
    end
    ii=ii+1;
end

time = -0.25 + (1.75/225):1.75/225:1.5; % always 640
timesmooth = time(1+param.factor:length(time)-param.factor);

save(fullfile(parent,'others','dsm','eegdsm_libsvm_sensorsearch.mat'),'eegdsm',...
    'time','timesmooth','param');

