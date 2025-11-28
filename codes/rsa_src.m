% Representation Similarity Analysis in source space
% Written by Siddharth Talwar
% Last edited on 26-11-2025
% Loads all parrticipant specific decoding AUC, converts them into DSMs,
% correlates with other external models, perform sign permutation testing
% with cluster correction and outputs a .nii file.
% Uses SPM12 and bspmview to view the output.

%% Load data and make DSMs
clear
close all
clc
parent  = 'xxx/Temporal-dynamics-of-natural-sound-representations-in-the-brain-of-sighted-and-blind';

cd(parent);
deriv   = fullfile(parent,'derivatives');

addpath 'xxx/fieldtrip-20191024';
addpath('xxx/Temporal-dynamics-of-natural-sound-representations-in-the-brain-of-sighted-and-blind/codes');
ft_defaults

%folder info
param = [];
param.factor       = 0;
param.n_category   = 24;
param.iter         = 100;
param.pseudo_n     = 3;
param.trial_radius = 0;
pop=0;
acc=[];

timewin = [0:0.02:1.25];
timewin(end) = [];

tp = 5; %timepoints of interest (see timewin)

for popn = {'Blind','Control'}
    pop = pop+1;
    popn      = char(popn);
    dir_deriv = dir(fullfile(deriv,popn));
    names     = {dir_deriv([dir_deriv.isdir]).name};
    dir_deriv(ismember(names,{'.','..'})) = [];
    
    sub_ii = 1;
   for sub = 1:length(dir_deriv)
        
        if ismember(sub,[1 8:15]) && pop==1
            load(fullfile(deriv,popn,['sub',sprintf('%03d',sub)],...
                ['binarydecoding_svm_',num2str(param.n_category),'cat'],...
                ['iter',num2str(param.iter),'_pt',num2str(param.pseudo_n)],...
                'acc_libsvm_src_tmplatemodel.mat'),'auc_mean');
            
        elseif ismember(sub,[2:7 16:18]) && pop==1
            
            load(fullfile(deriv,popn,['sub',sprintf('%03d',sub)],...
                ['binarydecoding_svm_',num2str(param.n_category),'cat'],...
                ['iter',num2str(param.iter),'_pt',num2str(param.pseudo_n)],...
                'acc_libsvm_src_indivmodel.mat'),'auc_mean');
            
        elseif ismember(sub,[1 5 7:10 12:14 16:18]) && pop==2
            load(fullfile(deriv,popn,['sub',sprintf('%03d',sub)],...
                ['binarydecoding_svm_',num2str(param.n_category),'cat'],...
                ['iter',num2str(param.iter),'_pt',num2str(param.pseudo_n)],...
                'acc_libsvm_src_tmplatemodel.mat'),'auc_mean');
            
        elseif ismember(sub,[2:4 6 11 15]) && pop==2
            load(fullfile(deriv,popn,['sub',sprintf('%03d',sub)],...
                ['binarydecoding_svm_',num2str(param.n_category),'cat'],...
                ['iter',num2str(param.iter),'_pt',num2str(param.pseudo_n)],...
                'acc_libsvm_src_indivmodel.mat'),'auc_mean');
        end
      
        
        for roi = 1:size(auc_mean,1)
                
                [srcdsm{pop}(roi,:,:,:,sub_ii), ~, comparison] =...
                    make_dsm(squeeze(mean(auc_mean(roi,:,tp),3)'),24,0); 

            % dsm = roi x ncat x ncat x timepoints x subjects
        end
    
        sub_ii
        sub_ii = sub_ii+1;
        
    end
    
end

%% RSA section by section

% close all

parent = '/media/siddharth/DATA/CPP/Projects/Aud_Cat/codes/Temporal-dynamics-of-natural-sound-representations-in-the-brain-of-sighted-and-blind';

eeg    = load(fullfile(parent,'others','dsm','eeg_dsm_libsvm.mat'));

mtf    = load(fullfile(parent,'others','dsm','mtf_dsm.mat'));
categ8 = load(fullfile(parent,'others','dsm','8cat_dsm.mat'));
beh    = load(fullfile(parent,'others','dsm','beh_dsm.mat'));
dnn    = load(fullfile(parent,'others','dsm','dnn_dsm.mat'));
sem    = load(fullfile(parent,'others','dsm','w2v_dsm.mat'));

model1 = srcdsm; %always src
model2 = categ8.dsm; %one model at a time
model3 = mtf.dsm; % to partial out

% Corr with eeg(model1) and model2, controlling for model 3
cor=[]; pop=1; rho=[];

if iscell(model2)==0
    y          = tril(model2 + 99,-1);
    temp_model = reshape(y(y>0),[((size(model2,1)*size(model2,1))/2 ) - (size(model2,1)/2),1]);
    temp_model = (temp_model-99);
end

y          = tril(model3 + 99 ,-1);
temp_model2 = reshape(y(y>0),[((size(model3,1)*size(model3,1))/2 ) - (size(model3,1)/2),1]);
temp_model2 = (temp_model2-99);

for popn = {'Blind','Control'}
    popn = char(popn);
    
    for roi = 1 : size(model1{pop},1)
        
        for sub = 1 : size(model1{pop},5) % for every subject
            
            if iscell(model2)==1
                
                %%%% change this based on model
%                 y          = tril(mean(1-model2{pop}(:,:,:)+7,3) + 99,-1); %beh
                y          = tril(mean(model2{pop}(:,:,:),3) + 99,-1); % sem
                temp_model = reshape(y(y>0),[((size(model2{pop},1)*size(model2{pop},1))/2 ) - (size(model2{pop},1)/2),1]);
                temp_model = (temp_model-99);
            end
            
            for tp = 1 : size(model1{pop},4) % every eeg timepoint
                
                y             = tril( squeeze( model1{pop}(roi,:,:,tp,sub)) + 99,-1);
                temp_eegmodel = reshape(y(y>0),[((size(model1{pop},2)*size(model1{pop},2))/2 ) - (size(model1{pop},2)/2),1]);
                temp_eegmodel = temp_eegmodel-99;
                
                cor{pop}(roi,tp,sub) = corr(temp_eegmodel,temp_model,'type','Spearman');
                
                rho{pop}(roi,tp,sub) = partialcorr(temp_eegmodel,temp_model,temp_model2,'type','Spearman');
                 
            end
        end
        
    end
    pop = pop+1;
end
clc

%% new stats
tp = 1;
vec1 = squeeze(mean(rho{1}(:,tp,:),2));
vec2 = squeeze(mean(rho{2}(:,tp,:),2));
permnum = 5000;
thr = 0.05;
[src,pval] = fn_signperm_src(vec1,[],permnum,thr);

%% open in bspmview
% src.pow(aparc<=0) = 0;
addpath('/media/siddharth/DATA/Toolbox/fieldtrip-20191024/external/spm12/templates');

mri           = ft_read_mri('T1.nii');
cfg           = [];
cfg.parameter = 'all';
srcimg        = ft_sourceinterpolate(cfg,src,mri);

addpath('xxx/spm12');
addpath('xxx/bspmview-master');

cfg           = [];
cfg.filename  = 'xxx';
cfg.filetype  = 'nifti';
cfg.parameter = 'pow';
ft_sourcewrite(cfg,srcimg);

% open in bspmview