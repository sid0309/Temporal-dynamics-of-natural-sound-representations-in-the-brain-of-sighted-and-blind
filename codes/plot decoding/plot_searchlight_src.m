clear
close all
clc
parent  = '/media/siddharth/DATA/CPP/Projects/Aud_Cat/codes/Temporal-dynamics-of-natural-sound-representations-in-the-brain-of-sighted-and-blind';

cd(parent);
deriv   = fullfile(parent,'derivatives');

addpath '/media/siddharth/DATA/Toolbox/fieldtrip-20191024';
ft_defaults

addpath(genpath('/media/siddharth/DATA/Toolbox/fieldtrip-20191024/template'));

param = [];
param.factor       = 0;
param.n_category   = 24;
param.iter         = 100;
param.pseudo_n     = 3;
param.trial_radius = 0;

timewin = [0:0.02:1.25];% 8:15 16:21 22:27 29:35 37:41 43:45 49:50

% load participant data
acc=[];pop=0;
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
        
        acc{pop}(:,:,sub_ii) = mean(auc_mean,2); % source x time x sub
        sub_ii = sub_ii+1;
        
    end
    
end

%% localize the timewindow

timewin(end) = [];
tp = 49:50;
mat1 = squeeze(mean(acc{1}(:,tp,:),2))-0.5;
mat2 = squeeze(mean(acc{2}(:,tp,:),2))-0.5;
permnum = 5000;
thr = 0.05;
[src,pval] = fn_signperm_src(mat1,mat2,permnum,thr);

%% Prepare img for bidspmview

addpath('/media/siddharth/DATA/Toolbox/fieldtrip-20191024/external/spm12/templates')

mri = ft_read_mri('T1.nii');
cfg=[];
cfg.parameter = 'all';
srcintimg = ft_sourceinterpolate(cfg,src,mri);

addpath('/media/siddharth/DATA/Toolbox/spm12');
addpath('/media/siddharth/DATA/Toolbox/bspmview-master');

cfg=[];
cfg.filename  = 'tmp';
cfg.filetype  = 'nifti';
cfg.parameter = 'pow';
ft_sourcewrite(cfg,srcintimg);

% bspmview to plot