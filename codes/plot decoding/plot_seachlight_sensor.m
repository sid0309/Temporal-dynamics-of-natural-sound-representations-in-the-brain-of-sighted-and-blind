%% Plotting Sensor Searchlight
% Written by Siddharth Talwar
% Last edited on 26-11-2025
% This scripts calls the output from fn_svm_decode_libsvm_chan from
% individual paticipant folders. Implements 2-D sign permutation test
% across channels.

clear
clc
parent  = 'xxx/Temporal-dynamics-of-natural-sound-representations-in-the-brain-of-sighted-and-blind'; %change path

deriv   = fullfile(parent,'derivatives');
addpath 'xxx/fieldtrip-20191024'
ft_defaults

%input folder info
param.n_category   = 24;% category folder name
param.iter         = 100;% num iterations
param.factor       = 0; % smoothing factor
param.pseudo_n     = 3; % number of pseudotrials
param.popn         = {'Blind','Control'}; % folder names
param.excludecat   = [];

%% Plot blind vs control accuracy comparison

acc     = {};
accmeansub = {};
pop = 0;

for popn = param.popn
    pop = pop+1;
    popn      = char(popn);
    dir_deriv = dir(fullfile(deriv,popn));
    names     = {dir_deriv([dir_deriv.isdir]).name};
    dir_deriv(ismember(names,{'.','..'})) = [];
    sub_ii = 1;
    
    for sub = 1:length(dir_deriv)
        
        load(fullfile(deriv,popn,['sub',sprintf('%03d',sub)],...
            ['binarydecoding_svm_',num2str(param.n_category),'cat'],...
            ['iter',num2str(param.iter),...
            '_pt',num2str(param.pseudo_n)],'acc_libsvm_chan.mat'),...
            'auc_chan_comp_tp','comparison');
        
        auc_chan_comp_tp = squeeze(mean(auc_chan_comp_tp,2)); % mean of comparisons
        acc{pop}(:,:,sub_ii) = auc_chan_comp_tp; % chan x time x sub
        
        sub_ii = sub_ii+1;
        
    end
end

timesl = -0.25 + (1.75/size(acc{1},2)):1.75/size(acc{1},2):1.5; % always 640

%% layout

%load dummy var
load(fullfile(parent,'preprocess','Control','sub001','preproc.mat'));
cfg = [];
cfg.layout = '/media/siddharth/DATA/Toolbox/fieldtrip-20191024/template/layout/biosemi128.lay';
lay = ft_prepare_layout(cfg);

%% stats and plot

%put different timepoints to plot them individually. Time is averaged
%between these samples.

tp = 56:68;%33:45; %52:69 70:84 %87:101 %104:121 125:137

vec{1} = squeeze(mean(acc{1}(:,tp,:),2))-0.5; 
vec{2} = squeeze(mean(acc{2}(:,tp,:),2))-0.5;

for popn = 1:2
    
    sigind{popn} = [];
    
    % sign permutation test
    disp('Sign permutation test');
    [pval,nulldist,sigind{popn}] =  fn_signperm_chan(vec{popn},[],1,1000,0.05,lay);
    
    topo             = preprocdata;
    topo             = rmfield(topo,'trialinfo');
    topo             = rmfield(topo,'sampleinfo');
    topo.trial       = {};
    topo.trial{1}    = squeeze(mean(mean(acc{popn}(:,tp,:),2),3));
    topo.time(2:end) = [];
    topo.time{1}     = timesl(tp(1));
    tl               = ft_timelockanalysis([],topo);
    
    cfg              = [];
    cfg.parameter    = 'avg';
    cfg.layout       = lay;
    cfg.comment      = 'no';
    cfg.marker       = 'off';
    cfg.highlight    = 'on';
    cfg.highlightsymbol = '.';
    cfg.highlightsize   = 15;
    cfg.highlightchannel = sigind{popn};
    cfg.zlim      = [0.5 max(squeeze(mean(mean(acc{1}(:,tp,:),2),3)))];
    figure;
    ft_topoplotER(cfg, tl);
    colorbar;
    title([num2str(floor(timesl(tp(1))*1000)),' to ',...
        num2str(floor(timesl(tp(end))*1000)),' ms'],'FontSize',20);
    set(gca, 'FontName', 'Arial');
    colormap((colormap('parula')));
    
end

%% blind - sighted
popn=3;
topo.trial{1} = squeeze(mean(mean(acc{1}(:,tp,:),2),3)) -...
    squeeze(mean(mean(acc{2}(:,tp,:),2),3));

tl = ft_timelockanalysis([],topo);

disp('Sign permutation test');
[pval,nulldist,sigind{popn}] =  fn_signperm_chan(vec{1},vec{2},2,1000,0.05,lay);

figure;
cfg.highlightchannel = sigind{popn};%
cfg.zlim      = [0 0.0455];
ft_topoplotER(cfg, tl);
colorbar;
title([num2str(floor(timesl(tp(1))*1000)),' to ',...
    num2str(floor(timesl(tp(end))*1000)),' ms'],'FontSize',20);
set(gca, 'FontName', 'Arial');
