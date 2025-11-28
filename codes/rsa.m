%% Represenational Similarity Analysis
% Written by Siddharth Talwar
% Last edited on 26-11-2025
% This script loads EEG and different external models, performs statistics
% using sign permutation testing with cluster-wise correction and plots the
% results. Insert model names as needed.

clear
clc
close all

parent = 'xxx/Temporal-dynamics-of-natural-sound-representations-in-the-brain-of-sighted-and-blind';

eeg    = load(fullfile(parent,'others','dsm','eeg_dsm_libsvm.mat'));

mtf    = load(fullfile(parent,'others','dsm','mtf_dsm.mat'));
pitch  = load(fullfile(parent,'others','dsm','pitch_dsm.mat'));
hnr    = load(fullfile(parent,'others','dsm','hnr_dsm.mat'));
cog    = load(fullfile(parent,'others','dsm','cog_dsm.mat'));

categ8 = load(fullfile(parent,'others','dsm','8cat_dsm.mat'));
categ4 = load(fullfile(parent,'others','dsm','4cat_dsm.mat'));
categ2 = load(fullfile(parent,'others','dsm','2cat_dsm.mat'));
beh    = load(fullfile(parent,'others','dsm','beh_dsm.mat'));

dnn    = load(fullfile(parent,'others','dsm','dnn_dsm.mat'));

sem    = load(fullfile(parent,'others','dsm','w2v_dsm.mat'));

%%% Insert models here

model1        = eeg.avg_eegdsm; %always EEG
model2        = {mtf.dsm categ8.dsm beh.dsm_sr sem.dsm_w2v dnn.dsm(:,:,4)}; % models to correlate
model2_labels = {'mtf','8-category','behavior','semantic w2v','dnn l4'};
model3        = mtf.dsm; % model to partialout

time = -0.5 + (2.5/size(model1{1},3) ):2.5/size(model1{1},3) :2;

%% Corr with eeg(model1) and model2, controlling for model 3

cor=[]; pop=1; rho=[];

for mdl = 1:length(model2)
    
    if isnumeric(model2{mdl}) == 1 % to correlate
        y          = tril(model2{mdl} + 99,-1);
        temp_model2 = reshape(y(y>0),[numel(y(y>0)),1]);
        temp_model2 = (temp_model2-99);
    end
    
    if isnumeric(model3) == 1 % to partial out
        for ii = 1:size(model3,3)
            y           = tril(model3(:,:,ii) + 99 ,-1);
            temp_model3(:,ii) = reshape(y(y>0),[numel(y(y>0)),1]);
            temp_model3(:,ii) = (temp_model3(:,ii)-99);
        end
    end
    
    for pop = 1:length(model1)
        
        for sub = 1 : size(model1{pop},4) % for every subject
            
            if iscell(model2{mdl}) == 1 % to correlate
                
                if isequal(model2{mdl},beh.dsm_sr) == 1 % subject specific models
                    y          = tril(mean(1-model2{mdl}{pop}(:,:,sub)+7,3) + 99,-1); %for indiv dsm
                elseif isequal(model2{mdl},sem.dsm_w2v) == 1 % averaged across subjects
                    y          = tril(mean(model2{mdl}{pop},3) + 99,-1); %sem
                end
                
                temp_model2 = reshape(y(y>0),[numel(y(y>0)),1]);
                temp_model2 = (temp_model2-99);
                
            end
            
            if iscell(model3) == 1 % to partial out
                if isequal(model2{mdl},beh) == 1 % subject specific models
                    y          = tril(mean(1-model2{mdl}{pop}(:,:,sub)+7,3) + 99,-1); %for indiv dsm
                elseif isequal(model2{mdl},sem) == 1 % averaged across subjects
                    y          = tril(mean(model2{mdl}{pop},3) + 99,-1); %sem
                end
                
                temp_model3 = reshape(y(y>0),[numel(y(y>0)),1]);
                temp_model3 = (temp_model3-99);
                
            end
            
            for tp = 1 : size(model1{pop},3) % every eeg timepoint
                
                tmp = model1{pop}(:,:,tp,sub);
                y             = tril(tmp + 99,-1);
                temp_eegmodel = reshape(y(abs(y)>0),[numel(y(abs(y)>0)),1]);
                temp_eegmodel = (temp_eegmodel-99);
                
                cor{pop,mdl}(tp,sub) = corr(temp_eegmodel,temp_model2,...
                    'type','Spearman');
                
                rho{pop,mdl}(tp,sub) = partialcorr(temp_eegmodel,temp_model2,...
                    temp_model3,'type','Spearman');
                
            end
        end
        
        % mean and std dev
        cor_mean{pop,mdl}     = mean(cor{pop,mdl},2);
        cor_std{pop,mdl}  = std(squeeze(cor{pop,mdl})')./...
            sqrt(size(cor{pop,mdl},2));
        
        rho_mean{pop,mdl}     = mean(rho{pop,mdl},2);
        rho_std{pop,mdl}  = std(squeeze(rho{pop,mdl})')./...
            sqrt(size(rho{pop,mdl},2));
        
        % for plotting std dev
        curve1 = cor_mean{pop,mdl}' + cor_std{pop,mdl};
        curve2 = cor_mean{pop,mdl}' - cor_std{pop,mdl};
        xcor{pop,mdl} = [time, fliplr(time)];
        inBetweencor{pop,mdl} = [curve1, fliplr(curve2)];
        
        curve1 = rho_mean{pop,mdl}' + rho_std{pop,mdl};
        curve2 = rho_mean{pop,mdl}' - rho_std{pop,mdl};
        xrho{pop,mdl} = [time, fliplr(time)];
        inBetweenrho{pop,mdl} = [curve1, fliplr(curve2)];
        
    end
end
%% Stats on correlations

sigind   = [];
startind = 128; % only poststim
permnum  = 100;
thr      = 0.05;

for pop = 1 : 2 %3rd is the difference b/w 2 popn
    vec{pop} = cor{pop,mdl}(startind:end,:);
    
    for mdl = 1:length(model2)
        % sign permutation test
        [pval{pop,mdl},nulldist{pop,mdl},sigind{pop,mdl},cl_thr{pop,mdl}] = ...
            fn_signperm(vec{pop} ,[], 1,permnum,thr);
    end
    
end
% difference between group 1 and 2
for mdl = 1:length(model2)
    [pval{3,mdl},nulldist{3,mdl},sigind{3,mdl},cl_thr{3,mdl}] = ...
        fn_signperm(cor{1,mdl}(startind:end,:) ,cor{2,mdl}(startind:end,:), 2,permnum,thr);
end

%% Plot correlations

colr{1} = [255 158 74]/256;
colr{2} = [105 170 173]/256;

figure
for mdl = 1:length(model2)
    subplot(3,2,mdl);
    hold on;
    for pop = 1:2
        plot(time,(cor_mean{pop,mdl}),'LineWidth',2,'Color',colr{pop},...
            'Linestyle','-');     
    end
    
    for pop= 1:2
        patch(xcor{pop,mdl}, inBetweencor{pop,mdl}, colr{pop},'facealpha',0.25,'edgealpha',0.1,'edgecolor',colr{pop});
    end

    plot(time,zeros(1,size(time,2)),'k--');
    
    plot(time(startind + sigind{3,mdl}),repmat(-0.07,[1, length(sigind{3,mdl})]),'o',...
        'Color',[0 0 0]/256,'MarkerSize',5);
    plot(time(startind + sigind{1,mdl}),repmat(-0.08,[1, length(sigind{1,mdl})]),'o',...
        'Color',[255 158 74]/256,'MarkerSize',5);
    plot(time(startind + sigind{2,mdl}),repmat(-0.09,[1, length(sigind{2,mdl})]),'o',...
        'Color',[105 170 173]/256,'MarkerSize',5);
    
    xlabel('time', 'FontSize', 12);
    legend('Blind','Sighted');
    legend('boxoff');
    ylabel('Spearman correlation', 'FontSize', 16);
    xticks([-0.5:0.25:2]);
    plot(time,zeros(1,size(time,2)),'k--');
    xlim([-0.25 1.5]);
    
    title(['EEG and ',model2_labels{mdl}],'FontSize', 14);
    
    ax = gca;
    ax.FontSize = 12;
    set(gca, 'FontName', 'Avenir');
    set(gca,'linewidth',2)
    
    box off;
    hold off;
end