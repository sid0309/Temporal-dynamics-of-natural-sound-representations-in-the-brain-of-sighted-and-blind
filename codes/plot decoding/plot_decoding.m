%% Plotting timepoint by timepoint decoding across all sensors
% Written by Siddharth Talwar
% Last edited on 26-11-2025
% This scripts calls the output from fn_svm_decode_libsvm from
% individual paticipant folders. Implements 1-D sign permutation test.

clear
clc
parent  = 'xxx/Temporal-dynamics-of-natural-sound-representations-in-the-brain-of-sighted-and-blind'; %change path

deriv   = fullfile(parent,'derivatives');
addpath 'xxx/fieldtrip-20191024'
ft_defaults

%input folder info
param.n_category   = 24;%  category folder name
param.iter         = 100;% num iterations
param.factor       = 2; % smoothing factor (number of bins on either side)
param.pseudo_n     = 3; % number of pseudotrials
param.trial_radius = 0; % number of neighbourhood timepoints
param.popn         = {'Blind','Control'}; % {popn1,popn2}

%% Plot blind vs control accuracy comparison

auc     = {};
accmeansub = {};
pop = 0;

time = -0.5 + (2.5/640):2.5/640:2; % always 640
timesmooth = time(1+param.factor:length(time)-param.factor);

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
            '_pt',num2str(param.pseudo_n)],'acc_libsvm.mat'),...
            'auc_comp_iter','comparison');
        
        auc{pop}(:,:,sub_ii) = auc_comp_iter; % comp x time x sub

        sub_ii = sub_ii+1;
    end
    
    accstd{pop}      = std(squeeze(mean(auc{pop},1))'); % std dev
    accmeansub{pop}  = mean(auc{pop},3); % mean of sub: comp x time
    accmeancomp{pop} = mean(accmeansub{pop},1); % mean of comp: 1 x time
    
    % smoothen the accuracy graph by averaging window given by param.factor
    count = 1;
    
    for comp = 1 : size(accmeansub{pop},1)
        ii = 1;
        
        for tp = 1+param.factor : size(accmeansub{pop},2)-param.factor
            smooth_acc{pop}(count,ii,:) = mean(auc{pop}...
                (comp,tp-param.factor:tp+param.factor,:),2);
            ii=ii+1;
        end
        count = count+1;
    end
    smooth_meanacc{pop} = mean(smooth_acc{pop},3);
    smooth_meanaccstd{pop}  = std(squeeze(mean(smooth_acc{pop},1))')./...
        sqrt(size(smooth_acc{pop},3));
    
    % for plotting std dev
    curve1 = mean(smooth_meanacc{pop},1) + smooth_meanaccstd{pop};
    curve2 = mean(smooth_meanacc{pop},1) - smooth_meanaccstd{pop};
    x{pop} = [timesmooth, fliplr(timesmooth)];
    inBetween{pop} = [curve1, fliplr(curve2)];
    
    ii=ii+1;
end

%% Stats sign permutation test
% stats done after smoothing the accuracy by 2 bins on either side

startind = 128;
vec{1} = squeeze(mean(smooth_acc{1}(:,startind:end,:),1))-0.5;%onlypoststim
vec{2} = squeeze(mean(smooth_acc{2}(:,startind:end,:),1))-0.5;

thr = [0.05 0.05 0.05];

for pop = 1 : 2 %3rd is the difference b/w 2 popn
    
    disp('Sign permutation test');
    
    tic
    [pval,nulldist,sigind{pop}] = fn_signperm(vec{pop},[],1,5000,thr(pop));
    toc
    
end

pop=3;
tic
[pval,nulldist,sigind{pop}] = fn_signperm(vec{1},vec{2},2,5000,thr(pop));
toc

%% Plotting

disp('Plotting');

figure;
hold on;
plot(timesmooth,mean(smooth_meanacc{1}(:,:),1),'LineWidth',2,...
    'Color',[255 158 74]/256);

plot(timesmooth,mean(smooth_meanacc{2}(:,:),1),'LineWidth',2,...
    'Color',[105 170 173]/256);

patch(x{1}, inBetween{1}, [255 158 74]/256,'facealpha',0.25,...
    'edgealpha',0.1,'edgecolor',[255 158 74]/256);

patch(x{2}, inBetween{2},[105 170 173]/256,'facealpha',0.25,...
    'edgealpha',0.1,'edgecolor',[105 170 173]/256);

plot(timesmooth,repmat(0.5,...
    [(size(auc_comp_iter,2) - (param.factor*2)),1]),'k--');

plot(timesmooth(startind-1 + sigind{1}),repmat(0.48,...
    [1, length(sigind{1})]),'o','Color',[255 158 74]/256,'MarkerSize',5);
plot(timesmooth(startind-1 + sigind{2}),repmat(0.47,...
    [1, length(sigind{2})]),'o','Color',[105 170 173]/256,'MarkerSize',5);
plot(timesmooth(startind-1 + sigind{3}),repmat(0.49,...
    [1, length(sigind{3})]),'o','Color',[0 0 0],'MarkerSize',5);

xticks([-0.5:0.25:2]);
box off;
xlim([-0.25 1.5]);
legend('Blind','Sighted');
legend('boxoff');
title(['Mean of binary comparisons of ',num2str(param.n_category),' categories'],'FontSize',12);
xlabel('time','FontSize',12);
ylabel('auc','FontSize',12);
set(gca, 'FontName', 'Avenir');
set(gca,'linewidth',2)
