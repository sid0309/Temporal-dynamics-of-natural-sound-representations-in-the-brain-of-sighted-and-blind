%% Plotting
clear
clc
parent  = '/media/siddharth/DATA/CPP/Projects/Aud_Cat/codes/Temporal-dynamics-of-natural-sound-representations-in-the-brain-of-sighted-and-blind'; %change path

deriv   = fullfile(parent,'derivatives');
addpath '/media/siddharth/DATA/Toolbox/fieldtrip-20191024'
ft_defaults

%input
param.n_category   = 24;% category folder name
param.iter         = 100;% num iterations
param.factor       = 0; % smoothing factor
param.pseudo_n     = 3; % number of pseudotrials
param.trial_radius = 0; % number of neighbourhood timepoints
param.popn         = {'Blind','Control'}; % {popn1,popn2}
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
            '_pt',num2str(param.pseudo_n),...
            '_tr',num2str(param.trial_radius)],'acc_chan_libsvm.mat'),...
            'auc_comp_iter_chan','comparison');
        
        auc_chan_comp_tp = squeeze(mean(auc_chan_comp_tp,2)); % mean of comparisons
        acc{pop}(:,:,sub_ii) = auc_chan_comp_tp; % comp x time x sub
        
        sub_ii = sub_ii+1;
        
    end
end

timesl = -0.25 + (1.75/size(acc{1},2)):1.75/size(acc{1},2):1.5; % always 640
% timesmoothsl = time(1+param.factor:length(time)-param.factor);

%% layout

%load dummy var
load('/media/siddharth/DATA/CPP/Projects/Aud_Cat/preprocess/Pilot/sub001/preproc.mat');
cfg = [];
cfg.layout = '/media/siddharth/DATA/Toolbox/fieldtrip-20191024/template/layout/biosemi128.lay';
lay = ft_prepare_layout(cfg);

%% stats

load('/media/siddharth/DATA/CPP/Projects/Aud_Cat/Manuscript/plots_apr24/fig_generaldecode/newstats/clustersindx_gendecode.mat')
tp = 56:68;%33:45; %52:69 70:84 %87:101 %104:121 125:137

vec{1} = squeeze(mean(acc{1}(:,tp,:),2))-0.5; % only poststim
vec{2} = squeeze(mean(acc{2}(:,tp,:),2))-0.5;

% %arrange sc as per eb
% scindx = [10 2 3 4 11 15 13 7 17 5 9 1 6 8 12 18 14 16];
% vec{2} = vec{2}(:,scindx);
% 
% % paired subtraction
% vec{3} = vec{1} - vec{2};


for popn = 1:2
    sigind{popn} = [];
    %3rd is the difference b/w 2 popn
    
    % sign perm1-real(dsm)utation test
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
    %     cfg.colormap  = 'jet';
    %     subplot(2,1,popn);
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
% colormap((colormap('turbo')));

%% plot
%1:7, 8:14, 15:21, 21:27, 28:34, 35:40, 41:47 , 48:53, 54:60
aa=1;
figure;
tp = [{1:7}, {8:14}, {15:21}, {21:27}, {28:34}, {35:40}, {41:47} , {48:53}, {54:60}];
tp = [{2:10}, {12:36}, {40:57}];
tp = [{1:11}, {12:24}, {25:37}, {38:60}];
zax = [0.604 0.5759 0.5761 0.5684];

% sigind = [];
for popn = 1:3
    
    for tpind = 1:4
        %160-293 ms (2:10), 308-695 ms (12:36), 730-1007 ms (40:57)
        
        vec{1} = squeeze(mean(acc{1}(:,tp{tpind},:),2))-0.5; % only poststim
        vec{2} = squeeze(mean(acc{2}(:,tp{tpind},:),2))-0.5;
        
        %arrange sc as per eb
        scindx = [10 2 3 4 11 15 13 7 17 5 9 1 6 8 12 18 14 16];
        vec{2} = vec{2}(:,scindx);
        
        % paired subtraction
        vec{3} = vec{1} - vec{2};
        %3rd is the difference b/w 2 popn
        
        if popn<3
        % sign perm1-real(dsm)utation test
        disp('Sign permutation test');
            [pval,nulldist,sigind] =  fn_signperm_chan(vec{popn} ,5000,0.01,lay);
        
        
        topo{popn}             = preprocdata;
        topo{popn}             = rmfield(topo{popn},'trialinfo');
        topo{popn}             = rmfield(topo{popn},'sampleinfo');
        topo{popn}.trial       = {};
        topo{popn}.trial{1}    = squeeze(mean(mean(acc{popn}(:,tp{tpind},:),2),3));
        topo{popn}.time(2:end) = [];
        topo{popn}.time{1}     = timesmooth(tp{tpind}(1));
        tl               = ft_timelockanalysis([],topo{popn});
        
        cfg              = [];
        cfg.parameter    = 'avg';
        cfg.layout       = lay;
        cfg.comment      = 'no';
        cfg.marker       = 'off';
        cfg.highlight    = 'on';
        cfg.highlightsymbol = '.';
        cfg.highlightsize   = 8;
        cfg.highlightchannel = sigind;
        cfg.zlim      = [0.5 zax(tpind)]; %0.5736  0.6094 
        %     cfg.colormap  = 'jet';
        %     subplot(2,1,popn);
        subplot(3,4,aa);
        ft_topoplotER(cfg, tl);
        colorbar;
        title([num2str(floor(timesmooth(tp{tpind}(1))*1000)),' to ',...
            num2str(floor(timesmooth(tp{tpind}(end))*1000)),' ms'],'FontSize',8);
        set(gca, 'FontName', 'Arial');
%         colormap(flipud(colormap('hot')));
        aa=aa+1;
        
        elseif popn==3
            % blind - sighted
            topo{popn} = topo{1};
            topo{popn}.trial{1} = squeeze(mean(mean(acc{1}(:,tp{tpind},:),2),3)) -...
                squeeze(mean(mean(acc{2}(:,tp{tpind},:),2),3));
            
            zaxmin = squeeze(mean(mean(acc{2}(:,tp{tpind},:),2),3)) -...
                squeeze(mean(mean(acc{1}(:,tp{tpind},:),2),3));
            
            tl = ft_timelockanalysis([],topo{popn});
            
            disp('Sign permutation test');
                [pval,nulldist,sigind] =  fn_signperm_chan(vec{3} ,5000,0.05,lay);
            
            
            cfg.highlightchannel = sigind;%
            cfg.zlim      = [0 max(topo{popn}.trial{1})];% 0.0513
            subplot(3,4,aa);
            ft_topoplotER(cfg, tl);
            colorbar;
            title([num2str(floor(timesmooth(tp{tpind}(1))*1000)),' to ',...
                num2str(floor(timesmooth(tp{tpind}(end))*1000)),' ms'],'FontSize',8);
            set(gca, 'FontName', 'Arial');
%             colormap(flipud(colormap('hot')));
            aa=aa+1;
        end
    end
end

%% blind - sighted

topo.trial{1} = squeeze(mean(pat{1}(:,tp,:),3)) -...
    squeeze(mean(pat{2}(:,tp,:),3));
tl               = ft_timelockanalysis([],topo);
figure;
ft_topoplotER(cfg, tl);
colorbar;
sz([num2str(floor(timesmooth(tp(1))*1000)),' to ',...
    num2str(floor(timesmooth(tp(end))*1000)),' ms'],'FontSize',20);
set(gca, 'FontName', 'Arial');

%% sensor search topoplots movie, DSM and MDS
param.factor=2;
time = -0.5 + (2.5/640):2.5/640:2; % always 640
timesmooth = time(1+param.factor:length(time)-param.factor);

time_sl = -0.25 + (1.75/225):1.75/225:1.5; % always 640

% Change address 3 places below
%%%%%%%%%%%%%%%%%%%%%%%
% load eeg dsm (change address)
eeg = load(fullfile('/media/siddharth/DATA/CPP/Projects/Aud_Cat/others'...
    ,'eegdsm/eegdsm_libsvm.mat'));

% load fieldtrip dummy var (fieldtrip address, also for layout)
load('/media/siddharth/DATA/CPP/Projects/Aud_Cat/preprocess/Pilot/sub001/preproc.mat');
cfg = [];
cfg.layout = '/media/siddharth/DATA/Toolbox/fieldtrip-20191024/template/layout/biosemi128.lay';
lay = ft_prepare_layout(cfg);

% load sources
% load('/media/siddharth/DATA/CPP/Projects/Aud_Cat/others/srcloc/sourceprojnew_mne.mat')
%%%%%%%%%%%%%%%%%%%%%%%%%

% when to when
% startind = find(timesmooth>0);
startind = 126;
tprange = 126:509;

color = {[219 44 19]/256, [233 163 34]/256, [56 74 195]/256, [100 178 198]/256, ...
    [188 44 125]/256, [225 108 170]/256,...
    [89 157 112]/256,[137 242 150]/256, };

% for colorbar limit
maxim = max(max(mean(acc{1},3)));
minim = min(min(mean(acc{1},3)));
diffmaxim = max(max(mean(acc{1},3) - mean(acc{2},3)));
diffminim = min(min(mean(acc{1},3) - mean(acc{2},3)));
popname = {'eb','sc'};

% calculate mds before to fix axis
Y1 = []; Y2 = [];
for tp = tprange
    dsm1       =  mean(mean(eeg.avg_eegdsm{1}(:,:,[tp-1:tp+1],:),3),4);
    dsm2       =  mean(mean(eeg.avg_eegdsm{2}(:,:,[tp-1:tp+1],:),3),4);
    dsm3       = (dsm1 - dsm2);
    dsm3(find(dsm3<0)) = 0;
    Y1(:,:,tp) = cmdscale(dsm1,2);
    Y2(:,:,tp) = cmdscale(dsm2,2);
    Y3(:,:,tp) = cmdscale(dsm3,2);
end

xmax1 = max(max(Y1(:,1,:)));xmax2 = max(max(Y2(:,1,:)));xmax3 = max(max(Y3(:,1,:)));
xmin1 = min(min(Y1(:,1,:)));xmin2 = min(min(Y2(:,1,:)));xmin3 = min(min(Y3(:,1,:)));
ymax1 = max(max(Y1(:,2,:)));ymax2 = max(max(Y2(:,2,:)));ymax3 = max(max(Y3(:,2,:)));
ymin1 = min(min(Y1(:,2,:)));ymin2 = min(min(Y2(:,2,:)));ymin3 = min(min(Y3(:,2,:)));

%% Prepare the new file for video
vidObj = VideoWriter('topo_libsvm_sl.avi');
vidObj.FrameRate = 4;
open(vidObj);

for tp = 126:509
   
    tmp = time_sl-time(tp);
    sl_tp = find(tmp>0);
    sl_tp = (sl_tp(1));
    
    % timelock analysis at each tp
    topo             = preprocdata;
    topo             = rmfield(topo,'trialinfo');
    topo             = rmfield(topo,'sampleinfo');
    topo.trial       = {};
    topo.time(2:end) = [];
    topo.time{1}     = topo.time{1}(:,sl_tp);
    topo.trial{1}    = squeeze(mean(mean(acc{1}(:,sl_tp,:),3),2));
    tl1              = ft_timelockanalysis([],topo);
    topo.trial{1}    = squeeze(mean(mean(acc{2}(:,sl_tp,:),3),2));
    tl2              = ft_timelockanalysis([],topo);
    topo.trial{1}    = squeeze(mean(mean(acc{1}(:,sl_tp,:),3),2)) - squeeze(mean(mean(acc{2}(:,sl_tp,:),3),2));
    tl3              = ft_timelockanalysis([],topo);
    
    % topoplot
    cfg              = [];
    cfg.parameter    = 'avg';
    cfg.layout       = lay;
    cfg.comment      = 'no';
    cfg.zlim         = [minim maxim];
    cfg.marker       = 'off';
    cfg.highlight    = 'off';
    cfg.colormap     = 'jet';
    figure('units','normalized','outerposition',[0 0 1 1]);
    
    im1 = subplot(3,3,1);
    ft_topoplotER(cfg, tl1);
    title('EB');
    
    %1st colorbar
    cb = colorbar('horiz');
    set(cb,'position',[.188 .645 .1 .015 ]);
    x = suptitle([num2str(floor(timesmooth(tp)*1000)),' ms']);
    set(x,'Position',[0.518 -0.05]);
    set(gca, 'FontName', 'Arial');
    
    cfg.zlim         = [minim maxim];
    im3 = subplot(3,3,3);
    ft_topoplotER(cfg, tl2);
    title('SC');
    
    %3rd colorbar
    cb = colorbar('horiz');
    set(cb,'position',[.748 .645 .1 .015 ]);
    x = suptitle([num2str(floor(timesmooth(tp)*1000)),' ms']);
    set(x,'Position',[0.518 -0.05]);
    set(gca, 'FontName', 'Arial');
    
       cfg.zlim         = [diffminim diffmaxim];
    im2 = subplot(3,3,2);
    ft_topoplotER(cfg, tl3);
    title('EB - SC');
    
    %2nd colorbar
    cb = colorbar('horiz');
    set(cb,'position',[.4677 .645 .1 .015 ]);
    x = suptitle([num2str(floor(timesmooth(tp)*1000)),' ms']);
    set(x,'Position',[0.518 -0.05]);
    set(gca, 'FontName', 'Arial');
    %%%
    
    subplot(3,3,4);
    ii=1;
    for category = 1:3:24
        scatter(Y1(category:category+2,1,tp),Y1(category:category+2,2,tp),150,...
            'MarkerFacecolor',color{ii},'MarkerEdgeColor','k');
        hold on;
        ii = ii+1;
    end
    xlim([xmin1 xmax1]);ylim([ymin1 ymax1]);
    box off;
    set(gca,'XColor', 'none','YColor','none');
    
    subplot(3,3,5);
    ii=1;
    for category = 1:3:24
        scatter(Y3(category:category+2,1,tp),Y3(category:category+2,2,tp),150,...
            'MarkerFacecolor',color{ii},'MarkerEdgeColor','k');
        hold on;
        ii = ii+1;
    end
    xlim([xmin3 xmax3]);ylim([ymin3 ymax3]);
    box off;
    set(gca,'XColor', 'none','YColor','none');
    
    subplot(3,3,6);
    ii=1;
    for category = 1:3:24
        scatter(Y2(category:category+2,1,tp),Y2(category:category+2,2,tp),...
            150,'MarkerFacecolor',color{ii},'MarkerEdgeColor','k');
        hold on;
        ii = ii+1;
    end
    xlim([xmin2 xmax2]);ylim([ymin2 ymax2]);
    box off;
    set(gca,'XColor', 'none','YColor','none');
    set(gcf,'color','w');
    
    subplot(3,3,7);
    ax = gca;
    dsmhandle = imagesc(mean(eeg.avg_eegdsm{1}(:,:,tp,:),4)/0.9427,'Parent', ax);
    ax.XAxis.Visible = 'off';
    ax.YAxis.Visible = 'off';
    ax.CLim = [0.5 1];
    cb = colorbar('vert');
    set(cb,'position',[.3 0.15 .01 .15  ]);
    axis equal
    axis tight
    
    subplot(3,3,8);
    ax = gca;
    dsmhandle = imagesc((mean(eeg.avg_eegdsm{1}(:,:,tp,:),4) -...
        mean(eeg.avg_eegdsm{2}(:,:,tp,:),4))/0.3063,'Parent', ax);
    ax.XAxis.Visible = 'off';
    ax.YAxis.Visible = 'off';
    ax.CLim = [0 1];
    cb = colorbar('vert');
    set(cb,'position',[.58 0.15 .01 .15  ]);
    axis equal
    axis tight
    
    subplot(3,3,9);
    ax = gca;
    dsmhandle = imagesc(mean(eeg.avg_eegdsm{2}(:,:,tp,:),4)/0.7871,'Parent', ax);
    ax.XAxis.Visible = 'off';
    ax.YAxis.Visible = 'off';
    ax.CLim = [0.5 1];
    cb = colorbar('vert');
    set(cb,'position',[.86 0.15 .01 .15  ]);
    axis equal
    axis tight
 
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
    close all;
end
% Close the file.
close(vidObj);

