%% Preprocess EEG
% Written by Siddharth Talwar
% Last edited on 26-11-2025
% Remarks: Check for "%change" and Change addresses only in Initialize 
% section. 
% Sequence of process: import, filtering, epoching, baseline subtraction,
% downsampling, semi-automatic artifact rejection, rereferencing
% Run separatey for each population (popn).

%% Initialize

clear
close all
clc

popn        = 'Control'; %change or blinds or others

% add fieldtrip to path
restoredefaultpath % comment this line if you have saved other paths
addpath '/mnt/DATA/siddharth/toolbox/fieldtrip-20191024' %change
ft_defaults

% add projectname/parent folder and subfolders to path
parent      = '/mnt/DATA/siddharth/Aud_Cat'; %change
cd(parent);
raw         = fullfile(parent,'raw',popn,'eeg');
preproc     = fullfile(parent,'preprocess');

dir_raw     = dir([raw,'/*.bdf']);
names       = {dir_raw([dir_raw.isdir]).name};
dir_raw(ismember(names,{'.','..'}))=[];

dir_preproc = dir(preproc);
names       = {dir_preproc([dir_preproc.isdir]).name};
dir_preproc(ismember(names,{'.','..'}))=[];

addpath '/mnt/DATA/siddharth/Aud_Cat/codes/eegmvpa'; %change

%% Preprocessing Loop

for sub = 1:length(dir_raw)
    
    % mkdir in destination preprocess folder
    mkdir(fullfile(preproc,popn,['sub',sprintf('%03d',sub)])); %change
    
    % find dir of eeg files
    rawfile         = fullfile(raw,['Sub',sprintf('%03d',sub),'.bdf']);
    
    % import events and data, choose channel info
    cfg             = [];
    cfg.datafile    = rawfile;
    cfg.headerfile  = rawfile;
    cfg.eventfile   = rawfile;
    
    %%%%%%% Conditions specific to the recorded dataset and recording notes
    
    cfg.channel     = [1:69 130 71:127 129]; % from A1 to EXG2 for D32 and C6 %change DEFAULT
    
    if (sub == 8 && strcmp(popn,'Blind') == 1) % for one exceptional blind subject etra external b2
        
        cfg.channel = [1:33 133 35:63 134 65:69 130 71:127 129]; %129-D32, 130-C6, 131-M1, 132-M2, 133-B2, 134-B32
        
    elseif (sub == 3 && strcmp(popn,'Control') == 1) || (sub == 9 && strcmp(popn,'Blind') == 1) || (sub == 11 && strcmp(popn,'Blind') == 1)
        
        cfg.channel = [1:63 133 65:69 130 71:127 129]; %129-D32, 130-C6, 131-M1, 132-M2, 133-B32
        
    elseif (sub == 10 && strcmp(popn,'Blind') == 1 || sub >= 4 && sub <= 10 && strcmp(popn,'Control') == 1) || (sub >= 12 && strcmp(popn,'Blind') == 1 && sub < 16)  % case consistent after including 12th blind and 4th control
        
        cfg.channel = [1:63 133 65:69 130 71:97 134 99:127 129]; %129-D32, 130-C6, 131-M1, 132-M2, 133-B32, 134-D2
        
    elseif  sub == 11 && strcmp(popn,'Control') == 1  % Extra D1 used as external, 11th control only
        
        cfg.channel = [1:63 133 65:69 130 71:96 135 134 99:127 129]; %129-D32, 130-C6, 131-M1, 132-M2, 133-B32, 134-D2, 135-D1
        
    elseif (sub == 12 && strcmp(popn,'Control') == 1) || (sub == 16 && strcmp(popn,'Blind') == 1) % Extra B19 used as external
        
        cfg.channel = [1:50 135  52:63 133 65:69 130 71:97 134 99:127 129]; %129-D32, 130-C6, 131-M1, 132-M2, 133-B32, 134-D2, 135-B19
    
    elseif (sub >= 13 && strcmp(popn,'Control') == 1) || (sub >= 17 && strcmp(popn,'Blind') == 1)  % Extra B19 and D1 used as external
        
        cfg.channel = [1:50 135  52:63 133 65:69 130 71:96 136 134 99:127 129]; %129-D32, 130-C6, 131-M1, 132-M2, 133-B32, 134-D2, 135-B19, 136-D1   
        
    end

    % for one subject, A3 and A4 got switched
    if  sub == 15 && strcmp(popn,'Blind') == 1
        cfg.channel = [1:2 4 3 5:63 133 65:69 130 71:97 134 99:127 129]; %129-D32, 130-C6, 131-M1, 132-M2, 133-B32, 134-D2
    end
    
    data            = ft_preprocessing(cfg);
    events          = ft_read_event(cfg.eventfile);
    
    % switching channel info as per biosemi system in LLN
    data.label{70}  = 'C6'; % for all
    data.label{128} = 'D32'; % for all
    
    % for some exceptions because of unstable electrodes
    data.label{34} = 'B2';
    data.label{64} = 'B32';
    data.label{98} = 'D2';
    data.label{97} = 'D1';
    data.label{3}  = 'A3';
    data.label{4}  = 'A4';
    data.label{51} = 'B19';
    %%%%%%% End of Conditions specific to the recorded dataset and recording notes
    
    % filter- highpass at 0,1, lowpass at 200 and bandstop at 50, 100, 150
    cfg                     = [];
    cfg.lpfilter            = 'yes';
    cfg.lpfreq              = 100;
    cfg.lpfilttype          = 'but';
    cfg.lpfiltord           = 4;
    cfg.hpfilter            = 'yes';
    cfg.hpfreq              = 0.1;
    cfg.hpfilttype          = 'but';
    cfg.hpfiltord           = 4;
    cfg.bsfilter            = 'yes';
    cfg.bsfreq              = [49.8 50.2; 99.8 100.2];
    cfg.bsfilttype          = 'but';
    cfg.bsfiltord           = 4;
    fi                      = ft_preprocessing(cfg,data);
    
    % epoching [-0.5 2.5] w.r.t stimulus onset
    cfg                     = [];
    cfg.trialdef.prestim    = 0.50;  % in seconds
    cfg.trialdef.poststim   = 2;    % in seconds
    cfg.trialdef.eventtype  = 'STATUS';
    cfg.trialdef.eventvalue = 1:24; % trigger values in expt
    cfg.dataset             = rawfile;
    cfg.headerfile          = rawfile;
    cfg_event               = ft_definetrial(cfg);
    
    cfg                     = [];
    cfg.trl                 = cfg_event.trl;
    tr_fi                   = ft_redefinetrial(cfg,fi);
    
    % baseline subtraction w.r.t. prestim
    cfg                     = [];
    cfg.demean              = 'yes';
    cfg.baselinewindow      = [-inf 0];
    bl_tr_fi                = ft_preprocessing(cfg,tr_fi);
    
    %resampling/downsampling
    cfg                     = [];
    cfg.resamplefs          = 256;
    ds_bl_tr_fi             = ft_resampledata(cfg,bl_tr_fi);
    
    % detrend for particular subjects because of unstable acquisition
    % reference
    if sub == 15 && strcmp(popn,'Blind') == 1
        cfg = [];
        cfg.detrend = 'yes';
        ds_bl_tr_fi  = ft_preprocessing(cfg,ds_bl_tr_fi);
    end
    
    % semi-automatic artifact rejection
    % trials > 250 uV : automatically deleted
    % trials > 150 and < 250 : plots. press 'd' to discard,any other key to
    % keep the trial
    % low pass to smoothen the spikes (ONLY FOR AR)
    ii=1;
    cfg = [];
    cfg.lpfilter            = 'yes';
    cfg.lpfreq              = 30;
    cfg.lpfilttype          = 'but';
    cfg.lpfiltord           = 4;
    temp                    = ft_preprocessing(cfg,ds_bl_tr_fi);
    
    artfct_tr = [];
    for tr = 1:length(temp.trial)
        
        if any( any( abs( temp.trial{tr} ) > 250 ) ) %change automatic delete threshold
            
            artfct_tr(ii) = tr;
            ii            = ii+1;
            continue
            
        end
        if any( any( abs( temp.trial{tr} ) > 150 ) ) %change plot threshold
            
            figure;
            plot(temp.trial{tr}');
            ylim([-200 200]); %y axis scale
            title([num2str(tr),' ',num2str(temp.trialinfo(tr))]);
            inp = input('Reject - keep or discard? ','s');
            
            if strcmp(inp,'d')
                
                artfct_tr(ii) = tr;
                ii            = ii+1;
                
            end
            close all
        end
    end
    
    % delete noisy trials
    alltr                 = 1 : length(ds_bl_tr_fi.trial);
    alltr(artfct_tr)      = [];
    
    cfg                   = [];
    cfg.trials            = alltr;
    ar_ds_bl_tr_fi        = ft_preprocessing(cfg,ds_bl_tr_fi);
    
    % rereferencing
    cfg                   = [];
    cfg.reref             = 'yes';
    cfg.refchannel        = 'all';
    re_ar_ds_bl_tr_fi     = ft_preprocessing(cfg,ar_ds_bl_tr_fi);
    
    preprocdata           = re_ar_ds_bl_tr_fi;
    
    % visualize data
    tl                    = ft_timelockanalysis([],preprocdata);
    
    figure;
    plot(tl.avg');
    
    % find minimum number of trials per sound
    [sort_tr,~,n_tr] = make_pseudo(preprocdata,1,[]);
    for ii = 1:length(sort_tr)
        minimum(ii) = min(size(sort_tr{ii},3));
    end
    
    minimum = min(minimum);
    disp(['Min trials for a sound are ', num2str(minimum),'. Press any key to save.']);
    pause;
    
    save(fullfile(preproc,popn,['sub',sprintf('%03d',sub)],'preproc.mat'),'preprocdata');
    clear data fi tr_fi bl_tr_fi ds_bl_tr_fi ar_ds_bl_tr_fi re_ar_ds_bl_tr_fi preprocdata trials data
    
end