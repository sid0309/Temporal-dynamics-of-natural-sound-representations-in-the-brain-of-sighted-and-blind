function fn_computesourcetime_search(param)

% Source Time Series Decoding
% Written by Siddharth Talwar
% Last edited on 26-11-2025
% Inputs/param
% 1) param.popn         : population or folder name
% 2) param.sub          : subject indices
% 3) param.n_category   : number of categories
% 4) param.n_iterations : number of iterations
% 5) param.pseudo_k     : number of trials to be averaged to make pseudotr
% 6) param.pseudo_n     : number of pseudotrials reqd for decoding
%                         Put EITHER pseudo_k OR pseudo_n.
% 7) param.timewindow   : [begin end] time in data, [] for all timepoints
%
% Outputs
% 1) acc_mean           : Classfication accuracy for each sourcepoint and
%                         time. num_src x time, averaged iterations.
% 2) auc_mean           : Classfication AUC for each sourcepoint and
%                         time. num_src x time, averaged iterations.
% 3) param              : input parameters
%
% NOTES
% - This function computes participant specific BEM sourcemodel and
% leadfield matrix using Fieldtrip. A template is used for participants who
% don't have their headmodels or mri. Fieldtrip uses OpenMNE for
% computation, so make sure to download OpenMNE and launch MATLAB in
% the same environment (python/conda).
% - The spatial filter of LCMV is calculated across all trials and length.
% param.timewindow extracts the timewindow of interest after the
% computation of the spatial filter. Consequenty, data within 20 ms is
% averaged.
% - Calls another function fn_svm_decode_libsvm_src to decode the resuting
% source timeseries.
% - All outputs are saved in the derivatives folder, in the corresponding
%   analysis folder and sub-folder for each subject
% - The function uses parallel toolbox. Make sure there are enough workers.
%%
% Paths (change)

% change
parent  = 'xxx/Temporal-dynamics-of-natural-sound-representations-in-the-brain-of-sighted-and-blind'; %change
addpath 'xxx/fieldtrip-20191024' %fieldtrip change
ft_defaults
addpath('xxx/libsvm-master/matlab') %libsvm change

% don't change
% cd(parent);
preproc = fullfile(parent,'preprocess');
deriv   = fullfile(parent,'derivatives');
dir_preproc = dir(fullfile(preproc,param.popn));
names       = {dir_preproc([dir_preproc.isdir]).name};
dir_preproc(ismember(names,{'.','..'})) = [];

%% Make template grid (comment this section if you can load a template)

template = ft_read_mri(fullfile(parent,'others','src','T1.nii'));
template.coordsys = 'spm'; % so that FieldTrip knows how to interpret the coordinate system

% segment the template brain and construct a volume conduction model (i.e. head model):
% this is needed to describe the boundary that define which dipole locations are 'inside' the brain.
cfg          = [];
template_seg = ft_volumesegment(cfg, template);

cfg          = [];
cfg.method   = 'openmeeg';
template_headmodel = ft_prepare_headmodel(cfg, template_seg);
template_headmodel = ft_convert_units(template_headmodel, 'cm'); % Convert the vol to cm, because the CTF convenction is to express everything in cm.

% construct the dipole grid in the template brain coordinates
% the negative inwardshift means an outward shift of the brain surface for inside/outside detection

cfg              = [];
cfg.method       = 'basedonmri';
cfg.threshold    = 0.1;
cfg.resolution   = 0.1;
cfg.smooth       = 5;
cfg.tight        = 'yes';
% cfg.inwardshift  = 0.5;
% cfg.moveinward   = 0.5;
cfg.mri          = template;
cfg.headmodel    = template_headmodel;
template_grid    = ft_prepare_sourcemodel(cfg);

% make a figure with the template head model and dipole grid
% figure
% hold on
% ft_plot_headmodel(template_headmodel, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.1; camlight;
% ft_plot_mesh(template_grid.pos(template_grid.inside,:));

%% Loop for subjects

% stuff to load
load(fullfile(parent,'others','src','templategrid_gray_10mm.mat')); % output file from above section
tmplate = load(fullfile(parent,'others','src','fwd_model_search.mat'));
neigh = 5; % searchlight neighbours

for sub = param.sub
    clc
    
    % create the subject specific grid using the template grid
    if exist(fullfile(deriv,param.popn,['sub',sprintf('%03d',sub)],'fwdmodel')) == 7
        fwdmodel=load(fullfile(deriv,param.popn,['sub',sprintf('%03d',sub)],'fwdmodel','fwdmodel_meter.mat'));
        
        cfg              = [];
        cfg.warpmni      = 'yes';
        cfg.template     = template_grid;
        cfg.nonlinear    = 'yes';
        cfg.mri          = ft_convert_units(fwdmodel.mri_resliced,'cm');
        cfg.unit         ='cm';
        cfg.inwardshift  = 0.5;
        cfg.moveinward   = 0.5;
        cfg.headmodel    = ft_convert_units(fwdmodel.headmodel,'cm');
        fwdmodel.Sourcemodel = ft_prepare_sourcemodel(cfg);
        
        fwdmodel.Sourcemodel = ft_convert_units(fwdmodel.Sourcemodel,'m');
        
        % leadfield to be computed in meters
        cfg = [];
        cfg.sourcemodel = ft_convert_units(fwdmodel.Sourcemodel,'m');    %% where are the sources?
        cfg.headmodel   = ft_convert_units(fwdmodel.headmodel,'m');      %% how do currents spread?
        cfg.elec        = ft_convert_units(fwdmodel.elec_realigned,'m'); %% where are the sensors?
        fwdmodel.Leadfield = ft_prepare_leadfield(cfg);
        
        indiv = 1;
        
    else
        % create the subject specific grid, using the template grid that has just been created
        fwdmodel = tmplate;
        
        mri = ft_read_mri('/media/siddharth/DATA/Toolbox/fieldtrip-20191024/external/spm8/templates/T1.nii');
        mri.coordsys = 'acpc'; % checked before
        
        
        cfg              = [];
        cfg.warpmni      = 'yes';
        cfg.template     = template_grid;
        cfg.nonlinear    = 'yes';
        cfg.mri          = ft_convert_units(mri,'cm');
        cfg.unit         ='cm';
        cfg.inwardshift  = 0.5;
        cfg.moveinward   = 0.5;
        cfg.headmodel    = ft_convert_units(fwdmodel.headmodel,'cm');
        fwdmodel.Sourcemodel = ft_prepare_sourcemodel(cfg);
        
        fwdmodel.Sourcemodel = ft_convert_units(fwdmodel.Sourcemodel,'m');
        
        cfg = [];
        cfg.sourcemodel = ft_convert_units(fwdmodel.Sourcemodel,'m');    %% where are the sources?
        cfg.headmodel   = ft_convert_units(fwdmodel.headmodel,'m');      %% how do currents spread?
        cfg.elec        = ft_convert_units(fwdmodel.elec_realigned,'m'); %% where are the sensors?
        fwdmodel.Leadfield       = ft_prepare_leadfield(cfg);
        indiv = 0;
    end
    
    % import data
    load(fullfile(parent,'preprocess',param.popn,dir_preproc(sub).name,...
        'preproc.mat'),'preprocdata');
    
    % num of pseudotrials to name a folder
    if isempty(param.pseudo_n) == 1
        [~,~,n_tr] = make_pseudo(preprocdata,param.pseudo_k,...
            param.pseudo_n);
    else
        n_tr = param.pseudo_n;
    end
    
    % make output folders
    add  = fullfile(deriv,param.popn,['sub',sprintf('%03d',sub)],...
        ['binarydecoding_svm_',num2str(param.n_category),'cat']);
    mkdir(add);
    dest = fullfile(add,['iter',num2str(param.n_iterations)],'sourcetime');
    mkdir(dest);
    
    % timelockanalysis
    cfg            = [];
    cfg.covariance = 'yes';
    TL             = ft_timelockanalysis(cfg,preprocdata);
    
    % compute spatial filter using all the data
    cfg                    = [];
    cfg.channel            = fwdmodel.elec_realigned.label;
    cfg.elec               = fwdmodel.elec_realigned;
    cfg.method             = 'lcmv';
    cfg.grid               = fwdmodel.Leadfield;
    cfg.headmodel          = fwdmodel.headmodel;
    cfg.lcmv.keepfilter    = 'yes';
    cfg.sourcemodel        = fwdmodel.Sourcemodel;
    cfg.keepfilter         = 'yes';
    cfg.lcmv.lambda        = '5%';
    %         cfg.lcmv.fixedori      = 'yes';
    %         cfg.reducerank         = 3;
    S                      = ft_sourceanalysis(cfg, TL);
    
    % find eeg channels
    chansel = ft_channelselection('EEG', preprocdata.label);   % find the names
    chansel = match_str(preprocdata.label, chansel);           % find the indices
    
    % extract filter and pos of source points inside the brain
    bf_filt = S.avg.filter(S.inside==1);
    bf_pos  = S.pos(S.inside==1,:);
    
    sph = find(fwdmodel.Sourcemodel.inside==1);
    
    % extract timewindow of interest
    if isempty(param.timewindow) ~= 1
        cfg              = [];
        cfg.toilim       = param.timewindow;
        preprocdata = ft_redefinetrial(cfg,preprocdata);
    end
    
    % compute number of comparisons
    comparison = flip(combnk(1:param.n_category,2));
    
    if comparison(1,1)~=1 % make sure the combinations start with the 1st
        comparison = flip(comparison);
    end
    
    time = preprocdata.time{1};
    
    timewin = [param.timewindow(1):0.02:param.timewindow(2)];
    
    % average window 20 ms
    trial = cell(1,length(preprocdata.trial));
    for ii = 1:length(preprocdata.trial)
        kk=1;
        for jj = 1:length(timewin)-1
            
            tmp = find(preprocdata.time{1} > timewin(jj) & preprocdata.time{1} < timewin(jj+1));
            trial{ii}(:,kk) = mean(preprocdata.trial{ii}(:,tmp),2);
            kk=kk+1;
        end
    end
    
    preprocdata.trial = trial;
    trial=[];
    
    acc_mean =...
        zeros(length(sph),size(comparison,1),...
        size(preprocdata.trial{1},2));
    auc_mean = acc_mean;
    
    disp(['Computing source-time activity from each cluster for sub ',...
        num2str(sub),' and decoding ']);
    lineLength = 1;
    
    sm = fwdmodel.Sourcemodel.pos;
    
    % Compute sourcetime series and decoding
    parfor roi = 1:length(sph)
        
        for aa = 1:length(sph)
            distance(aa) = pdist([sm(sph(roi,:));sm(sph(aa,:))]);
        end
        
        [distance,srcind] = sort(distance);
        current_sph = srcind(1:neigh); % by number of points
        %         current_sph = srcind(distance<(1.5/1000)); % by distance
        
        virtualchanneldata= cell(1,length(current_sph));
        
        % create time series from the each point in sphere
        for source = 1:length(current_sph)
            
            virtualchanneldata{source} = [];
            virtualchanneldata{source}.label = {'cortex'};
            virtualchanneldata{source}.time = preprocdata.time;
            
            % compute time series for sphere
            for i=1:length(preprocdata.trial) % loop for trial
                
                virtualchanneldata{source}.trial{i} =...
                    bf_filt{current_sph(source)} * preprocdata.trial{i}(chansel,:);
                %                 [u, s, v] = svd(virtualchanneldata{source}.trial{i}, 'econ');
                %                 virtualchanneldata{source}.trial{i} = (u(:,1)'* virtualchanneldata{source}.trial{i}); % divide by 1000 to scale down
            end
            
            timeseries = cat(2, virtualchanneldata{source}.trial{:});
            
            [u, s, ~] = svd(timeseries, 'econ');
            
            for i=1:length(preprocdata.trial) % loop for trial
                virtualchanneldata{source}.trial{i} =...
                    (u(:,1)'*bf_filt{current_sph(source)}*  preprocdata.trial{i}(chansel,:));
            end
        end
        
        % put all source time series in one structure
        roi_timeseries= virtualchanneldata{1};
        
        for tr = 1:length(virtualchanneldata{1}.trial)
            roi_timeseries.trial{tr} = [];
            for source = 1:length(virtualchanneldata)
                roi_timeseries.trial{tr} =...
                    cat(1,roi_timeseries.trial{tr},virtualchanneldata{source}.trial{tr});
            end
        end
        roi_timeseries.trialinfo = preprocdata.trialinfo; %trial trigger structure intact
        roi_timeseries.num_s     = neigh;
        roi_timeseries.label     = repmat({'cortex'},[roi_timeseries.num_s,1]);
        
        [acc_mean(roi,:,:,:),auc_mean(roi,:,:,:)] = ...
            fn_svm_decode_libsvm_src(param,roi_timeseries,comparison); % roi x comp x time
        
        %         fprintf(repmat('\b',1,lineLength));
        %         lineLength = fprintf('%2.1f percent',...
        %             roi./length(sph)*100);
        
    end
    
    if indiv == 1
        save(fullfile(dest,...
            ['auc_sourcetime_svd_indivfwdmodel_search_','acc_libsvm_src_indivmodel.mat'),'acc_mean','auc_mean','param');
    elseif indiv == 0
        save(fullfile(dest,'acc_libsvm_src_tmplatemodel.mat'),'acc_mean','auc_mean','param');
    end
end
