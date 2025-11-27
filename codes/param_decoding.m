% Written by Siddharth Talwar
% Last edited on 26-11-2025

% Decode timepoint by timepoint using svm
% Each section below contains an example of running different analyses
% 1) to decode sounds considering all sensors as features
% 2) Sensor searchlight - to decode sounds at each seed channel and
%                         neighbouring sensors.
% 3) Temporal Generalization
% 4) Source searchlilght - to decode sounds at each source point within 20
%                          ms.

% Inputs/param-
% To be used after preprocessing. Data needs to be in preproc folder in
% fieldtrip format.
%
% Inputs -
% 1) param.n_category   : number of categories
% 2) param.timewindow   : [begin end] time in data, [] for all timepoints
% 3) param.fs           : sampling frequency of the dataset
% 4) param.ds           : downsample time by a factor (2,4,8 etc)
% 5) param.n_iterations : number of iterations (from creating pseudotrials
%                         to decoding.
% 6) param.pseudo_k     : number of trials to be averaged to make pseudotr
% 7) param.pseudo_n     : number of pseudotrials reqd for decoding
%                         Put EITHER pseudo_k OR pseudo_n.
% 8) param.popn         : population name / folder name
%
% NOTES
% - Each function takes the preprocessed file as input stored in preprocess 
% - All outputs are saved in the derivatives folder, in the corresponding
% analysis folder and sub-folder for each subject
% - The function uses parallel toolbox. Make sure there are enough workers.
% - Pseudotrials are created using an inbuilt function make_pseudo.
% - Cross validation folds are implemented using leave one trial out.
% - CHANGE paths in main function
% Provide pseudo_k (num trials to avg) or pseudo_n (num pseudotrials to
% calculate. Leave the other as '[]'.

clear
close all
clc

restoredefaultpath
parent  = '/media/siddharth/DATA/CPP/Projects/Aud_Cat/codes/Temporal-dynamics-of-natural-sound-representations-in-the-brain-of-sighted-and-blind/codes'; %change
addpath(parent); %change

%% Parameters

param                = [];
param.n_category     = 24; % 2, 4, 8 (category) or 24 (pairwise)
param.timewindow     = [-0.1 1.5]; %[begin end] in seconds
param.fs             = 256; % sampling frequency of preprocessed data
param.ds             = 2; % downsample factor
param.n_iterations   = 10; % number of iterations
param.pseudo_k       = []; % number of trials to avg OR
param.pseudo_n       = 3; % OR provide num of pseudotr needed
param.popn           = 'Blind'; % population/folder name
param.sub            = 1:18; % 1:n or specific subjects eg [2 3]

%% Example for pairwise or categorical timepoint by timepoint decoding considering all channels
% Input - same as above
% Outputs-
% 1) acc_comp_iter   : computes accuracy of classification between each
%                      comparison. comparison x time x iter
% 2) auc_comp_iter   : computes AUC of classification between each
%                      comparison. comparison x time x iter
% 3) comparison:     : indexed list of all binary comparisons
% 4) param           : input parameters

fn_svm_decode_libsvm(param);

%% Example for sensor searchlight pairwise or categorical decoding
% Input - same as above
% Outputs-
% 1) acc_comp_iter   : computes accuracy of classification between each
%                      comparison at each seed channel.
%                      channel X comparison x time x iter
% 2) auc_comp_iter   : computes AUC of classification between each
%                      comparison at each seed channel.
%                      channel X comparison x time x iter
% 3) comparison      : indexed list of all binary comparisons
% 4) param           : input parameters

fn_svm_decode_libsvm_chan(param)

%% Example for temporal generalization
% Input - same as above
% Outputs-
% 1) acc_comp_iter   : computes accuracy of classification between each
%                      comparison at each seed channel.
%                      channel X comparison x time x iter
% 2) auc_comp_iter   : computes AUC of classification between each
%                      comparison at each seed channel.
%                      channel X comparison x time x iter
% 3) comparison      : indexed list of all binary comparisons
% 4) param           : input parameters

fn_svm_decode_tempgen(param)

%% Example for source searchlight pairwise or categorical decoding

% Outputs
% 1) acc_mean           : Classfication accuracy for each sourcepoint and
%                         time. num_src x time, averaged iterations.
% 2) auc_mean           : Classfication AUC for each sourcepoint and
%                         time. num_src x time, averaged iterations.
% 3) param              : input parameters

% Input
param                = [];
param.n_category     = 2;       % 2, 4, 8 or 24 input
param.timewindow     = [0 0.2];   %[begin end] in seconds
param.n_iterations   = 1;     % number of iterations
param.pseudo_k       = [];      % number of trials to avg OR
param.pseudo_n       = 3;       % OR provide num of pseudotr needed
param.popn           = 'Blind'; % population/folder name
param.sub            = [1:18];       % 1:n or specific subjects eg [2 3]
fn_computesourcetime_search(param)
