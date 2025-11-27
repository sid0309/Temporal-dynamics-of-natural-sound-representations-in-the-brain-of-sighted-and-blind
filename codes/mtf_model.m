%% INTRODUCTION
% the following scripts provide the basic idea of the usage of [nsltools]
% you are welcomed to report any bugs, type-errors as well as suggestions to
% tschi@isr.umd.edu (Taishih Chi), sas@isr.umd.edu (Shihab Shamma)

%% SETUP
% make sure the nsltools is in the path list
% load colormap, filterbank and parameters
clear
clc
addpath '/media/siddharth/DATA/Toolbox/nsltools'
loadload;

% you may change the parameters, see WAV2AUD
% paras(1): frame jump, in ms, e.g., 4, 8, 16, 32 ms
% paras(2): Time constant, in ms, e.q., 4, 8, 16, 32, 64 ms
% paras(3): Nonlinear factor, e.g., .1
% paras(4): Octave shift, e.g., -1 for 8k, 0 for 16k (sampling rate)
% rv: rate vector, e.g., 2.^(1:5) or 2.^(1:.25:5); (Hz)
% sv: scale vector, e.g., 2.^(-2:3) or 2.^(-2:.25:3); (cycle/octave)
paras = [8 8 0 -1];
% WAVFORM

dir_sounds = dir('/media/siddharth/DATA/CPP/Projects/Aud_Cat/codes/Temporal-dynamics-of-natural-sound-representations-in-the-brain-of-sighted-and-blind/others/stimuli/*.wav');
%%
for sound = 1:length(dir_sounds)
    
    % sequential format:	see LOADFILE
    % .au format:		see AUREAD
    % .aiff format:		see AIFFREAD
    % x = loadfile('_done.au');
    [x] = audioread(fullfile(dir_sounds(sound).folder,dir_sounds(sound).name));
    
    % x = auread('_come.au');	% 8k sample rate
    
    % down sample to 8 kHz if necessary
    % if the sample rate is not the power of 2, you need to use RESAMPLE
    % change the paras(4) according to the sample rate
    % x = decimate(x, 5);
    %x = x(2:2:length(x));
    x = resample(x,8000,44100);
    
    % if you want to see a schematic plot, you may run
    % xh = schematc(x);
    
    % make the sequence as unit sequence, i.e., ~ N(0, 1) or
    % zero-mean and unit variance (optional but recommended)
    % besides, the length 8192, 16384, 32768, etc are preferrable
    x = unitseq(x);
    
    
    % SPECTRUM
    
    tp = 0:8000/(1):8000; % how many steps?????######
    tp(1) = 1;
    
    for ii = 1:length(tp)
        
        % y is an M-by-128 matrix, M is the # of frames
        % a short message will be displayed for every complete octave
        y{sound,ii} = wav2aud(x, paras);
        
        % plot the spectrogram
        %             aud_plot(y, paras);
        
        % CORTICAL REPRESENTATION
        % static representation (e.g., 60th frame of the spectrogram)
        % z = aud2cors(y, sv);
        % plot it
        % cor_plts(z, sv, paras);
        
        % dynamic 4-D representation
        % you can assign an output filename for saving the data see AUD2COR
        % The process will take a while. 15 seconds for 1 second input data on SPARC Ultra-5_10 machine with default parameters
        fcorname     = '_come.cor';
        cr{sound,ii} = aud2cor(y{sound,ii}, paras, rv, sv, fcorname);
%         crstatic{sound,ii}     = aud2cors(y{sound,ii}, sv);
        
%                 crnorm{sound,ii} = cr{sound,ii} ./ max(max(max(max(cr{sound,ii}))));
        %     (4D, scale-rate(up-down)-time-freq.)
    end
end

%% Make dsm

dsm=[];
for tp = 1:size(cr,2)
    
    for ii = 1:size(cr,1)
        
        for jj = 1:size(cr,1)
            
            tmp1 = mean(cr{ii,tp},3);
            tmp2 = mean(cr{jj,tp},3);
            
            z1 = reshape(tmp1,[size(tmp1,1)*size(tmp1,2)*...
                size(tmp1,3)*size(tmp1,4), 1]);
            
            z2 = reshape(tmp2,[size(tmp2,1)*size(tmp2,2)*...
                size(tmp2,3)*size(tmp2,4), 1]);
            
              dsm(ii,jj,tp) = pdist2(real(z1)',real(z2)','cosine'); % cosine dist

            if ii == jj
                dsm(ii,jj,tp) = 0;
            end
            
        end
    end
end

