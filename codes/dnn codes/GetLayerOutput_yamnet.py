# # #!/usr/bin/env python3
# # # -*- coding: utf-8 -*-
# # """
# # Created on Tue Sep 22 11:23:59 2020
# #!/usr/bin/env python
# # Instru to create the proper environment to run
# # the Yamnet network

# In[set environment]:
from __future__ import division, print_function

import sys

import pandas as pd
import numpy as np
import soundfile as sf
import resampy
import tensorflow as tf
from tensorflow.keras import Model, layers
import os
from copy import deepcopy #import copy as c
import h5py


# In[set paths etc.]:

workdir='/mnt/DATA/siddharth/myenv/workspace'
# dataset='formisano'
# dataset='giordano'
dataset='giordano'
models_dir='/mnt/DATA/siddharth/myenv/AcoSemDNN_Behav_fMRI_Repo/code/nlp_dnn_models' #from gio
dnn_dir=models_dir+'/audioset/yamnet/'
#dnn_weights_fn=dnn_dir+'yamnet.h5'
dnn_weights_fn='/mnt/DATA/siddharth/myenv/models/research/audioset/yamnet/'+'yamnet.h5' #default weights

stims_dir=workdir+'/stimuli/'
stims_list=stims_dir+'stimlist.csv'

out_dir=workdir+'/data/'

if not os.path.exists(out_dir):
    os.mkdir(out_dir)

wavs_dir=stims_dir

## In[read stims list]:
dat=pd.read_csv(stims_list,header=None)
dat=dat.drop(labels=0,axis=0)
stim_fns=np.asarray(dat[0])

out_fns=deepcopy(stim_fns);
for i in range(len(stim_fns)):
    out_fns[i]=out_dir+out_fns[i].replace('wav','hdf5')


# In[initialize yamnet]:

os.chdir(dnn_dir)
# from params import params as params
import params_yamnet_BLG as params
os.chdir('/mnt/DATA/siddharth/myenv/models/research/audioset/yamnet/')
import yamnet as yamnet_model
os.chdir(dnn_dir)

yamnet_params=params.Params()
#yamnet_params.patch_hop_seconds = 0.1  # 10 Hz scores frame rate.

yamnet = yamnet_model.yamnet_frames_model(yamnet_params)
yamnet.load_weights(dnn_weights_fn)
yam_layers=yamnet.layers
    



## In[]:

def pad_waveform(waveform,params):
    #do waveformpadding here
    min_waveform_seconds = (
        params.patch_window_seconds +
        params.stft_window_seconds - params.stft_hop_seconds)
    min_num_samples = tf.cast(min_waveform_seconds * params.sample_rate, tf.int32)
    num_samples = tf.shape(waveform)[0]
    num_padding_samples = tf.maximum(0, min_num_samples - num_samples)
    # In addition, there might be enough waveform for one or more additional
    # patches formed by hopping forward. If there are more samples than one patch,
    # round up to an integral number of hops.
    num_samples = tf.maximum(num_samples, min_num_samples)
    num_samples_after_first_patch = num_samples - min_num_samples
    
    hop_samples = tf.cast(params.patch_hop_seconds * params.sample_rate, tf.int32)
    num_hops_after_first_patch = tf.cast(tf.math.ceil(
        tf.cast(num_samples_after_first_patch, tf.float32) /
        tf.cast(hop_samples, tf.float32)), tf.int32)
    num_padding_samples += (
        hop_samples * num_hops_after_first_patch - num_samples_after_first_patch)
    num_padding_samples=int(num_padding_samples.numpy())
    waveform_padded=np.concatenate((waveform,np.zeros(num_padding_samples)))
    # print(num_samples.numpy())
    # print(num_samples_after_first_patch.numpy())
    # print(hop_samples.numpy())
    # print(min_waveform_seconds)
    # print(params.stft_hop_seconds)
    # print(params.stft_window_seconds)
    return waveform_padded



# In[]:


layer_nams=['layer1/relu',
            'layer2/pointwise_conv/relu',
            'layer3/pointwise_conv/relu',
            'layer4/pointwise_conv/relu',
            'layer5/pointwise_conv/relu',
            'layer6/pointwise_conv/relu',
            'layer7/pointwise_conv/relu',
            'layer8/pointwise_conv/relu',
            'layer9/pointwise_conv/relu',
            'layer10/pointwise_conv/relu',
            'layer11/pointwise_conv/relu',
            'layer12/pointwise_conv/relu',
            'layer13/pointwise_conv/relu',
            'layer14/pointwise_conv/relu',
            'global_average_pooling2d', #embedding
            'dense',    #logits
            'activation']#predictions
layer_nice=['layer01relu',
            'layer02relu',
            'layer03relu',
            'layer04relu',
            'layer05relu',
            'layer06relu',
            'layer07relu',
            'layer08relu',
            'layer09relu',
            'layer10relu',
            'layer11relu',
            'layer12relu',
            'layer13relu',
            'layer14relu',
            'embedding', #embedding
            'logits',    #logits
            'predictions']#predictions

for i  in range(len(out_fns)):
    wav_file=wavs_dir+stim_fns[i]
    
    print(wav_file)
    wav_data, sr = sf.read(wav_file, dtype=np.int16)
    assert wav_data.dtype == np.int16, 'Bad sample type: %r' % wav_data.dtype
    wave = wav_data / 32768.0  # Convert to [-1.0, +1.0]

    # Convert to mono and the sample rate expected by YAMNet.
    if len(wave.shape) > 1:
      wave = np.mean(wave, axis=1)
    if sr != yamnet_params.sample_rate:
      wave = resampy.resample(wave, sr, yamnet_params.sample_rate)
    
    wave_padded=pad_waveform(wave,yamnet_params)
    
    extractor = Model(inputs=yamnet.inputs,
                        outputs=[layer.output for layer in yamnet.layers])
    features = extractor(wave_padded)
    
    out=[];
    for ref in layer_nams:
        for l in range(len(yam_layers)):        
            if yam_layers[l].name.find(ref)>-1:#yam_layers[l].name==ref:
                tmp=np.asarray(features[l])
                out.append(tmp)
                #print(yam_layers[l].name)
                #print(tmp.shape)
    
    out_file=out_fns[i]
    print('saving data')
    if os.path.exists(out_file):
        os.remove(out_file)
    hf=h5py.File(out_file,'w')
    for j in range(len(out)):
        hf.create_dataset(layer_nice[j],data=out[j])
    hf.close()
    print('all done')
    
    
    

