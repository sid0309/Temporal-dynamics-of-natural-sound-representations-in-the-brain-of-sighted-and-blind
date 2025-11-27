# Temporal-dynamics-of-natural-sound-representations-in-the-brain-of-sighted-and-blind

Temporal-dynamics-of-natural-sound-representations-in-the-brain-of-sighted-and-blind

PLEASE READ BEFORE RUNNING SCRIPTS

The codes refers to the 4 analyses parts in the study -

1) Preprocessing (preprocessing.m)

2) Decoding using preprocessed data (Run sections within param_decoding.m to run the following functions)
	a) Sensor space (fn_svm_decode_libsvm.m)
	b) Sensor Searchlight (fn_svm_decode_libsvm_chan.m)
	c) Temporal Generalization (fn_svm_decode_tempgen.m)
	d) Source Searchlight (fn_computesourcetime_search.m)

3) Creating models
	a) Modulation Transfer Function (mtf_model.m, also need nsl toolbox)
	b) YAMNet - Deep Neural Networks (adapted from Giordano et al., 2023 Nature Neuroscience)
	c) EEG models - Convert decoding AUC into DSMs at source and sensor space (dsm_eeg.m)

4) Representation Similarity Analyses (rsa.m): Uses precomputed EEG and external DSMs

For questions contact - siddharthtalwar0309@gmail.com
