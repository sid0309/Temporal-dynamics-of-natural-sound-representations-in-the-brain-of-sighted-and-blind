# Temporal dynamics of natural sound representations in the brain of sighted and blind

PLEASE READ BEFORE RUNNING SCRIPTS

To visualize decoding and RSA results, download data from https://osf.io/xu5ye/files. <br/>
The link includes a **derivatives** folder that contains the precomputed decoding results. These files are used by the plotting scripts to generate the decoding plots. (see point 5 below) <br/>
For RSA, the **others folder** in the link above contains all precomputed EEG and external models. The rsa.m script uses these files directly to run the analysis.

If you are interested in specific stages of analyses, they are performed in the following order-

1) Preprocessing (preprocessing.m). Preprocessed data can be found at osf link.

2) Decoding using preprocessed data (Run sections within param_decoding.m to run the following functions)
	* Sensor space (fn_svm_decode_libsvm.m)
	* Sensor Searchlight (fn_svm_decode_libsvm_chan.m)
	* Temporal Generalization (fn_svm_decode_tempgen.m)
	* Source Searchlight (fn_computesourcetime_search.m)

3) Creating models
	* Modulation Transfer Function (mtf_model.m, needs nsl toolbox)
	* YAMNet - Deep Neural Networks (adapted from Giordano et al., 2023 Nature Neuroscience)
	* EEG models - Convert decoding AUC into DSMs at source and sensor space (dsm_eeg.m)

4) Representation Similarity Analyses (rsa.m): Uses precomputed EEG and external DSMs and plots the results

5) Plotting (plot_decoding, plot_seachlight_sensor, plot_searchlight_src): Uses precomputed decoding outputs in the derivatives folder

For questions contact - siddharthtalwar0309@gmail.com
