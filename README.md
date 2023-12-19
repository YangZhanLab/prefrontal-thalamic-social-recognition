This software package includes the codes used for population decoding of single-unit responses of freely-behaving mice to different social stimuli, as well as example data files for illustrating the running of each functions. The codes were written as MATLAB .m scripts, and tested with MATLAB R2023a under Windows 10 v22H2 operating system.

1. "mPFC_RE_decoding.m" was used to examine the decoding performance from mPFC and RE neurons. Population decoding was performed within each session, based on binned population vectors from the during stimulus period.  A naïve Bayesian decoder was trained to discriminate the firing of 4 different social stimuli. The bin length effect and population dimension effect on classification accuracies were compared for mPFC and RE. 

By default, this function will show the pre-computed multi-session accuracies for mPFC and RE, loaded from a file in the directory of "\data\decode_results_Eular_dist_4_class_all.mat".  To re-train and evaluate the model using the given example firing data, the running conditon statement in Line 34 should be set as "true". 

The example files include two sessions of single-unit firing for mPFC and RE each, saved respectively in the  directories of

              "\data\mPFC\produce_single_cell_data\"     and

              "\data\RE\produce_single_cell_data2020\".

It should be noted that, the pre-computed file will be overwritten after re-running the decoding. 

2. "single_trial_decoding.m" was used to examine the single-trial decoding performance between two experiment conditions, i.e., hm4d vs. control, based on a time-dependent cumulative Bayesian model. A spatial filter was first applied to the population vectors, and then the resultant latent representation was used for Bayesian likelihood estimation.  The two-class decoding performance, i.e., socal A vs. novel, were compared. The spatial filter was implemented as two subfunctions, called "SpRayleigh" and "TransformFea" respectively, and can be found in the latter part of the file. The SpRayleigh function was used to find the optimal spatial filter in a supervised manner by optimizing the Rayleigh quotient, and the TransformFea function was used to apply the trained filter to novel samples. 

By default, this function will show the pre-computed multi-session time curves of decoding accuracies of hm4d and control, saved in the files of "\data\decode_res_socialA_novel_bl23_high_var_cell_nbayes_accum_timemodel_pc1_triallim8.mat" and "\data\decode_res_socialA_novel_bl23_high_var_cell_nbayes_accum_timemodel_pc1_shuffle_triallim8.mat", where the former file gave the normal classification results, and the latter file gave the results when the model was trained on shuffled labels. To re-run the decoding on the given example files, the running conditon statement in Line 48 should be set as "true". 

The example files include two sessions of single-unit firing for hm4d and control each, saved respectively in the  directories of

              "\data\hm4d\produce_single_cell_data\"     and

              "\data\control\produce_single_cell_data\".

It should be noted that, the pre-computed files will be overwritten after re-running the decoding. 

