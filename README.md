# Fractional-dynamics-foster-deep-learning-of-COPDstage-prediction
## K-fold-validation and hold-out-validation  

"k_fold_validation/k_fold_fractional_model_train.ipynb" and "hold_out_validatio/hold_out_fractional_train.ipynb" describe our fractional dynamics deep learning model under k-fold and hold-out validation, where the inputs of these two files are "matrix A" files ("X_data_k_fold.npy" is the input file of "k_fold_fractional_model_train.ipynb"; and "vb_data.npy", "md1_data.npy", "md2_data.npy", and "cp_data.npy" are the input files of "hold_out_fractional_train.ipynb"). 

Please use this link to download our dataset (matrix A and raw signals) to execute our code: https://drive.google.com/drive/folders/1bSsQnvEm8DFJVicxSlkdPCwLciq2PE4v?usp=sharing

In this google drive "k-fold-validation" folder, "X_data_k_fold.npy" is the matrix A for all the samples, and "X_2_data.npy" is the raw signals. In the "hold-out-validation" folder, "vb_raw_data.npy", "md1_raw_data.npy", "md2_raw_data.npy", and "cp_raw_data.npy" are raw signals gathered from each institutions, respectively. "vb_data.npy", "md1_data.npy", "md2_data.npy", and "cp_data.npy" are the correlated "matrix A" files. 

Of note, "k_fold_fractional_model_train.ipynb" and "hold_out_fractional_train.ipynb" need 20 mins to execute. The other baseline models need at least 48 hours to be trained. 

The code package is available in Python. Run code in Google Colab is recommended.


## Fractional dynamics modeling of COPD relevant physiological signals
"Fractional dynamic modeling" folder contains the code to extract the fractional dynamic signatures from time-series. 

We provide an example edf file (signals) generated from a patient with mild symptoms (stage 1). Please use this link to download the file: https://drive.google.com/drive/folders/15VuS5EbrmnktFh9dsAX1NRQX-t-eGzV_?usp=sharing

The code package is written in Matlab. 

