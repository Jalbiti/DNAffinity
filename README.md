# DNAffinity

#### Machine Learning method to predict binding sites for different Transcription Factors using data from different experimental techniques.

DNAffinity presents a physics-based machine learning approach to predict *in vitro* 
transcription factor binding affinities from structural and mechanical DNA properties 
directly derived from atomistic molecular dynamics simulations. The method is able to 
predict affinities obtained with high-throughput techniques as different as uPBM, gcPBM 
and HT-SELEX. When complemented with chromatin structure information, our *in vitro* 
trained method provides also good estimates of *in vivo* binding sites in yeast.


The pipeline for a specific protein consists of:

1. Preprocessing the data (specific steps for each experimental technique; details on upbm and selex can be found under the **notebooks** folder)
2. Getting the preprocessed data file (e.g. the `Gata4_training.txt` file from the **test_data** folder) 
3. Running the corresponding regressor (e.g. `python upbm_regressor.py Gata4` for uPBM data, or `python selex_regressor.py Gata4 4` for selex data)

## Usage

---------------

This is an example of the main steps followed to predict affinities for a specific 
transcription factor. For testing purposes, just launch `{method}_regressor.py` to
run the calculations on the test case.

1. Make sure you have all the required files (compare formatting with the templates): 
   + `proteins/{protein}/{protein}_training.txt` for `upbm_regressor.py`
   + `proteins/{protein}/{cycle}.txt` for `selex_regressor.py`
   + `proteins/{protein}/{protein}_{concentration}.txt` and `proteins/{protein}/_freq_matrix_6.txt` for `gcpbm_regressor.py`
2. Launch `{method}_regressor.py` always passing the protein name as the 1st argument, and:
   + the cycle number as the 2nd argument for `selex_regressor.py`
   + the concentration as the 2nd argument for `gcpbm_regressor.py`

