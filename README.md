# DNAffinity

#### Machine Learning method to predict binding sites for different Transcription Factors using data from different experimental techniques.

DNAffinity presents a physics-based machine learning approach to predict *in vitro* transcription factor binding affinities from structural and mechanical DNA properties directly derived from atomistic molecular dynamics simulations. The method is able to predict affinities obtained with techniques as different as uPBM, gcPBM and HT-SELEX. When complemented with chromatin structure information, our *in vitro* trained method provides also good estimates of *in vivo* binding sites in yeast.


The pipeline for a specific protein consists on:
Markup : 1. Getting the data file (i.e under the **test_data** folder the `Gata4_alignment_weighted.txt` file) 
	 2. Running the corresponging regressor (i.e. `python upbm_regressor.py`)


