# DNAffinity

#### Machine Learning method to predict binding sites for different Transcription Factors using data from different experimental techniques.

DNAffinity presents a physics-based machine learning approach to predict *in vitro* transcription factor binding affinities from structural and mechanical DNA properties directly derived from atomistic molecular dynamics simulations. The method is able to predict affinities obtained with high-throughput techniques as different as uPBM, gcPBM and HT-SELEX. When complemented with chromatin structure information, our *in vitro* trained method provides also good estimates of *in vivo* binding sites in yeast.


The pipeline for a specific protein consists on:

1. Getting the data file (i.e under the **test_data** folder the `Gata4_alignment_weighted.txt` file) 
2. Running the corresponging regressor (i.e. `python upbm_regressor.py`)

This is a brief summary of the main functions:

* Data import: `readData`, `processReads`
* Data Processing: `underSampling`, `overSampling`
* Model Training: `peakDetection`, `peakScoring`
* Model Testing: `plotPeaks`
* Predicting: `syntheticNucMap`

## Usage
---------------

This is an example of the main steps followed to predict affinities for a specific Transcription Factor.

1- Import the functions from the corresponding scripts

2- Load in the data specifying the experimental technique
&nbsp;

    from loadData import readData
    readData(file, technique)

