# Optimal design of experiments by combining coarse and fine measurements
This repository provides the Matlab script and jupyter notebook that accompany our paper entitled "Optimal design of experiments by combining coarse and fine measurements" ([arXiv:1702.06001](https://arxiv.org/abs/1702.06001)). 

## Abstract
In many contexts it is extremely costly to perform enough high quality experimental measurements to accurately parameterize a predictive quantitative model. However, it is often much easier to carry out large numbers of experiments that indicate whether each sample is above or below a given threshold. Can many such binary or “coarse” measurements be combined with a much smaller number of high resolution or “fine” measurements to yield accurate models? Here, we demonstrate an intuitive strategy, inspired by statistical physics, wherein the coarse measurements are used to identify the salient features of the data, while the fine measurements determine the relative importance of these features. We illustrate our strategy by considering the problem of solubility prediction for small organic molecules from their 2D molecular structure.

## Installation
The packages needed to run the scripts are:

1. rdkit (http://www.rdkit.org)
2. glmnet (http://web.stanford.edu/~hastie/glmnet_matlab/) 

## Running the scripts 

The jupyter notebook generates the binary chemical fingerprints using SMILES strings as input. The Matlab script analyses the binary chemical fingerprints. 
