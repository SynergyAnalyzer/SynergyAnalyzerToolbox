# SynergyAnalyzer Toolbox for MATLAB, Version 1.2

December, 2023 

by Marta Russo, Alessandro Scano, Cristina Brambilla and Andrea d'Avella

The [Synergy Analyzer Toolbox](https://github.com/SynergyAnalyzer/SynergyAnalyzerToolbox.git) for MATLAB is an open source software to extract muscle and kinematic synergies; see [LICENSE](https://github.com/SynergyAnalyzer/SynergyAnalyzerToolbox/blob/main/LICENSE) for the terms of the license (GNU GPL v3).

The paper "SynergyAnalyzer: A Matlab toolbox implementing mixed-matrix factorization to identify kinematic-muscular synergies" is available at this [link](https://doi.org/10.1016/j.cmpb.2024.108217).

## Description

This toolbox has been developed to extract motor synergies. Three algorithms for the extraction have been implemented: NMF (non-negative matrix factorization, Lee & Seung 1999) , PCA (principal component analysis) and MMF (mixed-matrix factorization, Scano et al. 2022). Synergies can be extracted on EMG data, kinematic data and mixed data, i.e. both EMG and kinematics.
The toolbox also provides some basics features to preprocess, plot and average data.

## How to Install
Add the toolbox folder to your MATLAB path and access the MATLAB documentation for further details.

## How to Use
Please refer to the [Getting Started](https://github.com/SynergyAnalyzer/SynergyAnalyzerToolbox/blob/main/GETTINGSTARTED.md#getting-started) for a workflow example and to the [Data Requirements](https://github.com/SynergyAnalyzer/SynergyAnalyzerToolbox/blob/main/DATAREQUIREMENTS.md#data-requirements) for how to prepare your dataset.  

## How to Cite
To support this toolbox and its authors, please cite the works listed below. Thanks very much for your support.

Marta Russo, Alessandro Scano, Cristina Brambilla, Andrea d'Avella,
SynergyAnalyzer: A Matlab toolbox implementing mixed-matrix factorization to identify kinematic-muscular synergies,
Computer Methods and Programs in Biomedicine, Volume 251, 2024, 108217, https://doi.org/10.1016/j.cmpb.2024.108217.

## Examples
For more details check the [Getting Started](https://github.com/SynergyAnalyzer/SynergyAnalyzerToolbox/blob/main/GETTINGSTARTED.md#getting-started) Documentation page. 
Otherwise, try to run ```demo_nmf_emg.m``` for an example of extracting muscle synergies, or ```demo_mmf_mix.m``` for an example of kinematic-muscular synergies.
