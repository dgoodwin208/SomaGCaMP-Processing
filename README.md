# SomaGCaMP Processing scripts

This repository is a collection of scripts that we used to process the raw data in "Precision calcium imaging of dense neural populations via a cell body-targeted calcium indicator" by Shemesh et al, 2020. We call the soma-localized GCaMP somaGCaMP

# Requirements
Note: it requires that the user has the latest MATLAB CaImAn package installed, which is accessible at [this link](https://github.com/flatironinstitute/CaImAn-MATLAB), and maintained by the [FlatIron Institute](https://www.simonsfoundation.org/flatiron/). CaImAn is a robust CNMF-based package for the processing of GCaMP data, and was effectively used in SomaGCaMP analysis. We are grateful to Eftychios Pnevmatikakis and Andrea Giovannucci for their help in applying CaImAn to SomaGCaMP data.

# Overview: 
The user can start with postProcessCNMF.m: That will call processCNMF.m, which is what actually runs the CaImAn CNMF package with the parameters we used for the paper. 


