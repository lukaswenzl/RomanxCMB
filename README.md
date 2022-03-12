# RomanxCMB
Nancy Grace Roman Space Telescope x CMB lensing Forecast

This repository contains the main components of the code used in Wenzl et al. 2021 https://arxiv.org/abs/2112.07681
Feel free to make use of this code under the condition of citing this paper in any resulting publications. If you make use of modules that are based on external codes please also acknowledge them.
For additional assistance feel free to contact the corresponding author.

## Introduction

This code contains the main components to forecast parameter constraints on a 6x2pt analysis of Roman weak lensing and clustering data as well as CMB lensing data from Simons Observatory. The code is written as a range of modules.


## Requirements

* CosmoSIS https://bitbucket.org/joezuntz/cosmosis/
* cosmosis-standard-library https://bitbucket.org/joezuntz/cosmosis-standard-library/
* Bayesfast https://github.com/HerculesJack/bayesfast
* cosmolike (Non public code needed to calculate non-gaussian corrections to covariance matrix for new survey scenarios. Not included in this repo.)

External codes included in this repo due to required adaptations. Please make sure to cite the original source.

* HMcode (see mead folder) https://github.com/alexander-mead/HMcode
* SO noise models (see SO_noise_curves folder) https://github.com/simonsobs/so_noise_models
* Eisenstein & Hu (1998, astro-ph/9710252) fitting formula (see eisenstein_hu_cdm folder)

additional details can be found in the module.yaml file in each module!



## Capabilities

* Calculate observables of a 6x2pt analysis for a set of cosmological and nuisance parameters

* Create a fiducial datavector and estimate covariance matrix (gaussian contributions calculated with cosmosis. For the smaller non-gaussian contributions cosmolike is needed)

* Create forecast for parameter constraints of the analysis. Can use fisher, emcee, polychord, etc. included in cosmosis or use our optimized approach based on bayesfast (included script TODO)

This repo contains a range of modules that could be useful for other applications:

* Cosmosis wrapper for HMcode2020 (mead folder)
* Cosmosis wrapper for Eisenstein & Hu (1998, astro-ph/9710252) (eisenstein_hu_cdm folder)
* Modules to calculate the effects of a range of modified gravity models on the growth and CMB lensing kernel
* Modules to sample over sigma8 and sigma8_of_z  (sample_sigma8 and sample_sigma8_of_z folders)
* An implementation of the Intrinsic Alignment model in Krause, Eifler & Blazek 2016 
* various small helpful scrips in miscellaneous folder 

# How to use

## Installation

This code is written as a module for CosmoSIS. CosmoSIS can be installed in a variety of ways, a convenient option is the provided Docker installation, see https://cosmosis.readthedocs.io/en/latest/installation/docker.html

To install this code, navigate to your cosmosis folder, then clone this repo into the modules folder:

```bash
cd modules
git clone <link to this repo>
```



To compile the fortran and c codes contained in this repo together with the rest of cosmosis add the RomanxCMB folder: open the Makefile in the modules folder and add the Folder to it. The final file should look like this:

```
include ${COSMOSIS_SRC_DIR}/config/compilers.mk
include ${COSMOSIS_SRC_DIR}/config/subdirs.mk

SUBDIRS = RomanxCMB

-include Makefile.modules
```

Then recompile cosmosis in the cosmosis home directory with make

To avoid crashes for non-existent folders create the standard output folder manually:

```bash
cd .. #back to cosmosis home directory.
mkdir 6x2pt_Roman_SO
mkdir 6x2pt_Roman_SO_gaussian
```

This should make it possible to run the cosmosis pipeline to calculate a likelihood. Cosmosis comes with samplers like polychord that should work now. In the paper we use an optimized sampling approach to reduce computation time. To make this work you need to install bayesfast ( https://github.com/HerculesJack/bayesfast ) on your system and adapt the fisher and importance samplers. 
For the latter replace the files in cosmosis/samplers/fisher and cosmosis/samplers/importance with the contents of modified_fisher and modified_importance folders (these are slight changed that make the MOPED compression work and allow for thinned importance sampling).


## Example files

We provide a shared folder with the non-gaussian contributions to the covariance matrix to be able to run the code to build the full datavector. We also provide example datavectors to run the likelihood pipeline.
The folder can be found under https://drive.google.com/drive/folders/102bFeT5nDoF8aKmel2shVDkLQH5l53af?usp=sharing (created in 2022, should be functional for at least a few years.) 


## Run single sample

Ini files specify the full pipeline to calculate the observables based on a set of input parameters. Making use of a pre-calculated data vector and covariance matrix cosmosis also calculates the likelihood. For a forecast the fiducial data vector and the estimated covariance matrix are used.

To create the datavector for a <scenario> run the corresponding create_datavetor_<scenario>.ini file. Then, run the pipeline with the corresponding params_<scenario>.ini file.

In the paper we consider 3 survey scenarios: HLS optimistic corresponds to no extension; HLS conservative corresponds to the _pessim extension; Wide scenario corresponds to the _wide extension.

To make the code run without downloading the large files reference above the example below only calculate the gaussian covariance matrix.

First create a datavector.

```
cosmosis modules/RomanxCMB/create_datavector_gaussian_covariance.ini
```

Run the test sampler which calculates one likelihood and outputs all observables into a folder named 6x2pt_Roman_SO_gaussian

```
cosmosis modules/RomanxCMB/params.ini
```







