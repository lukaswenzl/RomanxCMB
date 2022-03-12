# RomanxCMB
Nancy Grace Roman Space Telescope x CMB lensing Forecast

This repository contains the main components of the code used in Wenzl et al. 2021 https://arxiv.org/abs/2112.07681
Feel free to make use of this code under the condition of citing this paper in any resulting publications. If you make use of modules that are based on external codes please also acknowledge them.
Also, feel free to contact the corresponding author for additional assistance.

## Introduction

This code contains the main components to forecast parameter constraints on a 6x2pt analysis of Roman weak lensing and clustering data as well as CMB lensing data from Simons Observatory. The code is written as a range of modules.


## Requirements

* CosmoSIS https://bitbucket.org/joezuntz/cosmosis/
* cosmosis-standard-library https://bitbucket.org/joezuntz/cosmosis-standard-library/
recommended installation for both: use the provided Docker installation, see https://cosmosis.readthedocs.io/en/latest/installation/docker.html

* Bayesfast https://github.com/HerculesJack/bayesfast
* cosmolike (Non public code needed to calculate non-gaussian corrections to covariance matrix for new survey scenarios. Not included in this repo.)

External codes included in this repo due to required adaptations. Please make sure to cite the original source.

* HMcode (see mead folder) https://github.com/alexander-mead/HMcode
* SO noise models (see cmb folder) https://github.com/simonsobs/so_noise_models
* Eisenstein & Hu (1998, astro-ph/9710252) fitting formula (see eisenstein_hu_cdm folder)

additional details can be found in the module.yaml files in each module!



## Capabilities

Create a datavector with gaussian covaraince matrix build directly with cosmosis

```
cosmosis modules/RomanxCMB/create_datavector_gaussian_covariance.ini
```

Create datavector using covariance matrix from cosmolike

```
cosmosis modules/RomanxCMB/create_datavector.ini
```

Run forecast

```
cosmosis modules/RomanxCMB/params.ini
```



# How to use

## Installation

This code is written as a module for CosmoSIS. 
To install it, navigate to your cosmosis folder, then clone this repo into the modules folder:

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


To avoid crashes for non-existent folders create the standard output folder manually:

```bash
cd .. #back to cosmosis home directory.
mkdir 6x2pt_Roman_SO
```

This should make it possible to run the cosmosis pipeline to calculate a likelihood. Cosmosis comes with samplers like polychord that should work now. In the paper we use an optimized sampling approach to reduce computation time. To make this work you need to install bayesfast ( https://github.com/HerculesJack/bayesfast ) on your system and adapt the fisher and importance samplers. 
For the latter replace the files in cosmosis/samplers/fisher and cosmosis/samplers/importance with the contents of modified_fisher and modified_importance folders (these are slight changed that make the MOPED compression work and allow for thinned importance sampling).


## Example files

covariance matrix and other large files are not included in the repo.
Output files and other large files are not included in the github repo. Examples can be downloaded from TODO

## Run single sample

Ini files specify the full pipeline to calculate the observables based on a set of input parameters. Making use of a pre-calculated data vector and covariance the cosmosis also calculates the likelihood. For our forecast the fiducial data vector and the estimated covariance matrix are used.







