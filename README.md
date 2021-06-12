# RomanxCMB
Nancy Grace Roman Space Telescope x CMB Forecast

work in progress

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


Not needed anymore:
To switch between creating datavector and running forecast: 
for creating datavector turn on:
save_2pt
add_covariance
for runninng forcast instead use:
2pt_like
likelihoods = 2pt


covariance matrix and other large files are not included


# How to use
* written as a module for CosmoSIS. Clone github into modules/

* Also create a folder 6x2pt_Roman_SO in your cosmosis home directory.
mkdir 6x2pt_Roman_SO


* you have to add the module the makefile
nano modules/Makefile.modules
then add:
SUBDIRS+=RomanxCMB
and save

when using bayesfast

* replace the cosmosis/samplers/fisher and cosmosis/samplers/importance with the contents of modified_fisher and modified_importance

(these are slight changed that make the MOPED compression work and allow for thinned importance sampling)



