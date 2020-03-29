# WFIRSTxCMB
WFIRST x CMB Forecast

work in progress

## Capabilities

Create a datavector with gaussian covaraince matrix build directly with cosmosis

```
cosmosis modules/WFIRSTxCMB/create_datavector_gaussian_covariance.ini
```

Create datavector using covariance matrix from cosmolike

```
cosmosis modules/WFIRSTxCMB/create_datavector.ini
```

Run forecast

```
cosmosis modules/WFIRSTxCMB/params.ini
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

* Also create a folder 6x2pt_WFIRST_SO in your cosmosis home directory.
mkdir 6x2pt_WFIRST_SO


* you have to add the module the makefile
nano modules/Makefile.modules
then add:
SUBDIRS+=WFIRSTxCMB
and save



