J/A+AS/122/399      K and evolutionary corrections (Poggianti, 1997)
================================================================================
K and evolutionary corrections from UV to IR
     Poggianti B.M.
    <Astron. Astrophys. Suppl. Ser. 122, 399 (1997)>
    =1997A&AS..122..399P      (SIMBAD/NED BibCode)
================================================================================
ADC_Keywords: Galaxies, photometry ; Models, evolutionary
Keywords: galaxies: evolution - galaxies: photometry -
          galaxies: distances and redshifts - galaxies: fundamental parameters -
          cosmology: miscellaneous

Description:
    K and evolutionary corrections are given for the E, Sa and Sc Hubble
    types for the following filters up to the redshift z=3:
    Johnson-Bessell & Brett photometric system:
         U, B, V, R, I, J, H, K filters
    Modified Thuan & Gunn system:
         gri filters
    Cousins system:
         Rc Ic filters
    Bj, Rf, In filters.

File Summary:
--------------------------------------------------------------------------------
 FileName   Lrecl  Records   Explanations
--------------------------------------------------------------------------------
ReadMe         80        .   This file
sed.dat       230      294  *Flux for different SED models (tables 3-29)
kcorrv.dat     40       36   K corrections in the V band from this paper
                              (table 30) and from Pence (1976ApJ...203...39P)
kcorr.dat      34      602  *K corrections in the different filters
                              (tables 31-24, 42, 45-47)
ecorr.dat      34      602  *Evolutionary corrections (tables 35-38, 43, 48-50)
filters.dat    15      515   Filter response functions (tables 39-41, 44, 51)
tables.tex     99     4802   LaTeX version of the tables
--------------------------------------------------------------------------------
Note on sed.dat: SED = spectral energy distributions;
  the fluxes are in arbitrary units.
Note on kcorr.dat and ecorr.dat: Corrections as a function of redshift (z):
  elliptical model with average solar metallicity with an e-folding time of
  1 Gyr (E) and 1.4 Gyr (E2), Sa and Sc models.
--------------------------------------------------------------------------------

Byte-by-byte Description of file: sed.dat
--------------------------------------------------------------------------------
   Bytes Format Units   Label   Explanations
--------------------------------------------------------------------------------
   1-  5  I5    0.1nm   Lam     Wavelength
   7- 13  F7.4  ---     logF03  Flux for model SED: Elliptical of age 15 Gyr
  15- 21  F7.4  ---     logF04  Flux for model SED: Sa of age 15 Gyr
  23- 29  F7.4  ---     logF05  Flux for model SED: Sa of age 15 Gyr
  31- 37  F7.4  ---     logF06  Flux for model SED: E  of age 13.2 Gyr (z=0.1)
  39- 45  F7.4  ---     logF07  Flux for model SED: E  of age 10.6 Gyr (z=0.3)
  47- 53  F7.4  ---     logF08  Flux for model SED: E  of age  8.7 Gyr (z=0.5)
  55- 61  F7.4  ---     logF09  Flux for model SED: E  of age  7.4 Gyr (z=0.7)
  63- 69  F7.4  ---     logF10  Flux for model SED: E  of age  5.9 Gyr (z=1.0)
  71- 77  F7.4  ---     logF11  Flux for model SED: E  of age  4.3 Gyr (z=1.5)
  79- 85  F7.4  ---     logF12  Flux for model SED: E  of age  3.4 Gyr (z=2.0)
  87- 93  F7.4  ---     logF13  Flux for model SED: E  of age  2.2 Gyr (z=3.0)
  95-101  F7.4  ---     logF14  Flux for model SED: Sa of age 13.2 Gyr
 103-109  F7.4  ---     logF15  Flux for model SED: Sa of age 10.6 Gyr
 111-117  F7.4  ---     logF16  Flux for model SED: Sa of age  8.7 Gyr
 119-125  F7.4  ---     logF17  Flux for model SED: Sa of age  7.4 Gyr
 127-133  F7.4  ---     logF18  Flux for model SED: Sa of age  5.9 Gyr
 135-141  F7.4  ---     logF19  Flux for model SED: Sa of age  4.3 Gyr
 143-149  F7.4  ---     logF20  Flux for model SED: Sa of age  3.4 Gyr
 151-157  F7.4  ---     logF21  Flux for model SED: Sa of age  2.2 Gyr
 159-165  F7.4  ---     logF22  Flux for model SED: Sc of age 13.2 Gyr
 167-173  F7.4  ---     logF23  Flux for model SED: Sc of age 10.6 Gyr
 175-181  F7.4  ---     logF24  Flux for model SED: Sc of age  8.7 Gyr
 183-189  F7.4  ---     logF25  Flux for model SED: Sc of age  7.4 Gyr
 191-197  F7.4  ---     logF26  Flux for model SED: Sc of age  5.9 Gyr
 199-205  F7.4  ---     logF27  Flux for model SED: Sc of age  4.3 Gyr
 207-213  F7.4  ---     logF28  Flux for model SED: Sc of age  3.4 Gyr
 215-221  F7.4  ---     logF29  Flux for model SED: Sc of age  2.2 Gyr
--------------------------------------------------------------------------------

Byte-by-byte Description of file: kcorrv.dat
--------------------------------------------------------------------------------
   Bytes Format Units   Label     Explanations
--------------------------------------------------------------------------------
   1-  4  F4.2  ---     z         Redshift
   6- 10  F5.3  mag     E         K correction for E (1Gyr) model
  12- 16  F5.3  mag     E/S0      ? K correction for E/S0 Pence observations
  18- 22  F5.3  mag     Sa        K correction for Sa model
  24- 28  F5.3  mag     Sab       ? K correction for Sab Pence observations
  30- 34  F5.3  mag     Sc        K correction for Sc model
  36- 40  F5.3  mag     Sbc       ? K correction for Sbc Pence observations
--------------------------------------------------------------------------------

Byte-by-byte Description of file: kcorr.dat ecorr.dat
--------------------------------------------------------------------------------
   Bytes Format Units   Label     Explanations
--------------------------------------------------------------------------------
   1-  4  F4.2  ---     z         Redshift
   6-  7  A2    ---     Filt      [UBVRIJHKgricjfn ] Filter
   8- 13  F6.3  mag     E         Correction in band Filt for E (1Gyr) model
  15- 20  F6.3  mag     E2        ?Correction in band Filt for E2 (1.4Gyr) model
  22- 27  F6.3  mag     Sa        ?Correction in band Filt for Sa model
  29- 34  F6.3  mag     Sc        Correction in band Filt for Sc model
--------------------------------------------------------------------------------

Byte-by-byte Description of file: filters.dat
--------------------------------------------------------------------------------
   Bytes Format Units   Label     Explanations
--------------------------------------------------------------------------------
   1-  2  A2    ---     Filt      [UB12VRIHJKgricjfn ] Filter
   4-  8  I5    0.1nm   Lambda    Wavelength
  10- 14  F5.3  ---     Resp      Response
--------------------------------------------------------------------------------

Acknowledgements: Bianca Poggianti <biancap@archimede.pd.astro.it>

History:
  * 06-Aug-1996: first version
  * 27-Mar-2000: evolutionary corrections in U and B added (ecorr.dat,
    table tab35-38)
  * 03-Jul-2000: addition of K and EC corrections for Cousins' filters
    Rc and Ic, and for filters corresponding to Bj, Rf and In colours.

References:
   Pence W., 1976ApJ...203...39P
================================================================================
(End)                                         Patricia Bauer [CDS]   06-Aug-1996
