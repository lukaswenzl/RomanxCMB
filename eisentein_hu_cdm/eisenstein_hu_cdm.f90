! taken from http://background.uchicago.edu/~whu/transfer/transferpage.html
! converted from f77 to f90 with https://fortran.uk/plusfortonline.php?id=140&id1=3&id2=1&id3=1&id4=-1&id5=1
!*==TFMDM.spg  processed by SPAG 6.72Dc at 19:10 on 19 Jul 2020
!	Driver and Transfer Function / Power Spectrum subroutines
!	       for Eisenstein & Hu astro-ph/97XXXXX.
!		3/5/98: typo in lambda age fixed 3/5/98
!		9/1/98: typo in normalization of open gravity wave models fixed
!		9/22/98:bug in sigmavtop (affecting for massive neutrinos) fixed
!
!	  Here is a sample input file that may be cut out and used
!
!X------------------------------------------------------------------------------
! 		1.0, 0.0, 0.0, 0.05	Omega_0,Omega_L,Omega_nu,Omega_b
! 		0.5, 2.726, 1.		h,T_cmb,N_nu
! 		0			z
! 		30,30			kmax,numk
!  		1			iz
! 	        1.0			n
! 		0			ipower
!X------------------------------------------------------------------------------
!
!       DRIVER Program Power
!
!	  INPUT:
!
!	    omega   --	matter density in critical units (no lambda)
!	    omega_l --  cosmological constant in critical units
!	    omega_nu--  massive neutrino density in critical units
!	    omega_b --  baryon density in critical units
!	    h       --  Hubble constant in 100 km/s/Mpc
!	    T_cmb   --  CMB Temperature
!	    N_nu    --  number of degenerate massive neutrinos
!	    z	    --  redshift at which to evaluate transfer function
!	    kmax    --  maximum k (h Mpc^{-1}) for transfer function tabulation
!	    numk    --  number of k's for evaluation
!
!
!	    tilt    --  n (Power spectrum tilt)
!	    ipower  --  if (n<1) and omega or omegal=0, option to
!			normalize (0) with power law inflationary
!			gravity waves (1) no gravity wavesr; else
!			no gravity waves
!	    iz      --  (0) evaluate quantities at default redshifts
!	                (1) evaluate quantities at specified redshift
!
!
!	  COMMON/GLOBALVARIABLES:  Variables in transfer function
!				   evaluation
! 	    (all real*8)
!
!       *   omhh         --  Omega*h*h
!	*   f_nu         --  Omega_nu/Omega
!	*   f_baryon     --  Omega_b/Omega
!	*   N_nu         --  number of degenerate massive neutrinos
!	    y_d          --  (1+z_eq)/(1+z_d)
!	    alpha_nu     --  small scale suppression
!	    beta_c       --  logarithmic correction
!	    sound_horizon--  sound horizon in Mpc
!	    theta_cmb    --  T_cmb/2.7K
!	*   omega        --  matter density in critical units
!	*   omegal       --  cosmological constant in critical units
!	    z_equality   --  redshift of matter-radiation equality
!
!	* = must be defined before call to TFset_parameters
!
!	COMMON/PASSSIGMA: variables needed for mass fluctuation
!	                  evalutation
!
!	    scale  -- filter scale in Mpc
!	    tilt   -- n, power spectrum tilt
!	    z      -- redshift of evaluation
!	    ipower -- switch (0) P_cbnu (1) P_cb
!
!	OUTPUT
!
!	    delta_H   -- COBE normalization for horizon-scale
!		         fluctuation amplitude
!	  Examples of Top-hat and Gaussian window mass fluctuations
!
!	    sigma_8   -- delta M/M on 8  h^-1 Mpc (using P_cbnu)
!	    sigma_50  -- delta M/M on 50 h^-1 Mpc (using P_cbnu)
!	    sigma_Lya -- delta M/M on 1/(24 sqrt(0mega*h*h)) Mpc scale
!			 (using P_cb) at z=3
!	    sigma_DLA -- delta M/M on 50km/s scale
!			 (using P_cbnu) at z=4
!	    omega_gas -- Upper limit on Omega_gas (z=4) via
!			 Press-Schechter with gaussian filter at
!			 50 km/s deltac=1.33
!
!	  Extra Freebies
!
!	    sigmav_40  -- velocity on 40 h^-1 Mpc tophat spheres
!			 (using P_cb)
!
!	    age        -- age in Gyr at z
!	    tau	       -- Thomson optical depth for full ionization to z
!	    growth     -- D1(z=0)/D1(z) growth from z to 0
!	    suppression-- alpha_nu*D_1**(-p_cb) suppression in CDM baryon transfer function
!			  as k->Infty at z
!
!      TRANSFER FUNCTION SUBROUTINES/FUNCTIONS
!
!	TF_setparameters [subroutine]
!
!	  Uses COMMON/GLOBALVARIABLES (above)
!	  Sets unstarred variables in GLOBALVARIABLES list above
!
!	TF_master(k) [real*8 function]
!
!	  Argument: k wavenumber in Mpc^{-1} (IMPORTANT: no h)
!	  Uses COMMON/GLOBALVARIABLES (above)
!
!	Growth [subroutine]
!
!	  Arguments:
!
!	      z       -- redshift to evaluate growth
!	      k       -- wavenumber in Mpc^{-1} (IMPORTANT: no h)
!	      DD_cb   -- on exit set to D_cb(z)/D_1(z=0)
!	      DD_cbnu -- on exit set to D_cbnu(z)/D_1(z=0)
!	      DD0     -- on exit set to D_1(z=0)
!
!	  Uses COMMON/GLOBALVARIABLES (above)
!
!	sound_horizon [real*8 function]: approximate sound horizon
!	  [a convenient approximation that does not require the
!	   parameter set up of TF_setparameters; unnecessary if
!	   TF_setparameters is called.
!
!	  Arguments:
!
!	      omhh     -- Omega*h*h
!	      f_baryon -- Omega_b/Omega
!
!       sigmatop [real*8 function]: evaluates integrand of sigma for
! 	  a top hat filter (unnormalized)
!
!	  Arguments:
!
!	      kl       -- ln(k Mpc^{-1})  [real*8]
!
!	  Uses COMMON/PASSSIGMA
!
!       sigmagaus [real*8 function]: evaluates integrand of sigma for
! 	  a gaussian filter (unnormalized)
!
!	  Uses COMMON/PASSSIGMA
!
!	sigmavtop [real*8 function]: evaluates integrand of sigmav for
!	  a top hat filter (unnormalized)
!
!	  Uses COMMON/PASSSIGMA
!	       COMMON/GLOBALVARIABLES
!
!	  Arguments:
!
!	      kl       -- ln(k Mpc^{-1})  [real*8]
!
!	  Uses COMMON/PASSSIGMA
!
!	rhombint [real*8 function]: Rhomberg integration
!				    (taken from CMBfast)
!
!	  Arguments:
!
!	      f       -- integrand externally defined
!	      a       -- lower limit of integration [real*8]
!	      b       -- upper limit of integration [real*8]
!	      tol     -- tolerance [real*8]
!
!	erfcc [real*8 function]: complementary error function
!				 (taken from N*merical R*cipes)
!
!	  Arguments:
!
!	      x       -- real*8
!
 
      MODULE ehu_cdm
 
        IMPLICIT NONE

        public calculate_transfer
        public GROWTH
        contains

        SUBROUTINE calculate_transfer(input_OMEga, input_OMEgal, omeganu, omegab, h, t_cmb, &
            & input_N_Nu, kmin, kmax, numk, transfer, k_arr, sigma_8) !todo Lukas
      !input_OMEga , input_OMEgal , omeganu , omegab, h , t_cmb , input_N_Nu, kmax , numk, Z 
      ! parameters: Omega_0, Omega_Lambda, Omega_nu, Omega_b, h, T_cmb, N_nu, k_max (h Mpc^{-1}) (default 10), numk #pts per decade (default 50), Z redshift for transfer function, 
  !*--TFMDM180
  !*** Start of declarations inserted by SPAG
        REAL*8 age , ALPha_nu , anorm , arg , BETa_c , dd0 , dd_cb ,      &
             & dd_cbnu , deltac , dum1 , dum2 , dum3 , F_Baryon , &
             & F_Nu , g , k
      !moved kmax
        INTEGER, INTENT(INOUT) :: numk
        INTEGER i , IPOwer , iz
        REAL*8, INTENT(INOUT) :: kmin, kmax, input_OMEga, input_OMEgal, omeganu, &
             & omegab, h, t_cmb, input_N_Nu, sigma_8
        REAL*8, dimension(:), allocatable, intent(INOUT) :: transfer, k_arr
        REAL*8 OMEga , omegae , omegagas , OMEgal ,    &
             & omegaz , OMHh , SCAle , sigmadla ,   &
             & sigmav_40 , sigma_50 ,    &
             & sigma_lya , SOUnd_horizon , tau
        !removed Lukas: ROMBINT, SIGMAGAUS, SIGMATOP, SIGMAVTOP, TF_MASTER, ERFCC
        REAL*8 THEta_cmb , TILt , tol , t_master ,    &
             & vcirc , xe , yp , Y_D , Z , Z_Equality
  !*** End of declarations inserted by SPAG
        !removed: EXTERNAL ROMBINT , SIGMATOP , SIGMAGAUS , SIGMAVTOP
        REAL*8 N_Nu
        COMMON /GLOBALVARIABLES/ OMHh , F_Nu , F_Baryon , N_Nu , Y_D ,    &
                               & ALPha_nu , BETa_c , SOUnd_horizon ,      &
                               & THEta_cmb , OMEga , OMEgal , Z_Equality
        COMMON /PASSSIGMA/ SCAle , TILt , Z , IPOwer
   
        

      !save inputs into global variables
        OMEga = input_OMEga 
        OMEgal = input_OMEgal 
        N_Nu = input_N_Nu
        !t_cmb = 2.726
        !kmax = 10.
        !numk = 50 
        Z = 0.
        allocate(transfer(numk))
        allocate(k_arr(numk))
        
   
  !	Translate Parameters into forms GLOBALVARIABLES form
   
        F_Nu = omeganu/OMEga
        F_Baryon = omegab/OMEga
        IF ( F_Nu.EQ.0 ) F_Nu = 1E-10
        IF ( F_Baryon.EQ.0 ) F_Baryon = 1E-10
        IF ( t_cmb.LE.0 ) t_cmb = 2.726
        IF ( N_Nu.LT.0 ) N_Nu = 0.
        IF ( h.GT.10 ) WRITE (6,*) 'WARNING: H_0/100km/s/Mpc needed'
        THEta_cmb = t_cmb/2.7
        OMHh = OMEga*h*h
   
        WRITE (6,*) '  '
        WRITE (6,*) 'Confirming Cosmological Parameters'
   
        WRITE (6,99001) OMEga , OMEgal , F_Nu , F_Baryon
  99001 FORMAT ('   Omega_0=',F5.3,'  Omega_Lambda=',F5.3,'  f_nu=',F5.3, &
               &'  f_baryon=',F5.3)
        WRITE (6,99002) h*100 , t_cmb , N_Nu
  99002 FORMAT ('   H_0    =',F5.1,'  T_cmb       =',F5.3,'K',' N_nu=',   &
              & F4.2)
   
   
        CALL TFSET_PARAMETERS
   
  !	sound_horizon=sound_horizon_fit(omhh,f_baryon)
  !  what redshift?
   
        !WRITE (6,*) '  '
        !WRITE (6,*) 'Redshift for Transfer Function'
        !READ * , Z
   
        OPEN (10,FILE='trans.dat')
   
        !WRITE (6,*) 'k_max (h Mpc^{-1}),#pts per decade (10,50)'
        !READ * , 
   
        IF ( kmax.LE.0 ) kmax = 10.
        IF ( numk.LE.0 ) numk = 50
   
        !kmin = 0.0001
        kmin = 1.E-5 !changed to lower value -Lukas
        !numk = numk*LOG10(kmax/kmin)
   
        WRITE (6,*) '  '
        WRITE (6,99003) Z
  99003 FORMAT ('Writing  [k(h Mpc^-1), T_master, T_cb, T_cbnu](z=',F5.2, &
               &')  to file: trans.dat...')
        WRITE (6,*) '  '
        DO i = 1 , numk
   
           k = 10.**(i*(LOG10(kmax/kmin)/numk))*kmin
           CALL GROWTH(Z,k*h,dd_cb,dd_cbnu,dd0)
           t_master = TF_MASTER(k*h)
  !
           WRITE (10,99004) k , t_master , t_master*dd_cb ,               &
                          & t_master*dd_cbnu
            print *, 'index ', i
            transfer(i) = t_master*dd_cbnu
            k_arr(i) = k
   
  99004    FORMAT (1X,20E13.5)
   
        ENDDO
   
   
  !	Power Spectrum Statistics
   
        WRITE (6,*) 'Calculating power spectrum statistics...'
   
        !WRITE (6,99005) Z
  !99005 FORMAT (' Use redshift z=',F5.2,                                  &
        !       &' globally (0) no, use defaults (1) yes')
        !READ * , iz
   
        !WRITE (6,*) 'spectral index?'
   
        !IF ( ((1.-OMEgal-OMEga)*OMEgal).GT.1E-5 ) THEN
   
  !	  COBE normalization for general lambda + open universe
  !	  available only for tilt=1
   
        !   TILt = 1.
   
        !ELSE
  !
        !   READ * , TILt
   
        !ENDIF

        !Lukas: hardcoding spectral index to 1
        iz = 1
        TILt = 1.
   
   
   
        IF ( ABS(1-OMEga-OMEgal).LT.1E-5 ) THEN
   
           IF ( TILt.LT.1 ) THEN
              WRITE (6,*) '  Gravity Waves (power law inflation)?' ,      &
                         &' (0) yes   (1) no'
              READ * , IPOwer
           ELSE
              IPOwer = 1
              WRITE (6,*) ' '
           ENDIF
   
           WRITE (6,*)                                                    &
                    &'     ----------------------------------------------'
           IF ( IPOwer.EQ.0 ) THEN
              anorm = 1.94E-5*OMEga**(-0.785-0.05*LOG(OMEga))             &
                    & *EXP(1.*(TILt-1.)+1.97*(TILt-1.)**2)
              WRITE (6,99006) anorm , TILt
  99006       FORMAT ('delta_H:   ',E12.5,'    n=',F5.3,                  &
                     &' (B&W, flat, gravity waves)')
              anorm = anorm**2*(2997.9/h)**(3.+TILt)
           ELSE
              anorm = 1.94E-5*OMEga**(-0.785-0.05*LOG(OMEga))             &
                    & *EXP(-0.95*(TILt-1.)-0.169*(TILt-1.)**2)
              WRITE (6,99007) anorm , TILt
  99007       FORMAT ('delta_H:   ',E12.5,'    n=',F5.3,                  &
                     &' (B&W, flat, no gravity waves)')
              anorm = anorm**2*(2997.9/h)**(3.+TILt)
           ENDIF
   
   
        ELSEIF ( OMEgal.EQ.0 ) THEN
   
           IF ( TILt.LT.1 ) THEN
              WRITE (6,*) '  Gravity Waves (power law inflation)?' ,      &
                         &'(0) yes   (1) no'
              READ * , IPOwer
           ELSE
              IPOwer = 1
              WRITE (6,*) ' '
           ENDIF
           WRITE (6,*)                                                    &
                    &'     ----------------------------------------------'
   
           IF ( IPOwer.EQ.0 ) THEN
              anorm = 1.95E-5*OMEga**(-0.35-0.19*LOG(OMEga)+0.15*(TILt-1))&
                    & *EXP(1.02*(TILt-1.)+1.70*(TILt-1.)**2)
              WRITE (6,99008) anorm , TILt
  99008       FORMAT ('delta_H:   ',E12.5,'    n=',F5.3,                  &
                     &' (H&W, open, min. gravity waves)')
              anorm = anorm**2*(2997.9/h)**(3.+TILt)
           ELSE
              anorm = 1.95E-5*OMEga**(-0.35-0.19*LOG(OMEga)-0.17*(TILt-1))&
                    & *EXP(-1.*(TILt-1.)-0.14*(TILt-1.)**2)
              WRITE (6,99009) anorm , TILt
  99009       FORMAT ('delta_H:   ',E12.5,'    n=',F5.3,                  &
                     &' (B&W, open, no gravity waves)')
              anorm = anorm**2*(2997.9/h)**(3.+TILt)
           ENDIF
   
        ELSE
   
           anorm = 2.422 - 1.166*EXP(OMEga) + 0.800*EXP(OMEgal)           &
                 & + 3.780*OMEga - 2.267*OMEga*EXP(OMEgal)                &
                 & + 0.487*OMEga**2 + 0.561*OMEgal +                      &
                 & 3.392*OMEgal*EXP(OMEga) - 8.568*OMEga*OMEgal +         &
                 & 1.080*OMEgal**2
           anorm = 1.E-5*anorm
           WRITE (6,99010) anorm , TILt
  99010    FORMAT ('delta_H:   ',E12.5,'    n=',F5.3,                     &
                  &' (B&W, open, lambda, no gravity waves)')
           anorm = anorm**2*(2997.9/h)**(3.+TILt)
  !
        ENDIF
   
        tol = 1.E-6
   
   
  !	Sigma_8: Mass fluctuations in a top hat of 8./h Mpc at
  !	         redshift z=0; total mass: CDM+baryon+neutrino
   
        SCAle = 8./h
        IF ( iz.EQ.0 ) Z = 0
        IPOwer = 0
   
        sigma_8 = SQRT                                                    &
                & (anorm*(ROMBINT(SIGMATOP,DLOG(0.001/SCAle),DLOG(0.1D0/  &
                & SCAle),tol)+ROMBINT(SIGMATOP,DLOG(0.1D0/SCAle),         &
                & DLOG(1D0/SCAle),tol)                                    &
                & +ROMBINT(SIGMATOP,DLOG(1D0/SCAle),DLOG(10.D0/SCAle),tol)&
                & +ROMBINT(SIGMATOP,DLOG(10.D0/SCAle),DLOG(100.D0/SCAle), &
                & tol)))
   
        WRITE (6,99011) sigma_8 , Z
  99011 FORMAT ('sigma_8:   ',E12.5,                                      &
              & '                                  z=',F5.2)
   
   
   
  !	Sigma_50: Mass fluctuations in a top hat of 50./h Mpc at
  !	         redshift z=0; total mass: CDM+baryon+neutrino
   
        SCAle = 50./h
        IF ( iz.EQ.0 ) Z = 0
        IPOwer = 0
   
        sigma_50 = SQRT                                                   &
                 & (anorm*(ROMBINT(SIGMATOP,DLOG(0.001/SCAle),DLOG(0.1D0/ &
                 & SCAle),tol)+ROMBINT(SIGMATOP,DLOG(0.1D0/SCAle),        &
                 & DLOG(1D0/SCAle),tol)                                   &
                 & +ROMBINT(SIGMATOP,DLOG(1D0/SCAle),DLOG(10.D0/SCAle),   &
                 & tol)+ROMBINT(SIGMATOP,DLOG(10.D0/SCAle),               &
                 & DLOG(100.D0/SCAle),tol)))
   
        WRITE (6,99012) sigma_50 , sigma_50/sigma_8 , Z
  99012 FORMAT ('sigma_50:  ',E12.5,'    sigma_50/sigma_8:',E12.5,' z=',  &
              & F5.2)
   
   
   
  !	Sigma_DLA:  Mass fluctuations in a gaussian of
  !	            scale corresponding to 50km/s at
  !	            redshift z=4; total mass: CDM+baryon+neutrino mass
   
   
  !	circular velocity cutoff
        vcirc = 50.
        IF ( iz.EQ.0 ) Z = 4
        SCAle = 0.3834*(vcirc)/(h*100.)                                   &
              & *(OMEga*(1.+Z)**(3.)+(1.-OMEga-OMEgal)*(1.+Z)**2+OMEgal)  &
              & **(-1./6.)
        IPOwer = 0
   
        sigmadla = SQRT                                                   &
                 & (anorm*(ROMBINT(SIGMAGAUS,DLOG(0.001/SCAle),DLOG(0.1D0/&
                 & SCAle),tol)+ROMBINT(SIGMAGAUS,DLOG(0.1D0/SCAle),       &
                 & DLOG(1D0/SCAle),tol)                                   &
                 & +ROMBINT(SIGMAGAUS,DLOG(1D0/SCAle),DLOG(10.D0/SCAle),  &
                 & tol)+ROMBINT(SIGMAGAUS,DLOG(10.D0/SCAle),              &
                 & DLOG(20.D0/SCAle),tol)))
   
        WRITE (6,99013) sigmadla , SCAle*h , Z
  99013 FORMAT ('sigma_DLA: ',E12.5,'    at',E12.5,' h^-1 Mpc,      z=',  &
              & F5.2)
   
  !	Upper limit on Omega_gas via Press-Schecter
  !	Press-Schechter threshold
   
        deltac = 1.33
        omegagas = OMEga*F_Baryon*ERFCC(deltac/(SQRT(2.)*sigmadla))
   
        WRITE (6,99014) omegagas , deltac , Z
  99014 FORMAT ('Omega_gas:<',E12.5,'    at ',F5.3,                       &
               &' gaus. deltac,        z=',F5.2)
   
  !	Sigma_Lya:  Mass fluctuations in a gaussian of
  !	            1/(24.04 sqrt(0mega h**2)) Mpc at
  !	            redshift z=3; CDM+baryon mass
   
        SCAle = 1./(24.04*SQRT(OMHh))
   
        IF ( iz.EQ.0 ) Z = 3
        IPOwer = 1
   
        sigma_lya = SQRT                                                  &
                  & (anorm*(ROMBINT(SIGMAGAUS,DLOG(0.001/SCAle),DLOG(0.1D0&
                  & /SCAle),tol)+ROMBINT(SIGMAGAUS,DLOG(0.1D0/SCAle),     &
                  & DLOG(1D0/SCAle),tol)                                  &
                  & +ROMBINT(SIGMAGAUS,DLOG(1D0/SCAle),DLOG(10.D0/SCAle), &
                  & tol)+ROMBINT(SIGMAGAUS,DLOG(10.D0/SCAle),             &
                  & DLOG(20.D0/SCAle),tol)))
        WRITE (6,99015) sigma_lya , SCAle*h , Z
  99015 FORMAT ('sigma_Lya: ',E12.5,'    at',E12.5,' h^-1 Mpc,      z=',  &
              & F5.2)
   
   
  !	Sigmav_40:  Rms velocity in top hat of 100h^{-1} Mpc, z=0
  !		     default; CDM+baryon
   
        SCAle = 40./h
        IF ( iz.EQ.0 ) Z = 0
        IPOwer = 1
   
        sigmav_40 = SQRT                                                  &
                  & (anorm*(ROMBINT(SIGMAVTOP,DLOG(0.001/SCAle),DLOG(0.1D0&
                  & /SCAle),tol)+ROMBINT(SIGMAVTOP,DLOG(0.1D0/SCAle),     &
                  & DLOG(1D0/SCAle),tol)                                  &
                  & +ROMBINT(SIGMAVTOP,DLOG(1D0/SCAle),DLOG(10.D0/SCAle), &
                  & tol)+ROMBINT(SIGMAVTOP,DLOG(10.D0/SCAle),             &
                  & DLOG(100.D0/SCAle),tol)))
        WRITE (6,99016) sigmav_40 , SCAle*h , Z
  99016 FORMAT ('sigma_v:    ',F6.2,' km/s    at ',F6.2,                  &
               &' h^-1 Mpc,           z=',F5.2)
   
        WRITE (6,*) '     ----------------------------------------------'
        WRITE (6,*) 'Extras:'
        WRITE (6,*) '  '
   
        IF ( iz.EQ.0 ) Z = 0.
   
  !	  yp = primordial Helium mass fraction
        yp = 0.24  !todo Lukas
  !	  xe = ionization fraction 1                 = fully ionized
  !			           (1-yp)/(1-yp/2)   = HeI
  !			           (1-3yp/4)/(1-yp/2)= HeII
        xe = (1.-yp)/(1.-yp/2)
  !	  neutral helium
        tau = 4.61E-2*(1.-yp/2.)*xe*omegab*h/OMEga**2
   
        IF ( OMEgal.EQ.0 ) THEN
           omegae = OMEga - 1.D-6
           arg = 2.D0*(1.D0-omegae)/omegae/(1.+Z) + 1.D0
           age = 9.7778/h*(SQRT(1.+omegae*Z)/(1.-omegae)/(1.+Z)           &
               & -omegae/2./(1.-omegae)**(1.5)                            &
               & *DLOG(arg+DSQRT(arg**2-1.D0)))
           WRITE (6,99020) Z , age
   
           tau = tau*(2.-3.*OMEga+SQRT(1.+OMEga*Z)*(OMEga*Z+3.*OMEga-2.))
   
           IF ( Z.GT.0 ) WRITE (6,99021) Z , tau
   
  !
        ELSEIF ( ABS(1.-OMEga-OMEgal).LT.1E-10 ) THEN
           g = SQRT(OMEga*(1.+Z)**3+1.-OMEga)
           omegaz = OMEga*(1.+Z)**3/g**2
  !
           age = 2./3.*9.7778/h/g/SQRT(1.-omegaz)                         &
               & *LOG((1.+SQRT(1.-omegaz))/SQRT(omegaz))
           WRITE (6,99020) Z , age
   
           tau = tau*OMEga*(SQRT(1.-OMEga+OMEga*(1.+Z)**3)-1.)
   
           IF ( Z.GT.0 ) WRITE (6,99021) Z , tau
  !
   
        ENDIF
   
        IF ( Z.GT.0 ) WRITE (6,99017) Z , 1./dd0
   
  99017 FORMAT ('Growth(z=',F5.2,'-0):',F5.2)
   
   
        WRITE (6,99018) SOUnd_horizon*h
  99018 FORMAT ('sound_horizon:  ',F7.2,' h^-1 Mpc')
   
   
        CALL GROWTH(Z,1.D5,dum1,dum2,dum3)
        WRITE (6,99019) ALPha_nu*dum1 , Z
  99019 FORMAT ('suppression:    ',F7.5,' z=',F5.2)
   
  99020 FORMAT ('Age(z=',F5.2,'):    ',F6.2,' Gyrs')
  99021 FORMAT ('tau(z=',F5.2,'):    ',F9.5,' (HeI)')
        
      END SUBROUTINE calculate_transfer
      
  !*==TFSET_PARAMETERS.spg  processed by SPAG 6.72Dc at 19:10 on 19 Jul 2020
   
   
        SUBROUTINE TFSET_PARAMETERS
        IMPLICIT NONE
  !*--TFSET_PARAMETERS557
  !*** Start of declarations inserted by SPAG
        REAL*8 ALPha_nu , BETa_c , F_Baryon , f_c , f_cb , F_Nu , f_nub , &
             & k_equality , obhh , OMEga , OMEgal , OMHh , p_c , p_cb ,   &
             & r_drag , r_equality , SOUnd_horizon , THEta_cmb , Y_D ,    &
             & z_drag
        REAL*8 Z_Equality
  !*** End of declarations inserted by SPAG
  !
        REAL*8 N_Nu
        COMMON /GLOBALVARIABLES/ OMHh , F_Nu , F_Baryon , N_Nu , Y_D ,    &
                               & ALPha_nu , BETa_c , SOUnd_horizon ,      &
                               & THEta_cmb , OMEga , OMEgal , Z_Equality
   
  ! Auxiliary variable
   
        obhh = OMHh*F_Baryon
  ! Main variables
   
        Z_Equality = 2.50E4*OMHh*THEta_cmb**(-4.) - 1.D0
        k_equality = 0.0746*OMHh*THEta_cmb**(-2.)
   
        z_drag = 0.313*OMHh**(-0.419)*(1.+0.607*OMHh**(0.674))
        z_drag = 1E0 + z_drag*obhh**(0.238*OMHh**(0.223))
        z_drag = 1291E0*OMHh**(0.251)/(1E0+0.659*OMHh**(0.828))*z_drag
   
        Y_D = (1.+Z_Equality)/(1.+z_drag)
   
        r_drag = 31.5*obhh*THEta_cmb**(-4.)*1000E0/(1E0+z_drag)
        r_equality = 31.5*obhh*THEta_cmb**(-4.)*1000E0/(1E0+Z_Equality)
   
        SOUnd_horizon = 2./3./k_equality*SQRT(6./r_equality)              &
                      & *LOG((SQRT(1.+r_drag)+SQRT(r_drag+r_equality))    &
                      & /(1.+SQRT(r_equality)))
   
        p_c = -(5.-SQRT(1.+24*(1.-F_Nu-F_Baryon)))/4.
        p_cb = -(5.-SQRT(1.+24*(1.-F_Nu)))/4.
        f_c = 1. - F_Nu - F_Baryon
        f_cb = 1. - F_Nu
        f_nub = F_Nu + F_Baryon
   
   
        ALPha_nu = (f_c/f_cb)*(2.*(p_c+p_cb)+5.)/(4.*p_cb+5)
        ALPha_nu = ALPha_nu*(1.-0.553*f_nub+0.126*f_nub**3)
        ALPha_nu = ALPha_nu/(1.-0.193*SQRT(F_Nu)+0.169*F_Nu)
        ALPha_nu = ALPha_nu*(1.+Y_D)**(p_c-p_cb)
        ALPha_nu = ALPha_nu*(1.+(p_cb-p_c)                                &
                 & /2.*(1.+1./(4.*p_c+3.)/(4.*p_cb+7.))/(1.+Y_D))
        BETa_c = 1./(1.-0.949*f_nub)
   
   
        END
  !*==TF_MASTER.spg  processed by SPAG 6.72Dc at 19:10 on 19 Jul 2020
   
   
        REAL*8 FUNCTION TF_MASTER(K)
   
        IMPLICIT NONE
  !*--TF_MASTER615
  !*** Start of declarations inserted by SPAG
        REAL*8 ALPha_nu , BETa_c , F_Baryon , F_Nu , gamma_eff , K ,      &
             & OMEga , OMEgal , OMHh , q , q_eff , q_nu , SOUnd_horizon , &
             & THEta_cmb , Y_D , Z_Equality
  !*** End of declarations inserted by SPAG
        REAL*8 N_Nu
        COMMON /GLOBALVARIABLES/ OMHh , F_Nu , F_Baryon , N_Nu , Y_D ,    &
                               & ALPha_nu , BETa_c , SOUnd_horizon ,      &
                               & THEta_cmb , OMEga , OMEgal , Z_Equality
   
        q = K*THEta_cmb**2/OMHh
        gamma_eff = (SQRT(ALPha_nu)+(1.-SQRT(ALPha_nu))                   &
                  & /(1.+(0.43*K*SOUnd_horizon)**4))
   
        q_eff = q/gamma_eff
        TF_MASTER = DLOG(DEXP(1.D0)+1.84*BETa_c*SQRT(ALPha_nu)*q_eff)
        TF_MASTER = TF_MASTER/                                            &
                  & (TF_MASTER+q_eff**2*(14.4+325./(1.+60.5*q_eff**1.11)))
   
        q_nu = 3.92*q*SQRT(N_Nu/F_Nu)
        TF_MASTER = TF_MASTER*(1.+(1.2*F_Nu**(0.64)*N_Nu**(0.3+0.6*F_Nu)) &
                  & /(q_nu**(-1.6)+q_nu**(0.8)))
   
   
        END
  !*==GROWTH.spg  processed by SPAG 6.72Dc at 19:10 on 19 Jul 2020
   
        SUBROUTINE GROWTH(Z,K,Dd_cb,Dd_cbnu,Dd0)
   
        IMPLICIT NONE
  !*--GROWTH646
  !*** Start of declarations inserted by SPAG
        REAL*8 ALPha_nu , BETa_c , d , Dd0 , Dd_cb , Dd_cbnu , F_Baryon , &
             & F_Nu , K , olz , OMEga , OMEgal , OMHh , oz , p_cb , q ,   &
             & SOUnd_horizon , THEta_cmb , Y_D , y_fs
        REAL*8 Z , Z_Equality
  !*** End of declarations inserted by SPAG
        REAL*8 N_Nu
   
        COMMON /GLOBALVARIABLES/ OMHh , F_Nu , F_Baryon , N_Nu , Y_D ,    &
                               & ALPha_nu , BETa_c , SOUnd_horizon ,      &
                               & THEta_cmb , OMEga , OMEgal , Z_Equality
   
        q = K*THEta_cmb**2/OMHh
   
        y_fs = 17.2*F_Nu*(1.+0.488*F_Nu**(-7./6.))*(N_Nu*q/F_Nu)**2
   
        oz = OMEga*(1.+Z)                                                 &
           & **3/(OMEgal+(1.-OMEgal-OMEga)*(1.+Z)**2+OMEga*(1.+Z)**3)
        olz = OMEgal/(OMEgal+(1.-OMEgal-OMEga)*(1.+Z)**2+OMEga*(1.+Z)**3)
   
        d = (1.+Z_Equality)/(1.+Z)                                        &
          & *5.*oz/2.*(oz**(4./7.)-olz+(1.+oz/2.)*(1.+olz/70.))**(-1.)
   
        Dd0 = d/((1.+Z_Equality)                                          &
            & *5.*OMEga/2.*(OMEga**(4./7.)-OMEgal+(1.+OMEga/2.)           &
            & *(1.+OMEgal/70.))**(-1.))
   
        p_cb = -(5.-SQRT(1.+24*(1.-F_Nu)))/4.
   
        Dd_cb = (1.+(d/(1.+y_fs))**(0.7))**(-p_cb/0.7)*d**(p_cb)
   
        Dd_cbnu = ((1.-F_Nu)**(-0.7/p_cb)+(d/(1.+y_fs))**(0.7))           &
                & **(-p_cb/0.7)*d**(p_cb)
   
        END
  !*==SOUND_HORIZON_FIT.spg  processed by SPAG 6.72Dc at 19:10 on 19 Jul 2020
   
        REAL*8 FUNCTION SOUND_HORIZON_FIT(Omhh,F_baryon)
        IMPLICIT NONE
  !*--SOUND_HORIZON_FIT686
  !*** Start of declarations inserted by SPAG
        REAL*8 F_baryon , obhh , Omhh
  !*** End of declarations inserted by SPAG
  !
        obhh = F_baryon*Omhh
        SOUND_HORIZON_FIT = 44.5*LOG(9.83/Omhh)/SQRT(1.+10.0*obhh**(0.75))
   
        END
  !*==SIGMATOP.spg  processed by SPAG 6.72Dc at 19:10 on 19 Jul 2020
   
   
        REAL*8 FUNCTION SIGMATOP(Kl)
   
        IMPLICIT NONE
  !*--SIGMATOP701
  !*** Start of declarations inserted by SPAG
        REAL*8 dd0 , dd_cb , dd_cbnu , k , Kl , SCAle , TILt ,&
             & x , Z
        ! removed TF_MASTER
        INTEGER IPOwer
  !*** End of declarations inserted by SPAG
        COMMON /PASSSIGMA/ SCAle , TILt , Z , IPOwer
   
        k = EXP(Kl)
        x = SCAle*k
   
        CALL GROWTH(Z,k,dd_cb,dd_cbnu,dd0)
   
        IF ( IPOwer.EQ.0 ) THEN
           SIGMATOP = k**(3.+TILt)*(dd0*dd_cbnu*TF_MASTER(k))             &
                    & **2*(3.*(x*DCOS(x)-DSIN(x))/x**3)**2
        ELSE
   
           SIGMATOP = k**(3.+TILt)*(dd0*dd_cb*TF_MASTER(k))               &
                    & **2*(3.*(x*DCOS(x)-DSIN(x))/x**3)**2
        ENDIF
   
        END
  !*==SIGMAGAUS.spg  processed by SPAG 6.72Dc at 19:10 on 19 Jul 2020
   
        REAL*8 FUNCTION SIGMAGAUS(Kl)
   
        IMPLICIT NONE
  !*--SIGMAGAUS729
  !*** Start of declarations inserted by SPAG
        REAL*8 dd0 , dd_cb , dd_cbnu , k , Kl , SCAle , TILt ,&
             & x , Z
        ! removed: TF_MASTER
        INTEGER IPOwer
  !*** End of declarations inserted by SPAG
        COMMON /PASSSIGMA/ SCAle , TILt , Z , IPOwer
   
        k = EXP(Kl)
        x = SCAle*k
   
        CALL GROWTH(Z,k,dd_cb,dd_cbnu,dd0)
   
        IF ( IPOwer.EQ.0 ) THEN
           SIGMAGAUS = k**(3.+TILt)*(dd0*dd_cbnu*TF_MASTER(k))            &
                     & **2*EXP(-x**2.)
        ELSE
           SIGMAGAUS = k**(3.+TILt)*(dd0*dd_cb*TF_MASTER(k))              &
                     & **2*EXP(-x**2.)
        ENDIF
   
        END
  !*==SIGMAVTOP.spg  processed by SPAG 6.72Dc at 19:10 on 19 Jul 2020
   
        REAL*8 FUNCTION SIGMAVTOP(Kl)
   
        IMPLICIT NONE
  !*--SIGMAVTOP756
  !*** Start of declarations inserted by SPAG
        REAL*8 ALPha_nu , BETa_c , d , dd0 , dd_cb , dd_cbnu , f ,        &
             & F_Baryon , F_Nu , g , h , k , Kl , OMEga , OMEgal ,        &
             & omegalz , omegaz , OMHh , p_cb , q
        REAL*8 SCAle , SOUnd_horizon , THEta_cmb , TILt , x , &
             & Y_D , y_fs , Z , Z_Equality
        ! removed: TF_MASTER
        INTEGER IPOwer
  !*** End of declarations inserted by SPAG
        REAL*8 N_Nu
        COMMON /GLOBALVARIABLES/ OMHh , F_Nu , F_Baryon , N_Nu , Y_D ,    &
                               & ALPha_nu , BETa_c , SOUnd_horizon ,      &
                               & THEta_cmb , OMEga , OMEgal , Z_Equality
        COMMON /PASSSIGMA/ SCAle , TILt , Z , IPOwer
   
        k = EXP(Kl)
        x = SCAle*k
        q = k*THEta_cmb**2/OMHh
   
        CALL GROWTH(Z,k,dd_cb,dd_cbnu,dd0)
        h = SQRT(OMHh/OMEga)
        g = SQRT(OMEga*(1.+Z)**3+(1.-OMEga-OMEgal)*(1.+Z)**2+OMEgal)
        omegaz = OMEga*(1.+Z)**3/g**2
        omegalz = OMEgal/g**2
        p_cb = -(5.-SQRT(1.+24*(1.-F_Nu)))/4.
        y_fs = 17.2*F_Nu*(1.+0.488*F_Nu**(-7./6.))*(N_Nu*q/F_Nu)**2
        d = (1.+Z_Equality)/(1.+Z)                                        &
          & *5.*omegaz/2.*(omegaz**(4./7.)-omegalz+(1.+omegaz/2.)         &
          & *(1.+omegalz/70.))**(-1.)
   
        f = omegaz**(0.6) + 1./70.*omegalz*(1.+omegaz/2.)
   
        IF ( IPOwer.EQ.0 ) THEN
   
           f = f*(1.+p_cb/((1.-F_Nu)**(-0.7/p_cb)+(d/(1.+y_fs))**0.7))
           SIGMAVTOP = k**(3.+TILt)*(dd0*dd_cbnu*TF_MASTER(k))            &
                     & **2*(3.*(x*DCOS(x)-DSIN(x))/x**3)**2
        ELSE
           f = f*(1.+p_cb/(1.+(d/(1.+y_fs))**0.7))
           SIGMAVTOP = k**(3.+TILt)*(dd0*dd_cb*TF_MASTER(k))              &
                     & **2*(3.*(x*DCOS(x)-DSIN(x))/x**3)**2
        ENDIF
   
        SIGMAVTOP = (f*100.*h*g/(1.+Z)/k)**2*SIGMAVTOP
   
        END
  !*==ROMBINT.spg  processed by SPAG 6.72Dc at 19:10 on 19 Jul 2020
   
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        REAL*8 FUNCTION ROMBINT(F,A,B,Tol)
  !  Rombint returns the integral from a to b of using Romberg integration.
  !  The method converges provided that f(x) is continuous in (a,b).
  !  f must be double precision and must be declared external in the calling
  !  routine.  tol indicates the desired relative accuracy in the integral.
  !
  ! Taken from CMBFAST; Seljak & Zaldarriaga (1996)
   
   
        !PARAMETER (MAXITER=30,MAXJ=5)
        IMPLICIT NONE
  !*--ROMBINT816
  !*** Start of declarations inserted by SPAG
        REAL*8 A , B , error , fourj , g , g0 , g1 , gmax , h , Tol
        INTEGER i , j , jmax , k , MAXITER , MAXJ , nint
  !*** End of declarations inserted by SPAG
        PARAMETER (MAXITER=30,MAXJ=5) !moved down by Lukas, commented out above
        DIMENSION g(MAXJ+1)
        REAL*8 F
        EXTERNAL F
  !
        h = 0.5D0*(B-A)
        gmax = h*(F(A)+F(B))
        g(1) = gmax
        nint = 1
        error = 1.0D20
        i = 0
   100  i = i + 1
        IF ( i.GT.MAXITER .OR. (i.GT.5 .AND. ABS(error).LT.Tol) ) THEN
           ROMBINT = g0
           IF ( i.GT.MAXITER .AND. ABS(error).GT.Tol ) WRITE (*,*)        &
               & 'Rombint failed to converge; integral, error=' ,         &
              & ROMBINT , error
        ELSE
  !  Calculate next trapezoidal rule approximation to integral.
           g0 = 0.0D0
           DO k = 1 , nint
              g0 = g0 + F(A+(k+k-1)*h)
           ENDDO
           g0 = 0.5D0*g(1) + h*g0
           h = 0.5D0*h
           nint = nint + nint
           jmax = MIN(i,MAXJ)
           fourj = 1.0D0
           DO j = 1 , jmax
  !  Use Richardson extrapolation.
              fourj = 4.0D0*fourj
              g1 = g0 + (g0-g(j))/(fourj-1.0D0)
              g(j) = g0
              g0 = g1
           ENDDO
           IF ( ABS(g0).GT.Tol ) THEN
              error = 1.0D0 - gmax/g0
           ELSE
              error = gmax
           ENDIF
           gmax = g0
           g(jmax+1) = g0
           GOTO 100
        ENDIF
        END
  !*==ERFCC.spg  processed by SPAG 6.72Dc at 19:10 on 19 Jul 2020
   
   
        REAL*8 FUNCTION ERFCC(X)
  ! calculates the complimentary error function
  ! Taken from N*merical R*cipes.
        IMPLICIT NONE
  !*--ERFCC872
   
        REAL*8 X
        REAL*8 t , z
        z = ABS(X)
        t = 1./(1.+0.5*z)
        ERFCC = t*EXP(-z*z-1.26551223+                                    &
              & t*(1.00002368+t*(.37409196+t*(.09678418+                  &
              & t*(-.18628806+t*(.27886807+                               &
              & t*(-1.13520398+t*(1.48851587+t*(-.82215223+t*.17087277))))&
              & )))))
        IF ( X.LT.0. ) ERFCC = 2. - ERFCC
        END

      END module ehu_cdm
  