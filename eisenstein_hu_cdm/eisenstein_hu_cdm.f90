!Modified by Lukas Wenzl
!*==GAUSS.spg  processed by SPAG 6.72Dc at 19:48 on 22 Jul 2020
!converted from F77 to F90 with https://fortran.uk/plusfortonline.php?id=140&id1=3&id2=1&id3=1&id4=-1&id5=1

!     Date: 23 Nov 2007
!     A. Kiakotou, O. Elgaroy & O. Lahav 2007, arXiv:0709.0253
!**************************************************************************
!     Modifying the Eisenstein & Hu (1998, astro-ph/9710252) fitting
!     formula for matter power spectrum, optimized now for
!     3 degenerate massive neutrinos.
!     The driver below generates a COBE-normalized power spectrum
!**************************************************************************
 
      MODULE ehu_cdm
        
        IMPLICIT NONE

        public POWER
        public GROWTH
        contains

        !make Omeaga, H, OMEGAb, OMEGAnu inputs for POWER and then you should actually be good, need to give which k you want
        !think about precision of variables!
        !delete the calculate_pk() function and use module or test function instead
!         SUBROUTINE calculate_pk()
!         !IMPLICIT NONE
!   !*--GAUSS14
!         INTEGER nstep , ik
!         DOUBLE PRECISION OMEga , H , OMEgab , OMEganu
   
   
!         COMMON /BLOCK1/ H , OMEgab , OMEganu
!         DOUBLE PRECISION pi , h0 , rkmin , rkmax
!         DOUBLE PRECISION rk(1000) , dlnk
!         DOUBLE PRECISION p(1000)
   
   
!         OPEN (UNIT=4,FILE='pk.out',STATUS='UNKNOWN')
   
   
!   !       The distance unit is 100 km/sec (or  Mpc/h)
   
!         pi = 4.D0*ATAN(1.D0)
   
!         PRINT * , 'Assuming a Flat Universe'
!         PRINT * , 'Omega_m = Omega_c + Omega_b + Omega_nu'
!         PRINT * , 'Om_lambda = 1 - Omega_m'
   
!         PRINT * , ' ENTER Omega_m,  Omega_b , Omega_nu , h :'
!         READ * , OMEga , OMEgab , OMEganu , H
   
!   !       HUBBLE FACTOR
!         h0 = 100.D0
   
!   !       TYPE*, 'INPUT  N_STEP'
!   !       ACCEPT*, NSTEP
!         nstep = 1000
!   !        print *, ' NSTEP = ', NSTEP
   
   
!         rkmin = 1.D-4
!         rkmax = 10.D0*2.D0
   
!   !         print *, 'rkmin, rkmax=', rkmin, rkmax
   
!         dlnk = LOG(rkmax/rkmin)/DFLOAT(nstep-1)
   
!   !         print *, 'dlnk=', dlnk
!         DO ik = 1 , nstep
!            rk(ik) = rkmin*EXP(dlnk*(DFLOAT(ik)-0.5D0))
!         ENDDO
   
!         CALL POWER(nstep,rk,p)
   
!         DO ik = 1 , nstep
!            WRITE (4,99001) rk(ik) , p(ik)
!   99001    FORMAT (E30.19,E30.19)
!         ENDDO
!      END SUBROUTINE calculate_pk
  !*==POWER.spg  processed by SPAG 6.72Dc at 19:48 on 22 Jul 2020
   
   
   
        SUBROUTINE POWER(input_omega, input_omegab, input_omeganu, input_h, nsval, sigma8, anorm,Nk,Rk,Pk, z, &
            & A_s)
      !input_omega, input_omegab, input_omeganu, input_h, sigma8, anorm, Nk, Rk in units Mpc/h, Pk at one redshift, z
        IMPLICIT NONE
  !*--POWER77
        INTEGER Nk
        DOUBLE PRECISION t_0 , o_cdm , t_cmb, N_nu
        DOUBLE PRECISION omegal
        DOUBLE PRECISION Rk(Nk)
        DOUBLE PRECISION tcbnu(Nk) , Pk(Nk)
        DOUBLE PRECISION OMEga , H , OMEgab , OMEganu

        DOUBLE PRECISION, INTENT(INOUT) :: input_omega, &
             & input_omegab, input_omeganu, input_h, sigma8, &
             & z, anorm, nsval, A_s
        
        REAL(8) :: OMHh , F_Nu , F_Baryon , Y_D ,    &
             & ALPha_nu , BETa_c , SOUnd_horizon ,      &
             & THEta_cmb, Z_Equality
        COMMON /GLOBALVARIABLES/ OMHh , F_Nu , F_Baryon , N_Nu , Y_D ,    &
             & ALPha_nu , BETa_c , SOUnd_horizon ,      &
             & THEta_cmb , OMEga , OMEgal , Z_Equality
        COMMON /BLOCK1/ H , OMEgab , OMEganu, t_cmb
   
        OMEga = input_omega
        OMEgab = input_omegab
        OMEganu = input_omeganu
        H = input_h
  !     redshift
   
      !   PRINT * , ' ENTER Redshift z :'
      !   READ * , z
   
  ! number of nu flavours
  !        print *, ' ENTER Number of degenerate neutrinos:'
  !        read *, N_nu
   
  !       Try N_nu = 3.04
        !n_nu = 3.0
        n_nu = 3.046
        PRINT * , 'Assumed number of degenerate massive neutrinos:' , n_nu
   
  !     CMB temperature
        !t_0 = 2.726D0
        t_0 = 2.7255D0
        t_cmb = t_0
  !     spectral index
        !nsval = 1.D0
   
   
  !     determine Omega_lambda and Omega_cdm
  !     Note that the program is restricted to flat cosmologies
   
   
        omegal = 1.D0 - OMEga
        o_cdm = OMEga - OMEganu - OMEgab
  !
  !     write parameters to a file
  !
      !   OPEN (7,FILE='inp.dat')
      !   WRITE (7,*) OMEga , omegal
      !   WRITE (7,*) OMEganu , OMEgab
      !   WRITE (7,*) H , t_0 , n_nu
      !   WRITE (7,*) z , nsval
      !   CLOSE (7)
   
   
  !     compute transfer function
  !
        CALL TFMDM(Nk,Rk,tcbnu,anorm,sigma8)
  !
  !        compute power spectrum
  !
        CALL CONVERT(Nk,Rk,tcbnu,anorm,nsval,H,Pk, A_s, OMEga)!-OMEganutodo is it correct to substract omeganu here?
  !
        !PRINT * , 'Output k, P(k) in file pk.out'
        END
  !*==CONVERT.spg  processed by SPAG 6.72Dc at 19:48 on 22 Jul 2020
   
   
   
   
        SUBROUTINE CONVERT(N,K,Tcbnu,Anorm,N_s,H,Pk, A_s, om0)
        IMPLICIT NONE
  !*--CONVERT146
        INTEGER N
        DOUBLE PRECISION K(N) , Tcbnu(N) , Pk(N) , N_s , H , Anorm , hn
        DOUBLE PRECISION PI
        DOUBLE PRECISION k_ov_h, A_s, om0
        PARAMETER (PI=3.141592654D0)
        INTEGER i
   
        hn = H**(3.D0+N_s)
        !print *, "small scale transfer: ", Tcbnu(1) !Does go to 1 for small k as desired

        DO i = 1 , N
           !Pk(i) = 2.D0*PI*PI*hn*Anorm*Tcbnu(i)*Tcbnu(i)*K(i)**N_s
           k_ov_h = K(i)
           !A_s = 2.1e-9
           !om0 = 0.3156  !unclear!! todo
           Pk(i) =A_s*(2d0*k_ov_h**2d0*2998d0**2d0/5d0/om0)**2d0 &
             *Tcbnu(i)**2d0*(k_ov_h*H/0.05d0)**(N_s-1d0) &
             *2d0*PI**2d0/k_ov_h**3d0
            !note removed growth factor, multiplying that in the module file!
           ! Eq.(74) of Komatsu et al., ApJS, 180, 330 (2009) with kWMAP=0.05/Mpc
            ! Remember that ak is in units of h/Mpc whereas "k" in Eq.(74) is in units
            ! of 1/Mpc. Be carefull you are using the same definition of As in other modules
            ! i.e same k_pivot
        ENDDO
        END
  !*==TFMDM.spg  processed by SPAG 6.72Dc at 19:48 on 22 Jul 2020
   
   
   
   
   
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
  !	    numk    --  number of k's per decade for evaluation
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
  !				 (taken from Numerical Recipes)
  !
  !	  Arguments:
  !
  !	      x       -- real*8
  !
   
        SUBROUTINE TFMDM(N2df,K2df,Tfunc,Anorm,Sigma_8)
        IMPLICIT NONE
  !*--TFMDM341
  !*** Start of declarations inserted by SPAG
        REAL*8 ALPha_nu , Anorm , BETa_c , dd0 , dd_cb , dd_cbnu ,        &
             & F_Baryon , F_Nu , k , kmax , kmin , OMEga , OMEgal , OMHh ,&
             & SCAle , Sigma_8 , SOUnd_horizon 
        !removed Lukas: ROMBINT, SIGMATOP, TF_MASTER
        REAL*8 THEta_cmb , TILt , tol , t_cmb , t_master , Y_D , Z ,      &
             & Z_Equality
        INTEGER i , IPOwer , iz , numk
  !*** End of declarations inserted by SPAG
        INTEGER N2df
        DOUBLE PRECISION K2df(N2df) , Tfunc(N2df)
        !removed: EXTERNAL ROMBINT , SIGMATOP , SIGMAGAUS , SIGMAVTOP
        REAL*8 N_Nu
        COMMON /GLOBALVARIABLES/ OMHh , F_Nu , F_Baryon , N_Nu , Y_D ,    &
                               & ALPha_nu , BETa_c , SOUnd_horizon ,      &
                               & THEta_cmb , OMEga , OMEgal , Z_Equality
        COMMON /PASSSIGMA/ SCAle , TILt , Z , IPOwer
        DOUBLE PRECISION omeganu , omegab , nsval , h
        COMMON /BLOCK1/ H , OMEgab , OMEganu, t_cmb
   
  !      Input Cosmological Parameters
   
      !   OPEN (7,FILE='inp.dat')
      !   READ (7,*) OMEga , OMEgal
      !   READ (7,*) omeganu , omegab
      !   READ (7,*) h , t_cmb , N_Nu
      !   READ (7,*) Z , nsval
      !   CLOSE (7)
   
        
        PRINT * , 'nsval:' , nsval
   
        !TILt = nsval
        IPOwer = 1
        !kmax = 30
        !numk = 30
   
   
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
   
   
        CALL TFSET_PARAMETERS
   
  !	sound_horizon=sound_horizon_fit(omhh,f_baryon)
   
        !kmin = 0.0001
        !numk = numk*LOG10(kmax/kmin)
   
        DO i = 1 , N2df
           k = K2df(i)
           CALL GROWTH(Z,k*h,dd_cb,dd_cbnu,dd0)
           t_master = TF_MASTER(k*h)
           Tfunc(i) = t_master
           !*dd_cbnu*dd0 removed the growth factor here and add it 
           ! note that dd0 is normalized for COBE not to for when using a primodial power spectrum
           ! it has dd0 -> 1 for z=0 small k
           ! when using priomodial power spectrum you want dd0 ->1 for high z and lower k (here it is about 1.25 for lambda cdm)
        ENDDO
   
      !everything below here can be commented out: todo
  !	Power Spectrum Statistics
   
        iz = 1
        TILt = 1
   
        IF ( ABS(1-OMEga-OMEgal).LT.1E-5 ) IPOwer = 1
   
        IF ( IPOwer.EQ.0 ) THEN
           Anorm = 1.94E-5*OMEga**(-0.785-0.05*LOG(OMEga))                &
                 & *EXP(1.*(TILt-1.)+1.97*(TILt-1.)**2)
   
           Anorm = Anorm**2*(2997.9/h)**(3.+TILt)
   
        ELSE
           Anorm = 1.94E-5*OMEga**(-0.785-0.05*LOG(OMEga))                &
                 & *EXP(-0.95*(TILt-1.)-0.169*(TILt-1.)**2)
   
           Anorm = Anorm**2*(2997.9/h)**(3.+TILt)
        ENDIF
   
  !
   
        IF ( OMEgal.EQ.0 ) THEN
   
           IF ( TILt.LT.1 ) THEN
              IPOwer = 1
           ELSE
              IPOwer = 1
   
           ENDIF
   
           IF ( IPOwer.EQ.0 ) THEN
              Anorm = 1.95E-5*OMEga**(-0.35-0.19*LOG(OMEga)+0.15*(TILt-1))&
                    & *EXP(1.02*(TILt-1.)+1.70*(TILt-1.)**2)
   
              Anorm = Anorm**2*(2997.9/h)**(3.+TILt)
           ELSE
              Anorm = 1.95E-5*OMEga**(-0.35-0.19*LOG(OMEga)-0.17*(TILt-1))&
                    & *EXP(-1.*(TILt-1.)-0.14*(TILt-1.)**2)
              Anorm = Anorm**2*(2997.9/h)**(3.+TILt)
           ENDIF
   
        ELSE
   
           Anorm = 2.422 - 1.166*EXP(OMEga) + 0.800*EXP(OMEgal)           &
                 & + 3.780*OMEga - 2.267*OMEga*EXP(OMEgal)                &
                 & + 0.487*OMEga**2 + 0.561*OMEgal +                      &
                 & 3.392*OMEgal*EXP(OMEga) - 8.568*OMEga*OMEgal +         &
                 & 1.080*OMEgal**2
           Anorm = 1.E-5*Anorm
           Anorm = Anorm**2*(2997.9/h)**(3.+TILt)
  !
        ENDIF
  !
   
        tol = 1.E-6
  !
  !  we save some time by not computing the following statistics
  !  if you want to include them, just uncomment the lines below
  !  OE 12/12-01
  !
  !	Sigma_8: Mass fluctuations in a top hat of 8./h Mpc at
  !	         redshift z=0; total mass: CDM+baryon+neutrino
   
        SCAle = 8./h
        IF ( iz.EQ.0 ) Z = 0
        IPOwer = 0

        Sigma_8 = SQRT                                                    &
                & (Anorm*(ROMBINT(SIGMATOP,DLOG(0.001/SCAle),DLOG(0.1D0/  &
                & SCAle),tol)+ROMBINT(SIGMATOP,DLOG(0.1D0/SCAle),         &
                & DLOG(1D0/SCAle),tol)                                    &
                & +ROMBINT(SIGMATOP,DLOG(1D0/SCAle),DLOG(10.D0/SCAle),tol)&
                & +ROMBINT(SIGMATOP,DLOG(10.D0/SCAle),DLOG(100.D0/SCAle), &
                & tol)))
   
        PRINT * , 'COBE normalized z, sigma_8: ' , Z , Sigma_8
  !	write(6,51) sigma_8,z
  !51	FORMAT('sigma_8:   ',E12.5,
  !     *      '                                  z=',F5.2)
  !
  !
  !
  !	Sigma_50: Mass fluctuations in a top hat of 50./h Mpc at
  !	         redshift z=0; total mass: CDM+baryon+neutrino
  !
  !	scale = 50./h
  !	if (iz.eq.0) z=0
  !	ipower=0
  !
  !	sigma_50=
  !     *	    sqrt(anorm*(
  !     *	    rombint(sigmatop,dlog(0.001/scale),dlog(0.1d0/scale),tol)+
  !     *	    rombint(sigmatop,dlog(0.1d0/scale),dlog(1d0/scale),tol)+
  !     *	    rombint(sigmatop,dlog(1d0/scale),dlog(10.d0/scale),tol)+
  !     *	    rombint(sigmatop,dlog(10.d0/scale),dlog(100.d0/scale),tol)))
  !
  !	write(6,52) sigma_50,sigma_50/sigma_8,z
  !52	FORMAT('sigma_50:  ',E12.5,'    sigma_50/sigma_8:',E12.5,
  !     *		' z=',F5.2)
  !
  !
  !
  !	Sigma_DLA:  Mass fluctuations in a gaussian of
  !	            scale corresponding to 50km/s at
  !	            redshift z=4; total mass: CDM+baryon+neutrino mass
  !
  !
  !	circular velocity cutoff
  !	vcirc = 50.
  !	if (iz.eq.0) z=4
  !	scale = 0.3834*(vcirc)/(h*100.)*(omega*(1.+z)**(3.)
  !     *		+(1.-omega-omegal)*(1.+z)**2+omegal)**(-1./6.)
  !	ipower=0
  !
  !     	sigmadla=sqrt(anorm*(
  !     *	rombint(sigmagaus,dlog(0.001/scale),dlog(0.1d0/scale),tol)+
  !     *	rombint(sigmagaus,dlog(0.1d0/scale),dlog(1d0/scale),tol)+
  !     *	rombint(sigmagaus,dlog(1d0/scale),dlog(10.d0/scale),tol)+
  !     *	rombint(sigmagaus,dlog(10.d0/scale),dlog(20.d0/scale),tol)))
  !
  !	write(6,54) sigmadla,scale*h,z
  !54	FORMAT('sigma_DLA: ',E12.5,'    at',E12.5,
  !     *		' h^-1 Mpc,      z=',F5.2)
  !
  !	Upper limit on Omega_gas via Press-Schecter
  !	Press-Schechter threshold
  !
  !	deltac = 1.33
  !	omegagas = omega*f_baryon*erfcc(deltac/(sqrt(2.)*sigmadla))
  !
  !	write(6,55) omegagas,deltac,z
  !55	FORMAT('Omega_gas:<',E12.5,'    at ',F5.3,
  !     *		' gaus. deltac,        z=',F5.2)
  !
  !	Sigma_Lya:  Mass fluctuations in a gaussian of
  !	            1/(24.04 sqrt(0mega h**2)) Mpc at
  !	            redshift z=3; CDM+baryon mass
  !
  !	scale = 1./(24.04*sqrt(omhh))
  !
  !	if (iz.eq.0) z=3
  !	ipower=1
  !
  !        sigma_Lya=
  !     *	    sqrt(anorm*(
  !     *	    rombint(sigmagaus,dlog(0.001/scale),dlog(0.1d0/scale),tol)+
  !     *	    rombint(sigmagaus,dlog(0.1d0/scale),dlog(1d0/scale),tol)+
  !     *	    rombint(sigmagaus,dlog(1d0/scale),dlog(10.d0/scale),tol)+
  !     *	    rombint(sigmagaus,dlog(10.d0/scale),dlog(20.d0/scale),tol)))
  !	write(6,53) sigma_Lya,scale*h,z
  !53	FORMAT('sigma_Lya: ',E12.5,'    at',E12.5,
  !     *		' h^-1 Mpc,      z=',F5.2)
  !
  !
  !	Sigmav_40:  Rms velocity in top hat of 100h^{-1} Mpc, z=0
  !		     default; CDM+baryon
  !
  !	scale = 40./h
  !	if (iz.eq.0) z=0
  !	ipower=1
  !
  !        sigmav_40=
  !     *	    sqrt(anorm*(
  !     *	    rombint(sigmavtop,dlog(0.001/scale),dlog(0.1d0/scale),tol)+
  !     *	    rombint(sigmavtop,dlog(0.1d0/scale),dlog(1d0/scale),tol)+
  !     *	    rombint(sigmavtop,dlog(1d0/scale),dlog(10.d0/scale),tol)+
  !     *	    rombint(sigmavtop,dlog(10.d0/scale),dlog(100.d0/scale),tol)))
  !	write(6,60) sigmav_40,scale*h,z
  !60	FORMAT('sigma_v:    ',F6.2,' km/s    at ',F6.2,
  !     *		' h^-1 Mpc,           z=',F5.2)c
  !
  !	write(6,*) '     ----------------------------------------------'
  !	write(6,*) 'Extras:'
  !	write(6,*) '  '
  !
  !	if (iz.eq.0) z = 0.
  !
  !	  yp = primordial Helium mass fraction
  !	  yp = 0.24
  !	  xe = ionization fraction 1                 = fully ionized
  !			           (1-yp)/(1-yp/2)   = HeI
  !			           (1-3yp/4)/(1-yp/2)= HeII
  !	  xe = (1.-yp)/(1.-yp/2)
  !	  neutral helium
  !	  tau= 4.61E-2*(1.-yp/2.)*xe*omegab*h/omega**2
  !
  !	if (omegal.eq.0) then
  !	  omegae= omega-1.d-6
  !	  arg = 2.d0*(1.d0-omegae)/omegae/(1.+z) + 1.d0
  !     	  age = 9.7778/h*(
  !     *	        sqrt(1.+omegae*z)/(1.-omegae)/(1.+z) -
  !     *          omegae/2./(1.-omegae)**(1.5)*dlog(
  !     *			arg + dsqrt(arg**2 - 1.d0)))
  !	  write(6,56) z,age
  !
  !	  tau= tau*(2.-3.*omega+sqrt(1.+omega*z)*(omega*z+3.*omega-2.))
  !
  !	  if (z.gt.0) write(6,57) z,tau
  !
  !	else
  !
  !	if (abs(1.-omega-omegal).lt.1E-10) then
  !	  g = sqrt(omega*(1.+z)**3+1.-omega)
  !	  omegaz = omega*(1.+z)**3/g**2
  !
  !	  age = 2./3.*9.7778/h/g/sqrt(1.-omegaz)
  !     *  	*log((1.+sqrt(1.-omegaz))/sqrt(omegaz))
  !	  write(6,56) z,age
  !
  !	  tau= tau*omega*(sqrt(1.-omega+omega*(1.+z)**3) -1.)
  !
  !	  if (z.gt.0) write(6,57) z,tau
  !
  !	end if
  !
  !	end if
  !
  !56 	FORMAT('Age(z=',F5.2,'):    ',F6.2,' Gyrs')
  !57 	FORMAT('tau(z=',F5.2,'):    ',F9.5,' (HeI)')
  !
  !c	if (z.gt.0) write(6,58) z,1./DD0
  !
  !58	FORMAT('Growth(z=',F5.2,'-0):',F5.2)
  !
  !
  !	write(6,59) sound_horizon*h
  !59	FORMAT('sound_horizon:  ',F7.2,' h^-1 Mpc')
  !
  !
  !	call Growth(z,1.d5,dum1,dum2,dum3)
  !	write(6,61) alpha_nu*dum1,z
  !61	FORMAT('suppression:    ',F7.5,' z=',F5.2)
        END
  !*==TFSET_PARAMETERS.spg  processed by SPAG 6.72Dc at 19:48 on 22 Jul 2020
   
   
        SUBROUTINE TFSET_PARAMETERS
        IMPLICIT NONE
  !*--TFSET_PARAMETERS648
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
        IF ( f_cb.EQ.0.D0 ) f_cb = f_cb + 1.D-10
        IF ( 4.*p_cb+5.0.EQ.0.0 ) p_cb = p_cb + 1.D-10
        IF ( 1.-.193*SQRT(F_Nu)+.169*F_Nu.EQ.0.0 ) F_Nu = F_Nu + 1.D-10
        IF ( 4.*p_c+3.0.EQ.0.0 ) p_c = p_c + 1.D-10
        IF ( 4.*p_cb+7.0.EQ.0.0 ) p_cb = p_cb + 1.D-10
        ALPha_nu = (f_c/f_cb)*(2.*(p_c+p_cb)+5.)/(4.*p_cb+5.)
        ALPha_nu = ALPha_nu*(1.-0.553*f_nub+0.126*f_nub**3)
        ALPha_nu = ALPha_nu/(1.-0.193*SQRT(F_Nu)+0.169*F_Nu)
        ALPha_nu = ALPha_nu*(1.+Y_D)**(p_c-p_cb)
        ALPha_nu = ALPha_nu*(1.+(p_cb-p_c)                                &
                 & /2.*(1.+1./(4.*p_c+3.)/(4.*p_cb+7.))/(1.+Y_D))
        BETa_c = 1.D0/(1.D0-0.949D0*f_nub)
   
   
        END
  !*==TF_MASTER.spg  processed by SPAG 6.72Dc at 19:48 on 22 Jul 2020
   
   
        REAL*8 FUNCTION TF_MASTER(K)
   
        IMPLICIT NONE
  !*--TF_MASTER708
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
        IF ( ALPha_nu.LT.0.D0 ) ALPha_nu = 0.D0     ! fix put in 13/03/02 OE
        gamma_eff = (SQRT(ALPha_nu)+(1.-SQRT(ALPha_nu))                   &
                  & /(1.+(0.43*K*SOUnd_horizon)**4))
        q_eff = q/gamma_eff
        TF_MASTER = DLOG(DEXP(1.D0)+1.84*BETa_c*SQRT(ALPha_nu)*q_eff)
        TF_MASTER = TF_MASTER/                                            &
                  & (TF_MASTER+q_eff**2*(14.4+325./(1.+60.5*q_eff**1.11)))
  !     AK modified
        q_nu = 4.9*q*SQRT(N_Nu/F_Nu)
  !     AK modified
        TF_MASTER = TF_MASTER*(1.+(1.2*F_Nu**(0.68)*N_Nu**(0.3+0.6*F_Nu)) &
                  & /(q_nu**(-1.3)+q_nu**(0.45)))
   
   
        END
  !*==GROWTH.spg  processed by SPAG 6.72Dc at 19:48 on 22 Jul 2020
   
        SUBROUTINE GROWTH(Z,K,Dd_cb,Dd_cbnu,Dd0)
   
        IMPLICIT NONE
  !*--GROWTH740
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
  !*==SOUND_HORIZON_FIT.spg  processed by SPAG 6.72Dc at 19:48 on 22 Jul 2020
   
        REAL*8 FUNCTION SOUND_HORIZON_FIT(Omhh,F_baryon)
        IMPLICIT NONE
  !*--SOUND_HORIZON_FIT780
  !*** Start of declarations inserted by SPAG
        REAL*8 F_baryon , obhh , Omhh
  !*** End of declarations inserted by SPAG
  !
        obhh = F_baryon*Omhh
        SOUND_HORIZON_FIT = 44.5*LOG(9.83/Omhh)/SQRT(1.+10.0*obhh**(0.75))
   
        END
  !*==SIGMATOP.spg  processed by SPAG 6.72Dc at 19:48 on 22 Jul 2020
   
   
        REAL*8 FUNCTION SIGMATOP(Kl)
   
        IMPLICIT NONE
  !*--SIGMATOP795
  !*** Start of declarations inserted by SPAG
        REAL*8 dd0 , dd_cb , dd_cbnu , k , Kl , SCAle , TILt ,&
             & x , Z
        !removed Lukas: TF_MASTER
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
  !*==SIGMAGAUS.spg  processed by SPAG 6.72Dc at 19:48 on 22 Jul 2020
   
        REAL*8 FUNCTION SIGMAGAUS(Kl)
   
        IMPLICIT NONE
  !*--SIGMAGAUS823
  !*** Start of declarations inserted by SPAG
        REAL*8 dd0 , dd_cb , dd_cbnu , k , Kl , SCAle , TILt ,&
             & x , Z
      ! removed Lukas: TF_MASTER
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
  !*==SIGMAVTOP.spg  processed by SPAG 6.72Dc at 19:48 on 22 Jul 2020
   
        REAL*8 FUNCTION SIGMAVTOP(Kl)
   
        IMPLICIT NONE
  !*--SIGMAVTOP850
  !*** Start of declarations inserted by SPAG
        REAL*8 ALPha_nu , BETa_c , d , dd0 , dd_cb , dd_cbnu , f ,        &
             & F_Baryon , F_Nu , g , h , k , Kl , OMEga , OMEgal ,        &
             & omegalz , omegaz , OMHh , p_cb , q
        REAL*8 SCAle , SOUnd_horizon , THEta_cmb , TILt , x , &
             & Y_D , y_fs , Z , Z_Equality
        ! removed Lukas: 
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
  !*==ROMBINT.spg  processed by SPAG 6.72Dc at 19:48 on 22 Jul 2020
   
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
  !*--ROMBINT910
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
  !*==ERFCC.spg  processed by SPAG 6.72Dc at 19:48 on 22 Jul 2020
   
   
        REAL*8 FUNCTION ERFCC(X)
  ! calculates the complimentary error function
  ! Taken from Numerical Recipes.
        IMPLICIT NONE
  !*--ERFCC966
   
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
   
   
  