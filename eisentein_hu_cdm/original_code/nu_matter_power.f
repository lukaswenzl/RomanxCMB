
c     Date: 23 Nov 2007
c     A. Kiakotou, O. Elgaroy & O. Lahav 2007, arXiv:0709.0253
C**************************************************************************
C     Modifying the Eisenstein & Hu (1998, astro-ph/9710252) fitting
C     formula for matter power spectrum, optimized now for 
c     3 degenerate massive neutrinos.  
c     The driver below generates a COBE-normalized power spectrum
C**************************************************************************

        PROGRAM GAUSS
        implicit none 
        integer           nstep,IR,ik, np  
        double precision  OMEGA,H,OMEGAB,OMEGANU


        COMMON /BLOCK1/OMEGA,H,OMEGAB,OMEGANU 
        double precision pi, H0,sig8,RSTAR,rkmin,rkmax   
        double precision RK(1000),dlnk,sigma,sigma2,uk, ar, rka, prk
        double precision P(1000) 
        

        OPEN(UNIT=4,FILE='pk.out', STATUS='UNKNOWN')


C       The distance unit is 100 km/sec (or  Mpc/h)

        PI = 4.d0*ATAN(1.d0)

        print *, 'Assuming a Flat Universe'
        print *, 'Omega_m = Omega_c + Omega_b + Omega_nu'
        print *, 'Om_lambda = 1 - Omega_m'

        print *, ' ENTER Omega_m,  Omega_b , Omega_nu , h :'
        read *, OMEGA, OMEGAB,OMEGANU,H

C       HUBBLE FACTOR
        H0 = 100.d0

C       TYPE*, 'INPUT  N_STEP'
C       ACCEPT*, NSTEP
        NSTEP = 1000
c        print *, ' NSTEP = ', NSTEP


          RKMIN = 1.d-4
          RKMAX= 10.d0*2.d0

c         print *, 'rkmin, rkmax=', rkmin, rkmax 

        DLNK=LOG(RKMAX/RKMIN)/DFLOAT(NSTEP-1)

c         print *, 'dlnk=', dlnk
        do ik = 1, nstep 
           rk(ik) = rkmin*exp(dlnk*(dfloat(ik)-0.5d0))
        end do

        call power(nstep,rk,P)

           DO IK=1,NSTEP              
             write(4,80) rk(ik), P(ik)
 80          format(E30.19,E30.19)
           END DO 



        
        END



      subroutine power(Nk,rk,Pk)
      implicit none 
      integer             Nk
      double precision    T_0, O_cdm, N_nu, z
      double precision    nsval 
      double precision    OMEGAL 
      double precision    rk(Nk)
      double precision    Tcbnu(Nk), Pk(Nk),anorm
      double precision    sigma8
      double precision    OMEGA,H,OMEGAB, OMEGANU 
      COMMON /BLOCK1/OMEGA,H,OMEGAB,OMEGANU 



c     redshift

      print *, ' ENTER Redshift z :'
      read *, z

c number of nu flavours
c        print *, ' ENTER Number of degenerate neutrinos:'
c        read *, N_nu

c       Try N_nu = 3.04
        N_nu =3.0
        print *, 'Assumed number of degenerate massive neutrinos:', N_nu

c     CMB temperature
      T_0 = 2.726d0 
c     spectral index
      nsval = 1.d0 


c     determine Omega_lambda and Omega_cdm
c     Note that the program is restricted to flat cosmologies 


      OMEGAL  = 1.d0 - OMEGA 
      O_cdm   = OMEGA - OMEGANU - OMEGAB 
c
c     write parameters to a file 
c
      open(7,file='inp.dat')
      write(7,*)OMEGA, OMEGAL
      write(7,*)OMEGANU, OMEGAB  
      write(7,*)H, T_0, N_nu 
      write(7,*)z, nsval     
      close(7)
      

c     compute transfer function
c
      call tfmdm(Nk,rk,Tcbnu,anorm,sigma8)
c


c        compute power spectrum
c
      call convert(Nk,rk,Tcbnu,anorm,nsval,h,Pk)
c
      PRINT *, 'Output k, P(k) in file pk.out'

      end   




      subroutine convert(N,k,Tcbnu,anorm,n_s,h,Pk)
      implicit none 
      integer  N 
      double precision  k(N),Tcbnu(N),Pk(N),n_s,h,anorm, hn 
      double precision pi 
      parameter (pi = 3.141592654d0) 
      integer  i 

      hn = h**(3.d0+n_s)
      do i = 1, N 
         Pk(i)=2.d0*pi*pi*hn*anorm*Tcbnu(i)*Tcbnu(i)*k(i)**n_s
      end do
      return 
      end





c	Driver and Transfer Function / Power Spectrum subroutines 
c	       for Eisenstein & Hu astro-ph/97XXXXX.
c		3/5/98: typo in lambda age fixed 3/5/98
c		9/1/98: typo in normalization of open gravity wave models fixed 
c		9/22/98:bug in sigmavtop (affecting for massive neutrinos) fixed 
c
c	  Here is a sample input file that may be cut out and used
c
cX------------------------------------------------------------------------------
c 		1.0, 0.0, 0.0, 0.05	Omega_0,Omega_L,Omega_nu,Omega_b
c 		0.5, 2.726, 1.		h,T_cmb,N_nu
c 		0			z
c 		30,30			kmax,numk
c  		1			iz
c 	        1.0			n
c 		0			ipower
cX------------------------------------------------------------------------------
c 
c       DRIVER Program Power
c
c	  INPUT:
c
c	    omega   --	matter density in critical units (no lambda)
c	    omega_l --  cosmological constant in critical units
c	    omega_nu--  massive neutrino density in critical units
c	    omega_b --  baryon density in critical units
c	    h       --  Hubble constant in 100 km/s/Mpc
c	    T_cmb   --  CMB Temperature 
c	    N_nu    --  number of degenerate massive neutrinos
c	    z	    --  redshift at which to evaluate transfer function
c	    kmax    --  maximum k (h Mpc^{-1}) for transfer function tabulation
c	    numk    --  number of k's per decade for evaluation
c
c
c	    tilt    --  n (Power spectrum tilt)
c	    ipower  --  if (n<1) and omega or omegal=0, option to
c			normalize (0) with power law inflationary
c			gravity waves (1) no gravity wavesr; else
c			no gravity waves
c	    iz      --  (0) evaluate quantities at default redshifts
c	                (1) evaluate quantities at specified redshift
c
c
c	  COMMON/GLOBALVARIABLES:  Variables in transfer function
c				   evaluation
c 	    (all real*8)	
c 
c       *   omhh         --  Omega*h*h 
c	*   f_nu         --  Omega_nu/Omega
c	*   f_baryon     --  Omega_b/Omega
c	*   N_nu         --  number of degenerate massive neutrinos
c	    y_d          --  (1+z_eq)/(1+z_d)
c	    alpha_nu     --  small scale suppression
c	    beta_c       --  logarithmic correction
c	    sound_horizon--  sound horizon in Mpc
c	    theta_cmb    --  T_cmb/2.7K
c	*   omega        --  matter density in critical units 
c	*   omegal       --  cosmological constant in critical units
c	    z_equality   --  redshift of matter-radiation equality
c
c	* = must be defined before call to TFset_parameters
c
c	COMMON/PASSSIGMA: variables needed for mass fluctuation 
c	                  evalutation
c
c	    scale  -- filter scale in Mpc
c	    tilt   -- n, power spectrum tilt
c	    z      -- redshift of evaluation
c	    ipower -- switch (0) P_cbnu (1) P_cb
c
c	OUTPUT
c
c	    delta_H   -- COBE normalization for horizon-scale 
c		         fluctuation amplitude
c	  Examples of Top-hat and Gaussian window mass fluctuations
c
c	    sigma_8   -- delta M/M on 8  h^-1 Mpc (using P_cbnu)
c	    sigma_50  -- delta M/M on 50 h^-1 Mpc (using P_cbnu)
c	    sigma_Lya -- delta M/M on 1/(24 sqrt(0mega*h*h)) Mpc scale 
c			 (using P_cb) at z=3
c	    sigma_DLA -- delta M/M on 50km/s scale 
c			 (using P_cbnu) at z=4
c	    omega_gas -- Upper limit on Omega_gas (z=4) via
c			 Press-Schechter with gaussian filter at
c			 50 km/s deltac=1.33
c	
c	  Extra Freebies
c
c	    sigmav_40  -- velocity on 40 h^-1 Mpc tophat spheres
c			 (using P_cb)
c
c	    age        -- age in Gyr at z
c	    tau	       -- Thomson optical depth for full ionization to z
c	    growth     -- D1(z=0)/D1(z) growth from z to 0
c	    suppression-- alpha_nu*D_1**(-p_cb) suppression in CDM baryon transfer function
c			  as k->Infty at z
c
c      TRANSFER FUNCTION SUBROUTINES/FUNCTIONS
c	
c	TF_setparameters [subroutine]
c	
c	  Uses COMMON/GLOBALVARIABLES (above) 
c	  Sets unstarred variables in GLOBALVARIABLES list above
c
c	TF_master(k) [real*8 function]
c
c	  Argument: k wavenumber in Mpc^{-1} (IMPORTANT: no h)
c	  Uses COMMON/GLOBALVARIABLES (above) 
c	
c	Growth [subroutine]
c
c	  Arguments:
c
c	      z       -- redshift to evaluate growth
c	      k       -- wavenumber in Mpc^{-1} (IMPORTANT: no h)
c	      DD_cb   -- on exit set to D_cb(z)/D_1(z=0)
c	      DD_cbnu -- on exit set to D_cbnu(z)/D_1(z=0)
c	      DD0     -- on exit set to D_1(z=0)
c
c	  Uses COMMON/GLOBALVARIABLES (above) 
c
c	sound_horizon [real*8 function]: approximate sound horizon
c	  [a convenient approximation that does not require the
c	   parameter set up of TF_setparameters; unnecessary if
c	   TF_setparameters is called.
c	
c	  Arguments:
c
c	      omhh     -- Omega*h*h
c	      f_baryon -- Omega_b/Omega
c
c       sigmatop [real*8 function]: evaluates integrand of sigma for 
c 	  a top hat filter (unnormalized) 
c
c	  Arguments:
c
c	      kl       -- ln(k Mpc^{-1})  [real*8]		
c		
c	  Uses COMMON/PASSSIGMA
c
c       sigmagaus [real*8 function]: evaluates integrand of sigma for 
c 	  a gaussian filter (unnormalized) 
c
c	  Uses COMMON/PASSSIGMA
c
c	sigmavtop [real*8 function]: evaluates integrand of sigmav for
c	  a top hat filter (unnormalized)
c
c	  Uses COMMON/PASSSIGMA
c	       COMMON/GLOBALVARIABLES
c
c	  Arguments:
c
c	      kl       -- ln(k Mpc^{-1})  [real*8]		
c		
c	  Uses COMMON/PASSSIGMA
c
c	rhombint [real*8 function]: Rhomberg integration 
c				    (taken from CMBfast)
c
c	  Arguments:
c	  
c	      f       -- integrand externally defined 
c	      a       -- lower limit of integration [real*8]
c	      b       -- upper limit of integration [real*8]
c	      tol     -- tolerance [real*8]
c
c	erfcc [real*8 function]: complementary error function
c				 (taken from Numerical Recipes)
c
c	  Arguments:
c	
c	      x       -- real*8 
c

      subroutine TFmdm(N2dF,k2dF,Tfunc,anorm,sigma_8) 
	implicit real*8 (a-h,k,o-z)
        integer N2dF 
        double precision k2dF(N2dF), Tfunc(N2dF)  
	external rombint,sigmatop,sigmagaus,sigmavtop
 	real*8 N_nu
        common/GLOBALVARIABLES/omhh,f_nu,f_baryon,N_nu,y_d,alpha_nu,
     *     beta_c,sound_horizon,theta_cmb,omega,omegal,z_equality
	common/PASSSIGMA/scale,tilt,z,ipower
        double precision omeganu, omegab, nsval,h

c      Input Cosmological Parameters

        open(7,file='inp.dat')
        read(7,*)OMEGA, OMEGAL
        read(7,*)OMEGANU, OMEGAB  
        read(7,*)H, T_cmb, N_nu 
        read(7,*)z, nsval 
        close(7)



        tilt = nsval 
        ipower = 1  
        kmax = 30 
        numk = 30 


c	Translate Parameters into forms GLOBALVARIABLES form

	f_nu = omeganu/omega
	f_baryon = omegab/omega
	if (f_nu.eq.0) f_nu=1E-10
	if (f_baryon.eq.0) f_baryon=1E-10
	if (T_cmb.le.0) T_cmb=2.726
	if (N_nu.lt.0) N_nu=0.
	if (h.gt.10) write(6,*) 'WARNING: H_0/100km/s/Mpc needed'
	theta_cmb=T_cmb/2.7
	omhh = omega*h*h



        call TFset_parameters

c	sound_horizon=sound_horizon_fit(omhh,f_baryon) 

        kmin = 0.0001
        numk = numk*log10(kmax/kmin)


        do i = 1, N2dF           
           k = k2dF(i)  
	 call Growth(z,k*h,DD_cb,DD_cbnu,DD0)
	 T_master = TF_master(k*h)
         Tfunc(i) = T_master * DD_cbnu*DD0 
	end do


c	Power Spectrum Statistics

        iz = 1
        tilt =1

	if (abs(1-omega-omegal).lt.1E-5) then
           ipower=1
	 end if

	 if (ipower.eq.0) then 
            anorm = 1.94E-5*omega**(-0.785-0.05*log(omega))
     *			*exp(1.*(tilt-1.)+1.97*(tilt-1.)**2)

            anorm = anorm**2*(2997.9/h)**(3.+tilt)

	 else
            anorm = 1.94E-5*omega**(-0.785-0.05*log(omega))
     *                  *exp(-0.95*(tilt-1.)-0.169*(tilt-1.)**2)

            anorm = anorm**2*(2997.9/h)**(3.+tilt)
	 end if

	

         if (omegal.eq.0) then

	 if (tilt.lt.1) then 
            ipower = 1 
	 else 
	  ipower=1

	 end if

	 if (ipower.eq.0) then 
            anorm = 1.95E-5*omega**(-0.35-0.19*log(omega)+0.15*(tilt-1))
     *			*exp(1.02*(tilt-1.)+1.70*(tilt-1.)**2)

            anorm = anorm**2*(2997.9/h)**(3.+tilt)
	 else
            anorm = 1.95E-5*omega**(-0.35-0.19*log(omega)-0.17*(tilt-1))
     *			*exp(-1.*(tilt-1.)-0.14*(tilt-1.)**2)
            anorm = anorm**2*(2997.9/h)**(3.+tilt)
	 end if

	else

	anorm = 2.422-1.166*exp(omega)+0.800*exp(omegal)+3.780*omega
     *		-2.267*omega*exp(omegal)+0.487*omega**2+0.561*omegal
     *		+3.392*omegal*exp(omega)-8.568*omega*omegal
     *		+1.080*omegal**2
	anorm = 1.E-5*anorm
            anorm = anorm**2*(2997.9/h)**(3.+tilt)
	    
	end if 
	 

        
	tol = 1.E-6
c
c  we save some time by not computing the following statistics
c  if you want to include them, just uncomment the lines below 
c  OE 12/12-01 
c
c	Sigma_8: Mass fluctuations in a top hat of 8./h Mpc at
c	         redshift z=0; total mass: CDM+baryon+neutrino

	scale = 8./h
	if (iz.eq.0) z=0
	ipower=0

	sigma_8=
     *	    sqrt(anorm*(
     *	    rombint(sigmatop,dlog(0.001/scale),dlog(0.1d0/scale),tol)+
     *	    rombint(sigmatop,dlog(0.1d0/scale),dlog(1d0/scale),tol)+
     *	    rombint(sigmatop,dlog(1d0/scale),dlog(10.d0/scale),tol)+
     *	    rombint(sigmatop,dlog(10.d0/scale),dlog(100.d0/scale),tol)))

        print *, 'COBE normalized z, sigma_8: ', z, sigma_8
c	write(6,51) sigma_8,z
c51	FORMAT('sigma_8:   ',E12.5,
c     *      '                                  z=',F5.2)
c
c
c
c	Sigma_50: Mass fluctuations in a top hat of 50./h Mpc at
c	         redshift z=0; total mass: CDM+baryon+neutrino
c
c	scale = 50./h
c	if (iz.eq.0) z=0
c	ipower=0
c
c	sigma_50=
c     *	    sqrt(anorm*(
c     *	    rombint(sigmatop,dlog(0.001/scale),dlog(0.1d0/scale),tol)+
c     *	    rombint(sigmatop,dlog(0.1d0/scale),dlog(1d0/scale),tol)+
c     *	    rombint(sigmatop,dlog(1d0/scale),dlog(10.d0/scale),tol)+
c     *	    rombint(sigmatop,dlog(10.d0/scale),dlog(100.d0/scale),tol)))
c
c	write(6,52) sigma_50,sigma_50/sigma_8,z
c52	FORMAT('sigma_50:  ',E12.5,'    sigma_50/sigma_8:',E12.5,
c     *		' z=',F5.2)
c
c
c
c	Sigma_DLA:  Mass fluctuations in a gaussian of 
c	            scale corresponding to 50km/s at 
c	            redshift z=4; total mass: CDM+baryon+neutrino mass 
c
c
c	circular velocity cutoff
c	vcirc = 50.
c	if (iz.eq.0) z=4
c	scale = 0.3834*(vcirc)/(h*100.)*(omega*(1.+z)**(3.) 
c     *		+(1.-omega-omegal)*(1.+z)**2+omegal)**(-1./6.)
c	ipower=0
c
c     	sigmadla=sqrt(anorm*(
c     *	rombint(sigmagaus,dlog(0.001/scale),dlog(0.1d0/scale),tol)+
c     *	rombint(sigmagaus,dlog(0.1d0/scale),dlog(1d0/scale),tol)+
c     *	rombint(sigmagaus,dlog(1d0/scale),dlog(10.d0/scale),tol)+
c     *	rombint(sigmagaus,dlog(10.d0/scale),dlog(20.d0/scale),tol)))
c
c	write(6,54) sigmadla,scale*h,z
c54	FORMAT('sigma_DLA: ',E12.5,'    at',E12.5,
c     *		' h^-1 Mpc,      z=',F5.2)
c
c	Upper limit on Omega_gas via Press-Schecter
c	Press-Schechter threshold
c
c	deltac = 1.33
c	omegagas = omega*f_baryon*erfcc(deltac/(sqrt(2.)*sigmadla))
c
c	write(6,55) omegagas,deltac,z
c55	FORMAT('Omega_gas:<',E12.5,'    at ',F5.3,
c     *		' gaus. deltac,        z=',F5.2)
c
c	Sigma_Lya:  Mass fluctuations in a gaussian of 
c	            1/(24.04 sqrt(0mega h**2)) Mpc at 
c	            redshift z=3; CDM+baryon mass 
c
c	scale = 1./(24.04*sqrt(omhh))
c
c	if (iz.eq.0) z=3
c	ipower=1
c
c        sigma_Lya=
c     *	    sqrt(anorm*(
c     *	    rombint(sigmagaus,dlog(0.001/scale),dlog(0.1d0/scale),tol)+
c     *	    rombint(sigmagaus,dlog(0.1d0/scale),dlog(1d0/scale),tol)+
c     *	    rombint(sigmagaus,dlog(1d0/scale),dlog(10.d0/scale),tol)+
c     *	    rombint(sigmagaus,dlog(10.d0/scale),dlog(20.d0/scale),tol)))
c	write(6,53) sigma_Lya,scale*h,z
c53	FORMAT('sigma_Lya: ',E12.5,'    at',E12.5,
c     *		' h^-1 Mpc,      z=',F5.2)
c
c
c	Sigmav_40:  Rms velocity in top hat of 100h^{-1} Mpc, z=0
c		     default; CDM+baryon
c
c	scale = 40./h
c	if (iz.eq.0) z=0
c	ipower=1
c
c        sigmav_40=
c     *	    sqrt(anorm*(
c     *	    rombint(sigmavtop,dlog(0.001/scale),dlog(0.1d0/scale),tol)+
c     *	    rombint(sigmavtop,dlog(0.1d0/scale),dlog(1d0/scale),tol)+
c     *	    rombint(sigmavtop,dlog(1d0/scale),dlog(10.d0/scale),tol)+
c     *	    rombint(sigmavtop,dlog(10.d0/scale),dlog(100.d0/scale),tol)))
c	write(6,60) sigmav_40,scale*h,z
c60	FORMAT('sigma_v:    ',F6.2,' km/s    at ',F6.2,
c     *		' h^-1 Mpc,           z=',F5.2)c
c
c	write(6,*) '     ----------------------------------------------'
c	write(6,*) 'Extras:'
c	write(6,*) '  '
c
c	if (iz.eq.0) z = 0.
c
c	  yp = primordial Helium mass fraction
c	  yp = 0.24
c	  xe = ionization fraction 1                 = fully ionized
c			           (1-yp)/(1-yp/2)   = HeI
c			           (1-3yp/4)/(1-yp/2)= HeII
c	  xe = (1.-yp)/(1.-yp/2)
c	  neutral helium
c	  tau= 4.61E-2*(1.-yp/2.)*xe*omegab*h/omega**2
c
c	if (omegal.eq.0) then
c	  omegae= omega-1.d-6
c	  arg = 2.d0*(1.d0-omegae)/omegae/(1.+z) + 1.d0
c     	  age = 9.7778/h*( 
c     *	        sqrt(1.+omegae*z)/(1.-omegae)/(1.+z) -
c     *          omegae/2./(1.-omegae)**(1.5)*dlog(
c     *			arg + dsqrt(arg**2 - 1.d0))) 
c	  write(6,56) z,age
c
c	  tau= tau*(2.-3.*omega+sqrt(1.+omega*z)*(omega*z+3.*omega-2.))
c
c	  if (z.gt.0) write(6,57) z,tau
c       
c	else
c	
c	if (abs(1.-omega-omegal).lt.1E-10) then
c	  g = sqrt(omega*(1.+z)**3+1.-omega)
c	  omegaz = omega*(1.+z)**3/g**2
c	  
c	  age = 2./3.*9.7778/h/g/sqrt(1.-omegaz) 
c     *  	*log((1.+sqrt(1.-omegaz))/sqrt(omegaz))
c	  write(6,56) z,age
c
c	  tau= tau*omega*(sqrt(1.-omega+omega*(1.+z)**3) -1.)
c
c	  if (z.gt.0) write(6,57) z,tau
c	  
c	end if
c
c	end if
c
c56 	FORMAT('Age(z=',F5.2,'):    ',F6.2,' Gyrs')
c57 	FORMAT('tau(z=',F5.2,'):    ',F9.5,' (HeI)')
c
cc	if (z.gt.0) write(6,58) z,1./DD0
c
c58	FORMAT('Growth(z=',F5.2,'-0):',F5.2)
c
c
c	write(6,59) sound_horizon*h
c59	FORMAT('sound_horizon:  ',F7.2,' h^-1 Mpc')
c
c
c	call Growth(z,1.d5,dum1,dum2,dum3)
c	write(6,61) alpha_nu*dum1,z
c61	FORMAT('suppression:    ',F7.5,' z=',F5.2)
	end


       subroutine TFset_parameters
	implicit real*8 (a-h,k,o-z)
	
	real*8 N_nu
        common/GLOBALVARIABLES/omhh,f_nu,f_baryon,N_nu,y_d,alpha_nu,
     *     beta_c,sound_horizon,theta_cmb,omega,omegal,z_equality

c Auxiliary variable

        obhh = omhh*f_baryon
c Main variables

        z_equality = 2.50e4*omhh*theta_cmb**(-4.) - 1.D0
        k_equality = 0.0746*omhh*theta_cmb**(-2.)
          z_drag = 0.313*omhh**(-0.419)*(1.+0.607*omhh**(0.674))
          z_drag = 1e0 + z_drag*obhh**(0.238*omhh**(0.223))
        z_drag = 1291e0 * omhh**(0.251)/
     &           (1e0 + 0.659*omhh**(0.828)) * z_drag

	y_d = (1.+z_equality)/(1.+z_drag)

        R_drag = 31.5*obhh*theta_cmb**(-4.)*1000e0/(1e0 + z_drag)
        R_equality = 31.5*obhh*theta_cmb**(-4.)
     &               *1000e0/(1e0 + z_equality)

        sound_horizon = 2./3./k_equality*sqrt(6./R_equality)*
     &      log(( sqrt(1.+R_drag)+sqrt(R_drag+R_equality) )
     &       /(1.+sqrt(R_equality)))

	 p_c  = -(5.-sqrt(1.+24*(1.-f_nu-f_baryon)))/4.
	 p_cb = -(5.-sqrt(1.+24*(1.-f_nu)))/4.
	 f_c  = 1.-f_nu-f_baryon
	 f_cb = 1.-f_nu
	 f_nub= f_nu+f_baryon
         if(f_cb.eq.0.d0) f_cb=f_cb+1.d-10
         if(4.*p_cb+5.0.eq.0.0) p_cb=p_cb+1.d-10
         if(1.-.193*sqrt(f_nu)+.169*f_nu.eq.0.0)f_nu=f_nu+1.d-10
         if(4.*p_c+3.0.eq.0.0)p_c=p_c+1.d-10
         if(4.*p_cb+7.0.eq.0.0)p_cb=p_cb+1.d-10
	 alpha_nu= (f_c/f_cb)* (2.*(p_c+p_cb)+5.)/(4.*p_cb+5.)
	 alpha_nu= alpha_nu*(1.-0.553*f_nub+0.126*f_nub**3)
	 alpha_nu= alpha_nu/(1.-0.193*sqrt(f_nu)+0.169*f_nu)
	 alpha_nu= alpha_nu*(1.+y_d)**(p_c-p_cb)
	 alpha_nu= alpha_nu*(1.+ (p_cb-p_c)/2. *
     *			(1.+1./(4.*p_c+3.)/(4.*p_cb+7.))/(1.+y_d))
	 beta_c=1.d0/(1.d0-0.949d0*f_nub)

	 return

	end 


	real*8 function TF_master(k)

	implicit real*8 (a-h,k,o-z)
	real*8 N_nu
        common/GLOBALVARIABLES/omhh,f_nu,f_baryon,N_nu,y_d,alpha_nu,
     *     beta_c,sound_horizon,theta_cmb,omega,omegal,z_equality

	 q = k*theta_cmb**2/omhh
         if(alpha_nu.lt.0.d0) alpha_nu=0.d0       ! fix put in 13/03/02 OE
	 gamma_eff=(sqrt(alpha_nu) + (1.-sqrt(alpha_nu))/
     *		(1.+(0.43*k*sound_horizon)**4))
	 q_eff = q/gamma_eff
	 TF_master= dlog(dexp(1.d0)+1.84*beta_c*sqrt(alpha_nu)*q_eff)
	 TF_master = TF_master/(TF_master + q_eff**2*
     *	             (14.4 + 325./(1.+60.5*q_eff**1.11)))
c     AK modified
	 q_nu = 4.9*q*sqrt(N_nu/f_nu)
c     AK modified
	 TF_master = TF_master*
     *     (1.+(1.2*f_nu**(0.68)*N_nu**(0.3+0.6*f_nu))/
     *	   (q_nu**(-1.3)+q_nu**(0.45)))

	 return 

	end

	subroutine Growth(z,k,DD_cb,DD_cbnu,DD0)

	implicit real*8 (a-h,k,o-z)
	real*8 N_nu
        
        common/GLOBALVARIABLES/omhh,f_nu,f_baryon,N_nu,y_d,alpha_nu,
     *     beta_c,sound_horizon,theta_cmb,omega,omegal,z_equality       

        q = k*theta_cmb**2/omhh

        y_fs = 17.2*f_nu*(1.+0.488*f_nu**(-7./6.))*(N_nu*q/f_nu)**2

	 oz=omega*(1.+z)**3
     *     /(omegal+(1.-omegal-omega)*(1.+z)**2+omega*(1.+z)**3)
         olz=omegal
     *     /(omegal+(1.-omegal-omega)*(1.+z)**2+omega*(1.+z)**3)

	 D = (1.+z_equality)/(1.+z)*5.*oz/2.*(oz**(4./7.)-olz+(1.+oz/2.)*
     *		(1.+olz/70.))**(-1.)

	 DD0= D/((1.+z_equality)*5.*omega/2.*(omega**(4./7.)-omegal+
     *          (1.+omega/2.)*(1.+omegal/70.))**(-1.))

	 p_cb = -(5.-sqrt(1.+24*(1.-f_nu)))/4.

	 DD_cb = (1.+(D/( 1.+y_fs))**(0.7))**(-p_cb/0.7)*D**(p_cb)

	 DD_cbnu = ((1.-f_nu)**(-0.7/p_cb)
     *		+(D/(1.+y_fs))**(0.7))**(-p_cb/0.7)*D**(p_cb)

	return
	end 

        real*8 function sound_horizon_fit(omhh,f_baryon)
	implicit real*8 (a-h,k,o-z)
	
         obhh = f_baryon*omhh
         sound_horizon_fit = 44.5*log(9.83/omhh)
     &                      /sqrt(1.+10.0*obhh**(0.75))

        return
        end


	real*8 function sigmatop(kl)

	implicit real*8 (a-h,k,o-z)
	common/PASSSIGMA/scale,tilt,z,ipower

	 k = exp(kl)
	 x = scale*k

         call Growth(z,k,DD_cb,DD_cbnu,DD0)

	if (ipower.eq.0) then
         sigmatop = k**(3.+tilt)*(DD0*DD_cbnu*TF_master(k))**2
     *                *(3.*(x*dcos(x) - dsin(x))/x**3)**2
	else

         sigmatop = k**(3.+tilt)*(DD0*DD_cb*TF_master(k))**2
     *                *(3.*(x*dcos(x) - dsin(x))/x**3)**2
	end if

	return
	end

	real*8 function sigmagaus(kl)

	implicit real*8 (a-h,k,o-z)
	common/PASSSIGMA/scale,tilt,z,ipower

	 k = exp(kl)
	 x = scale*k

         call Growth(z,k,DD_cb,DD_cbnu,DD0)

	if (ipower.eq.0) then
         sigmagaus = k**(3.+tilt)*(DD0*DD_cbnu*TF_master(k))**2
     *                *exp(-x**2.) 
	else
         sigmagaus = k**(3.+tilt)*(DD0*DD_cb*TF_master(k))**2
     *                *exp(-x**2.) 
	end if

	return
	end

	real*8 function sigmavtop(kl)

	implicit real*8 (a-h,k,o-z)
	real*8 N_nu
        common/GLOBALVARIABLES/omhh,f_nu,f_baryon,N_nu,y_d,alpha_nu,
     *     beta_c,sound_horizon,theta_cmb,omega,omegal,z_equality
	common/PASSSIGMA/scale,tilt,z,ipower

	 k = exp(kl)
	 x = scale*k
	 q = k*theta_cmb**2/omhh

         call Growth(z,k,DD_cb,DD_cbnu,DD0)
	 h = sqrt(omhh/omega)
	 g = sqrt(omega*(1.+z)**3+(1.-omega-omegal)*(1.+z)**2+
     *		omegal)
	 omegaz = omega*(1.+z)**3/g**2
	 omegalz = omegal/g**2
	 p_cb = -(5.-sqrt(1.+24*(1.-f_nu)))/4.
	 y_fs = 17.2*f_nu*(1.+0.488*f_nu**(-7./6.))*(N_nu*q/f_nu)**2
	 D = (1.+z_equality)/(1.+z)*5.*omegaz/2.*(omegaz**(4./7.)
     *	        -omegalz+(1.+omegaz/2.)*
     *		(1.+omegalz/70.))**(-1.)

	 f = omegaz**(0.6) + 1./70.*omegalz *(1.+omegaz/2.)

	if (ipower.eq.0) then

	 f = f*(1.+p_cb/((1.-f_nu)**(-0.7/p_cb)+(D/(1.+y_fs))**0.7)) 
         sigmavtop = k**(3.+tilt)*(DD0*DD_cbnu*TF_master(k))**2
     *                *(3.*(x*dcos(x) - dsin(x))/x**3)**2
	else
	 f = f*(1.+p_cb/(1.+(D/(1.+y_fs))**0.7)) 
         sigmavtop = k**(3.+tilt)*(DD0*DD_cb*TF_master(k))**2
     *                *(3.*(x*dcos(x) - dsin(x))/x**3)**2
	end if

	 sigmavtop = (f*100.*h*g/(1.+z)/k)**2*sigmavtop

	return
	end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        real*8 function rombint(f,a,b,tol)
c  Rombint returns the integral from a to b of using Romberg integration.
c  The method converges provided that f(x) is continuous in (a,b).
c  f must be double precision and must be declared external in the calling
c  routine.  tol indicates the desired relative accuracy in the integral.
c
c Taken from CMBFAST; Seljak & Zaldarriaga (1996)


        parameter (MAXITER=30,MAXJ=5)
        implicit real*8 (a-h,o-z)
        dimension g(MAXJ+1)
        real*8 f
        external f
c
        h=0.5d0*(b-a)
        gmax=h*(f(a)+f(b))
        g(1)=gmax
        nint=1
        error=1.0d20
        i=0
10        i=i+1
          if (i.gt.MAXITER.or.(i.gt.5.and.abs(error).lt.tol))
     2      go to 40
c  Calculate next trapezoidal rule approximation to integral.
          g0=0.0d0
            do 20 k=1,nint
            g0=g0+f(a+(k+k-1)*h)
20        continue
          g0=0.5d0*g(1)+h*g0
          h=0.5d0*h
          nint=nint+nint
          jmax=min(i,MAXJ)
          fourj=1.0d0
            do 30 j=1,jmax
c  Use Richardson extrapolation.
            fourj=4.0d0*fourj
            g1=g0+(g0-g(j))/(fourj-1.0d0)
            g(j)=g0
            g0=g1
30        continue
          if (abs(g0).gt.tol) then
            error=1.0d0-gmax/g0
          else
            error=gmax
          end if
          gmax=g0
          g(jmax+1)=g0
        go to 10
40      rombint=g0
        if (i.gt.MAXITER.and.abs(error).gt.tol)
     2    write(*,*) 'Rombint failed to converge; integral, error=',
     3    rombint,error
        return
        end


      real*8 FUNCTION erfcc(x)
c calculates the complimentary error function
c Taken from Numerical Recipes.  
      implicit real*8 (a-h,o-z)

      REAL*8 x
      REAL*8 t,z
      z=abs(x)
      t=1./(1.+0.5*z)
      erfcc=t*exp(-z*z-1.26551223+t*(1.00002368+t*(.37409196+t*
     *(.09678418+t*(-.18628806+t*(.27886807+t*(-1.13520398+t*
     *(1.48851587+t*(-.82215223+t*.17087277)))))))))
      if (x.lt.0.) erfcc=2.-erfcc
      return
      END


