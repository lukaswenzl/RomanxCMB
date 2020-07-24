!Modified by Lukas Wenzl to become a test function of the eisenstein and hu cdm module

!this is meant as a standalone test program for the power spectrum calculation.


!based on "http://www.mpa-garching.mpg.de/~komatsu/crl/" January 1, 2012: E.Komatsu but applied to different power spectrum model.
!compute_pk for CosmoSIS

! A sample program for computing the linear power spectrum
! 
! - k is in units of h Mpc^-1
! - P(k) is in units of h^-3 Mpc^3

! gfortran test_calculate_pk.f90 eisenstein_hu_cdm.f90
! ./a.out



PROGRAM test_ehu_cdm_pk
    !USE interface_tools !maybe I need this?
    USE ehu_cdm
    IMPLICIT none
  
  !private
  !public compute_pk_edu_cdm
  
  !contains
  
  
    !subroutine compute_pk_edu_cdm(k_ov_h,A_s,n_s,h0,om0,ob0,pk) !todo add neutrinos!
        ! calculate the power spectrum for one specific k
      !  k_ov_h in h Mpc^-1

      integer, parameter :: dl=8
      real(dl) :: k_ov_h,A_s,h0,om0,ob0,n_s  ! INTENT(IN) 
      !real(dl) :: pk !INTENT(INOUT)
      DOUBLE PRECISION :: trans

      !to delete
      INTEGER :: numk
      !real(8) :: kmax, input_OMEga, input_OMEgal, omeganu, omegab, h
      real(8) :: kmin, kmax, input_omega, input_omeganu, input_omegab, input_h, nsval
      real(8) :: t_cmb, input_N_Nu
      real(8) :: K,Dd_cb,Dd_cbnu,Dd0
      real(8), dimension(:), allocatable :: transfer, k_arr

      INTEGER Nk
      DOUBLE PRECISION Rk(100), Pk(100)
      DOUBLE PRECISION anorm, z, sigma8

      h0 = 0.6727d0
      k_ov_h = 1 / h0
      A_s = 2.1E-9
      om0 = 0.3156d0
      ob0 = 0.0491685d0



      !ORIGINAL CODE:
      ! CALL eisensteinhu(k_ov_h*h0,om0*h0**2d0,ob0/om0,trans)
      ! Eq.(74) of Komatsu et al., ApJS, 180, 330 (2009) with kWMAP=0.05/Mpc
      ! Remember that ak is in units of h/Mpc whereas "k" in Eq.(74) is in units
      ! of 1/Mpc. Be carefull you are using the same definition of As in other modules
      ! i.e same k_pivot


      Nk = 100
      input_omega = 0.3 
      !input_omegal = 0.7
      !N_Nu = 1. !* input_N_Nu
      input_omeganu = 0.0
      input_omegab = 0.05
      input_h = 0.5
      t_cmb = 2.726
      kmax = 10.
      kmin = 1.E-5
      numk = 50 
      input_N_Nu = 3
      z = 0.
      nsval = 0.9645
      !CALL calculate_transfer(input_OMEga, input_OMEgal, omeganu, omegab, h, t_cmb,input_N_Nu, kmin, &
      !    & kmax, numk, transfer, k_arr, sigma_8)! note omega 0 is sum with neutrinos! todo

      !print *, 'Test of transfer gives ', transfer(5)
      !print *, 'Sigma_8  ', sigma_8

      CALL POWER(input_omega, input_omegab, input_omeganu, input_h, nsval, sigma8, anorm,Nk,Rk,Pk,z)
      print *, 'Test of power at k= ', Rk(5)
      print *, 'gives Pk = ', Pk(5)
      print *, 'Sigma_8  ', sigma8
      print *, 'normalization ', anorm ! I think this is A_s but which scale is unclear... cobe likely did not use k=0.05

      z = 21.1
      K = 1.
      CALL GROWTH(z,K,Dd_cb,Dd_cbnu,Dd0)

      print *, 'Dd_cbnu at redshift 21.1:  ', Dd_cbnu
      print *, 'Dd at redshift 21.1:  ', Dd0


  
      !pk=A_s*(2d0*k_ov_h**2d0*2998d0**2d0/5d0/om0)**2d0 &
      !     *trans**2d0*(k_ov_h*h0/0.05d0)**(n_s-1d0) &
      !     *2d0*3.14159265359d0**2d0/k_ov_h**3d0
  
      !  I commented growth_z0**2d0 so you have to put the z dependence by hands.
    !end subroutine compute_pk_edu_cdm
  
  
  END PROGRAM test_ehu_cdm_pk
  