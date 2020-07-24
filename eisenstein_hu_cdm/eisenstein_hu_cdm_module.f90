!gfortran -O3 -g -fPIC -I/cosmosis/cosmosis/datablock  -std=gnu -ffree-line-length-none -shared 
!-o eisenstein_hu_cdm_module.so spline.f90 interface_tools.f90 eisenstein_hu_cdm.f90 
!eisenstein_hu_cdm_module.f90 /cosmosis/cosmosis/datablock/libcosmosis_fortran.so
function setup(options) result(result)
  USE cosmosis_modules
  USE interface_tools
  implicit none
  integer(cosmosis_block), value :: options
  integer(cosmosis_status) :: status
  type(ini_settings), pointer :: settings
  type(c_ptr) :: result
  allocate(settings)

  status = 0
  status = status + datablock_get_double_default(options, option_section, "zmin", 0.0D-04, settings%zmin)
  status = status + datablock_get_double_default(options, option_section, "zmax", 3.0D+00, settings%zmax)
  status = status + datablock_get_int_default(options, option_section, "nz_steps", 800, settings%nz_steps)
  status = status + datablock_get_double_default(options, option_section, "kmin", 1.0D-05, settings%kmin)
  status = status + datablock_get_double_default(options, option_section, "kmax", 10D+00, settings%kmax)
  status = status + datablock_get_int_default(options, option_section, "nk_steps", 800, settings%nk_steps)

  settings%dz=(settings%zmax-settings%zmin)/(settings%nz_steps-1.0)
  settings%dk=(log(settings%kmax)-log(settings%kmin))/(settings%nk_steps-1.0)

  print*,""
  print*,"EH power spectrum options:"
  print*, "zmin",settings%zmin
  print*, "zmax",settings%zmax
  print*, "nz_steps",settings%nz_steps
  print*, "kmin",settings%kmin
  print*, "kmax",settings%kmax
  print*, "nk_steps",settings%nk_steps
  print*,""

  if (status .ne. 0) then
     write(*,*) "Failed setup of eisenstein hu cdm", status
     stop
  endif
  result = c_loc(settings)
end function setup



! ==========================================

function execute(block, config) result(status)
  use ehu_cdm
  use interface_tools
  use cosmosis_modules
  implicit none
  integer(cosmosis_block), value :: block
  integer(cosmosis_status) :: status
  type(c_ptr), value :: config
  type(ini_settings), pointer :: settings
  type(pk_settings) :: PK
  real(8) :: omega_baryon,omega_matter,w_de,omega_de,h0,n_s,n_run,A_s, omega_nu, sigma8, anorm
  real(8) :: omega_matter_and_nu, z_for_Pk_calc, dd_cb,dd_cbnu,dd0, D_0_new
  real(8) :: k, z,D_0,yval, ypval, yppval
  real(8), allocatable, dimension(:,:) :: P
  real(8), allocatable, dimension(:) :: dz,zbins,dz_interpolated

  real(dl), dimension(:), allocatable :: Pk_at_one_redshift, k_h_linear


  integer iz,n,i, j,n_growth

  status = 0
  call c_f_pointer(config, settings)


  !  LOAD COMSMOLOGICAL PARAMETERS

  status = status + datablock_get_double(block, cosmological_parameters_section, "OMEGA_B", omega_baryon)
  status = status + datablock_get_double(block, cosmological_parameters_section, "OMEGA_M", omega_matter)
  status = status + datablock_get_double_default(block, cosmological_parameters_section, "W",-1.0d0, w_de)
  status = status + datablock_get_double(block, cosmological_parameters_section, "h0", h0)
  status = status + datablock_get_double(block, cosmological_parameters_section, "n_s", n_s)
  status = status + datablock_get_double_default(block, cosmological_parameters_section, "n_run", 0.0d0, n_run)
  status = status + datablock_get_double(block, cosmological_parameters_section, "A_s", A_s)
  n_growth= datablock_get_array_length(block, GROWTH_PARAMETERS_SECTION, "d_z")
  status = status+ datablock_get_double_array_1d(block, GROWTH_PARAMETERS_SECTION, "d_z", dz,n_growth);
  status = status+  datablock_get_double_array_1d(block, GROWTH_PARAMETERS_SECTION, "z", zbins,n_growth);

  status = status + datablock_get_double(block, cosmological_parameters_section, "OMEGA_NU", omega_nu)


  

  if (status .ne. 0) then
    write(*,*) "Error in crl_eistenstein_hu"
    return
  endif

  !    INTERPOLATE GROWTH
  allocate(dz_interpolated(n_growth))

  call spline_cubic_set ( n_growth, zbins , dz, 2, 0.0, 2, 0.0, dz_interpolated )

  if ( (abs(zbins(1) - settings%zmin) .gt. 0.1d0) .or. (settings%zmax - zbins(n_growth)).gt. 0.1d0 ) then

     print*, "======================="
     print*,  " the chosen bounds of the growth module does not cover the entire redshift range requested for the power spectrum"
     print*, "zmin zmax from growth",zbins(1),zbins(n_growth)
     print*, "zmin zmax in Pk_nowiggle",settings%zmin,settings%zmax

     print*, "======================="
  end if



  ! create an array for P(k,z)

  ! first define the dimensions
  PK%num_k = settings%nk_steps
  PK%num_z = settings%nz_steps

  ! create empty containers, defined in interface.tools.f90

  call  allocate_matterpower(PK)
  allocate(k_h_linear(PK%num_k))

  z=settings%zmin
  k=log(settings%kmin)

  !  fill k and z arrays

  do i = 1, PK%num_z, 1
     PK%redshifts(i)= z
     z =z+settings%dz
  end do

  do i = 1, PK%num_k, 1

     PK%kh(i)= exp(k)
     k =k+settings%dk

  end do

  ! finally fill with P(k,z) with EH formula
  ! cycle over k

  omega_matter_and_nu = omega_matter+omega_nu
  z_for_Pk_calc = 0.
  allocate(Pk_at_one_redshift(PK%num_k))
  !    POWER(input_omega          , input_omegab, input_omeganu, input_h, nsval, sigma8, anorm,Nk      ,Rk,Pk                 ,z)
  CALL POWER(omega_matter_and_nu, omega_baryon, omega_nu     , h0, n_s ,&
          &  sigma8, anorm,PK%num_k,PK%kh,Pk_at_one_redshift,z_for_Pk_calc)
  !notes:
  !cosmosis: omega_m = 1-omega_lambda-omega_k-omega_nu
  ! so input_omega = omega_matter + omega_nu (we need to include neutrinos)
   print *, 'Test of power at k= ', PK%kh(5)
   print *, 'gives Pk = ', Pk_at_one_redshift(5)
  do i = 1, PK%num_k, 1
     !k=PK%kh(i)
     !call   compute_pknowiggle(k,A_s,n_s,h0,omega_matter,omega_baryon,PK%matpower(i,1))
     PK%matpower(i,1) = Pk_at_one_redshift(i)
  end do
  ! cycle over z just applying the growth

  !original code from ehu module: uses cosmosis growth function
   D_0=dz_growth(PK%redshifts(1),zbins,dz,dz_interpolated)
   print *, 'original D_0: ', D_0

   PK%matpower(:,1)=PK%matpower(:,1)*(D_0)**2
!   do i = 2, PK%num_z, 1
!      z=PK%redshifts(i)
!      PK%matpower(:,i)=PK%matpower(:,1)*(dz_growth(z,zbins,dz,dz_interpolated)/D_0)**2
!   end do

  
  do j = 1, Pk%num_k, 1
   z = PK%redshifts(1)
   CALL GROWTH(z,PK%kh(j)*h0,dd_cb,dd_cbnu,dd0)
   D_0_new = dd_cbnu*dd0
   print *, 'new D_0: ', D_0_new
   !print *, 'new dd_cbnu: ', dd_cbnu
   !print *, 'new dd0: ', dd0
   !PK%matpower(j,1)=PK%matpower(j,1)*(D_0_new)**2
   do i = 2, PK%num_z, 1
      z=PK%redshifts(i)
      CALL GROWTH(z,PK%kh(j)*h0,dd_cb,dd_cbnu,dd0)
      PK%matpower(j,i)=PK%matpower(j,1)*(dd_cbnu*dd0/D_0_new)**2
   end do
  end do


  ! save in datablock
  ! change by Lukas
  status = datablock_put_double_grid(block, "matter_power_lin", &
       "k_h", PK%kh, "z", PK%redshifts, "P_k", PK%matpower)


  !! deallocate everything here.

  call deallocate_matterpower(PK)

  if (allocated(dz_interpolated)) deallocate(dz_interpolated)

end function execute






function cleanup(config) result(status)
  use interface_tools
  use cosmosis_modules
  type(c_ptr), value :: config
  type(ini_settings), pointer :: settings
  integer(cosmosis_status) :: status

  !Free memory allocated in the setup function
  call c_f_pointer(config, settings)
  deallocate(settings)

  status = 0

end function cleanup

