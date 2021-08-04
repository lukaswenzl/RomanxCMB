module mead_settings_mod
	type mead_settings
		logical :: noisy
		real(8) :: kmin, kmax
		integer :: nk

		!real(8) :: numin, numax

		real(8) :: zmin, zmax
		integer :: nz

		logical :: feedback

		!boolean to switch between input p(k,z) and p(k, z=0)
		!only meant for testing
		logical :: power_input_zdep
		!boolean switch to minimize the nl samples to calculate
		logical :: optimize_nl_samples
		!boolean switch to us cosmosis growth instead of internal mead growth
		logical :: use_cosmosis_growth
		logical :: verbose


	end type mead_settings

end module mead_settings_mod

function setup(options) result(result)
	use mead_settings_mod
	use cosmosis_modules
	implicit none

	integer(cosmosis_block), value :: options
	integer(cosmosis_status) :: status
	type(mead_settings), pointer :: settings
	type(c_ptr) :: result
	status = 0
	
	allocate(settings)

	status = status + datablock_get(options, option_section, "zmin", settings%zmin)
	status = status + datablock_get(options, option_section, "zmax", settings%zmax)
	status = status + datablock_get(options, option_section, "nz", settings%nz)


	status = status + datablock_get(options, option_section, "kmin", settings%kmin)
	status = status + datablock_get(options, option_section, "kmax", settings%kmax)
	status = status + datablock_get(options, option_section, "nk", settings%nk)

	!status = status + datablock_get_double_default(options, option_section, "numin", 0.1D0, settings%numin)
	!status = status + datablock_get_double_default(options, option_section, "numax", 5.0D0, settings%numax)

	status = status + datablock_get_logical_default(options, option_section, "feedback", .true., settings%feedback)
	status = status + datablock_get_logical_default(options, option_section, "power_input_zdep", .true., settings%power_input_zdep)
	status = status + datablock_get_logical_default(options, option_section, "optimize_nl_samples", .false., settings%optimize_nl_samples)
	status = status + datablock_get_logical_default(options, option_section, "use_cosmosis_growth", .true., settings%use_cosmosis_growth)
	status = status + datablock_get_logical_default(options, option_section, "verbose", .false., settings%verbose)



	if (status .ne. 0) then
		write(*,*) "One or more parameters not found for hmcode"
		stop
	endif

	! WRITE(*,*) 'z min:', settings%zmin
	! WRITE(*,*) 'z max:', settings%zmax
	! WRITE(*,*) 'number of z:', settings%nz
	! WRITE(*,*)

	! WRITE(*,*) 'k min:', settings%kmin
	! WRITE(*,*) 'k max:', settings%kmax
	! WRITE(*,*) 'number of k:', settings%nk
	! WRITE(*,*)


	result = c_loc(settings)

end function setup


function execute(block,config) result(status)
	use mead_settings_mod
	use cosmosis_modules
	
	USE constants
    USE basic_operations
    USE array_operations
    USE cosmology_functions
    USE HMx
	!USE interpolate
    USE camb_stuff

	implicit none

	integer(cosmosis_block), value :: block
	integer(cosmosis_status) :: status
	type(c_ptr), value :: config
	type(mead_settings), pointer :: settings	
	integer, parameter :: LINEAR_SPACING = 0
	integer, parameter :: LOG_SPACING = 1
	character(*), parameter :: cosmo = cosmological_parameters_section
	character(*), parameter :: linear_power = matter_power_lin_section
	character(*), parameter :: nl_power = matter_power_nl_section
	character(*), parameter :: growth_section = "growth_parameters"

	integer :: i,j, z_index, nk_lin
	REAL, ALLOCATABLE :: k(:),  pk(:,:), ztab(:), atab(:)
	REAL, ALLOCATABLE :: k_lin(:),  pk_lin(:,:),pk_lin_k_only(:), z_lin(:), a_lin(:)
	REAL, ALLOCATABLE :: growth(:), rate(:), a_growth(:)
	TYPE(cosmology) :: cosm
	!CosmoSIS supplies double precision - need to convert
	real(8) :: Om_m, Om_lam, Om_b, h, w,wa, sig8, n_s, Om_nu, omnuh2
	real(8), ALLOCATABLE :: k_in(:), z_in(:), p_in(:,:), a_in(:)
	real(8), ALLOCATABLE :: k_out(:), z_out(:), a_out(:), p_out(:,:)
	real(8), ALLOCATABLE :: growth_in(:), rate_in(:), a_growth_in(:)
	integer :: n_growth

	!real(8) :: Halo_as, halo_eta0
	real(8) :: log10T_AGN
	INTEGER, PARAMETER :: icos_default = 1
	INTEGER :: icos
	INTEGER :: version != HMcode2020_feedback
	!LOGICAL, PARAMETER :: verbose = .TRUE.
	REAL :: kmin, kmax, zmin, zmax
	integer :: k_idx_cutoff
	REAL, ALLOCATABLE :: k_less(:), norm_lin(:)

	status = 0
	call c_f_pointer(config, settings)


	! feedback setting switches between HMcode versions
	if (settings%feedback ) then
		version = HMcode2020_feedback
		IF (settings%verbose) THEN
			WRITE (*, *) "HMcode2020_feedback"
		 END IF
	else 
		version = HMcode2020
		IF (settings%verbose) THEN
			WRITE (*, *) "HMcode2020"
		 END IF
	endif

	!Fill in the cosmology parameters. We need to convert from CosmoSIS 8-byte reals
	!to HMcode 4-byte reals, hence the extra bit
	status = status + datablock_get(block, cosmo, "omega_m", Om_m)
	status = status + datablock_get(block, cosmo, "omega_lambda", Om_lam)
	status = status + datablock_get(block, cosmo, "omega_b", Om_b)
	status = status + datablock_get(block, cosmo, "h0", h)
	status = status + datablock_get(block, cosmo, "sigma_8", sig8)
	status = status + datablock_get(block, cosmo, "n_s", n_s)
	status = status + datablock_get_double_default(block, cosmo, "w", -1.0D0, w)
	status = status + datablock_get_double_default(block, cosmo, "wa", 0.0D0, wa)
	status = status + datablock_get_double_default(block, cosmo, "omnuh2", 0.0D0, omnuh2)

	!The log10T_AGN is only used for HMcode2020_feedback, it is the free parameter for the baryon model
	status = status + datablock_get_double_default(block, cosmo, "log10T_AGN", 7.8D0, log10T_AGN)


	if (status .ne. 0 ) then
		write(*,*) "Error reading parameters for Mead code"
		return
	endif

	icos = icos_default
	CALL assign_cosmology(icos, cosm, settings%verbose)

	!note: in Mead2016 the code defined the internal omega_m without neutrinos. This I think has changed and here we
	!include neutrinos in omega_m so we can directly map cosm%Om_m=Om_m
	!also note: lower case is multiplied by h2 (eg om_m), upper case is not (eg Om_m)
	cosm%iw = iw_waCDM        ! Set to w_waCDM dark energy
    cosm%Om_v = 0.           ! Force vacuum energy density to zero (note that DE density is non-zero)
	cosm%Om_w=Om_lam
	cosm%w=w
	cosm%wa=wa

	cosm%Om_m=Om_m
    cosm%Om_b=Om_b
    cosm%h=h
    cosm%sig8=sig8
    cosm%ns=n_s
	cosm%m_nu= omnuh2 * 93.14

	cosm%Theat=10**log10T_AGN

    !And get the cosmo power spectrum, again as double precision
    !Also the P is 2D as we get z also
	status = status + datablock_get_double_grid(block, linear_power, &
        "k_h", k_in, "z", z_in, "p_k", p_in)
	a_in = 1/(1+z_in)

	if (settings%power_input_zdep) then
		allocate(k_lin(size(k_in)))
		allocate(a_lin(size(a_in)))
		allocate(z_lin(size(z_in)))
		allocate(pk_lin(size(k_in), size(z_in)))
		k_lin = k_in
		a_lin = a_in
		z_lin = z_in
		pk_lin = p_in
		!note we are putting in the power spectrum. In future updates this function will likely change
		!to the internal Delta^2, so for future updates a factor k^3 * 2pi^2 might be needed. 
		CALL init_external_linear_power_tables(cosm, k_lin, a_lin, pk_lin)
	else
		write(*,*) "This code should not be used its only for testing. Please set power_input_zdep=True for HMcode "
		allocate(k_lin(size(k_in)))
		allocate(pk_lin_k_only(size(k_in)))
        allocate(a_lin(1))
        a_lin = 1.         
		pk_lin_k_only = p_in(:, 1)
        CALL init_external_linear_power_tables(cosm, k_lin, a_lin, reshape(pk_lin_k_only, [nk_lin, 1]))
	endif

	if (status .ne. 0 ) then
		write(*,*) "Error reading P(k,z) for Mead code"
		return
	endif

	! Initialise cosmological model
	CALL init_cosmology(cosm)
	CALL print_cosmology(cosm)

	!if desired can pipe the growth from cosmosis into the mead module. For wCDM with neutrinos 
	!there shouldn't be any difference, but when using modified growth this can have an effect
	!note this needs to be after init_cosmology
	if(settings%use_cosmosis_growth) then
		status = status + datablock_get_double_array_1d(block, growth_section, "d_z", growth_in, n_growth)
		status = status + datablock_get_double_array_1d(block, growth_section, "f_z", rate_in, n_growth)
		status = status + datablock_get_double_array_1d(block, growth_section, "a", a_growth_in, n_growth)
		IF (settings%verbose) THEN
			WRITE (*, *) 'loaded growth from cosmosis pipeline'
		ENDIF
		growth = growth_in
		rate = rate_in
		a_growth = a_growth_in
		!WRITE (*, *) 'check is has growth is true before putting external growth: ', cosm%has_growth
		!WRITE (*, *) 'check normalization before putting external growth: ', growth(1)

		CALL init_external_growth(a_growth, growth, rate, cosm)
		!WRITE (*, *) 'check is has growth is true: ', cosm%has_growth
	endif

	if (.NOT. settings%optimize_nl_samples) then
		write(*,*) "WARNING: Bug in code leads to unphysical suppression for k< 1e-3, consider using a cutoff"
		!create a and k arrays based on parameters
		kmin = settings%kmin
		kmax = settings%kmax
		zmin = settings%zmin
		zmax = settings%zmax
		CALL fill_array_log(kmin, kmax, k, settings%nk)
		CALL fill_array(zmin, zmax, ztab, settings%nz)
		! need to calculate a = 1/(1+z)
		atab = 1 /(1+ztab) 

		CALL calculate_HMcode(k, atab, pk, settings%nk, settings%nz, cosm, version=version)

		!convert to double precision
		allocate(k_out(settings%nk))
		allocate(z_out(settings%nz))
		k_out = k
		z_out = ztab
		allocate(p_out(settings%nk,settings%nz))
		DO i=1,size(z_out)
			! To convert from the internal Delta^2 convention back to P(k) we use 
			! the included function
			p_out(:, i) = Pk_Delta(pk(:,i), k)
		END DO
	else
		!case optimize_nl_samples -> we will only calculate the nl power at key points and rescale the linear model

		IF (size(k_lin) .ne. settings%nk) THEN
			write(*,*) "ERROR: linear and non-linear power spectrum need same sampling to use optimize nl samples!"
			return
		ENDIF
		!create a and k arrays based on parameters
		kmin = settings%kmin
		kmax = settings%kmax
		zmin = settings%zmin
		zmax = settings%zmax
		CALL fill_array_log(kmin, kmax, k, settings%nk)
		CALL fill_array(zmin, zmax, ztab, settings%nz)
		! need to calculate a = 1/(1+z)
		atab = 1 /(1+ztab) 

		k_idx_cutoff = 90 
		!this index will depend on your sampling
		!should chose so between k = 1e-2 and 1e-3
		k_less = k(k_idx_cutoff:)
		!write(*,*) "k values:", k_less

		CALL calculate_HMcode(k_less, atab, pk, size(k_less), settings%nz, cosm, version=version)

		norm_lin = pk_lin(k_idx_cutoff,:) / Pk_Delta(pk(1, :), k_less(1))
		!write(*,*) "normalization difference:", norm_lin

		DO j=1,size(k_less)
			pk(j,:) = pk(j,:) * norm_lin
		END DO

		allocate(k_out(settings%nk))
		allocate(z_out(settings%nz))
		k_out = k
		z_out = ztab
		allocate(p_out(settings%nk,settings%nz))



		DO i=1,size(z_out)

			!check that norm_lin is not far from 1
			IF (norm_lin(i) .gt. 1.01 .or. norm_lin(i) .lt. 0.99) THEN
				write(*,*) "WARNING: Norm between lin and nl power large: ratio =", norm_lin(i)
				write(*,*) "at redshift =", z_out(i)
			ENDIF

			! To convert from the internal Delta^2 convention back to P(k) we use 
			! the included function
			DO j=k_idx_cutoff,settings%nk
				p_out(j, i) = Pk_Delta(pk(j-k_idx_cutoff+1,i), k_less(j-k_idx_cutoff+1))
				!write(*,*) "k difference:", k_less(j-k_idx_cutoff+1) -k(j) -> check, is zero
			END DO
		
			! enforce linar solution on large scales
			DO j=1,k_idx_cutoff
				p_out(j, i) = p_in(j,i) !Pk_Delta(pk(j,i), k_less(j))
			END DO
		END DO

	endif

	status = status + cosm%status

	!Convert k to k/h to match other modules
	!Output results to cosmosis
	status = status + datablock_put_double_grid(block,nl_power, "k_h", k_out, "z", z_out, "p_k", p_out)

	if (status .ne. 0 ) then
		write(*,*) "Mead pipeline failed on parameters"
		!print out parameters
		cosm%verbose= .true.
		CALL print_cosmology(cosm)
	endif

	if (status_shared .ne. 0 ) then
		status = status + 1
		write(*,*) "Integrations inside of Mead pipeline failed ", status_shared, " times on parameters: "
		!print out parameters
		cosm%verbose= .true.
		CALL print_cosmology(cosm)
	endif



	!Free memory
	deallocate(k)
	deallocate(ztab)
	deallocate(p_out)
	deallocate(k_in)
	deallocate(z_in)
	deallocate(p_in)
	deallocate(k_out)
	deallocate(z_out)

end function execute
