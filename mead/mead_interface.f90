module mead_settings_mod
	type mead_settings
		logical :: noisy
		real(8) :: kmin, kmax
		integer :: nk

		!real(8) :: numin, numax

		real(8) :: zmin, zmax
		integer :: nz

		logical :: feedback

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

	status = status + datablock_get_logical_default(options, option_section, "feedback", .false., settings%feedback)

	if (status .ne. 0) then
		write(*,*) "One or more parameters not found for hmcode"
		stop
	endif

	WRITE(*,*) 'z min:', settings%zmin
	WRITE(*,*) 'z max:', settings%zmax
	WRITE(*,*) 'number of z:', settings%nz
	WRITE(*,*)

	WRITE(*,*) 'k min:', settings%kmin
	WRITE(*,*) 'k max:', settings%kmax
	WRITE(*,*) 'number of k:', settings%nk
	WRITE(*,*)


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
    USE camb_stuff

	implicit none

	integer(cosmosis_block), value :: block
	integer(cosmosis_status) :: status
	type(c_ptr), value :: config
	type(mead_settings), pointer :: settings	
	integer, parameter :: LINEAR_SPACING = 0
	integer, parameter :: LOG_SPACING = 1
	character(*), parameter :: cosmo = cosmological_parameters_section
	character(*), parameter :: halo = halo_model_parameters_section
	character(*), parameter :: linear_power = matter_power_lin_section
	character(*), parameter :: nl_power = matter_power_nl_section

	!real(4) :: p1h, p2h,pfull, plin, z
	integer :: i,j, z_index, nk_lin
	REAL, ALLOCATABLE :: k(:),  pk(:,:), ztab(:), atab(:)
	REAL, ALLOCATABLE :: k_lin(:),  pk_lin(:), z_lin(:), a_lin(:)
	TYPE(cosmology) :: cosm
	!TYPE(tables) :: lut
	!CosmoSIS supplies double precision - need to convert
	real(8) :: Om_m, Om_lam, Om_b, h, w, sig8, n_s, Om_nu
	real(8), ALLOCATABLE :: k_in(:), z_in(:), p_in(:,:), a_in(:)
	real(8), ALLOCATABLE :: k_out(:), z_out(:), a_out(:), p_out(:,:)
	!real(8) :: Halo_as, halo_eta0
	real(8) :: log10T_AGN
	INTEGER, PARAMETER :: icos_default = 1
	INTEGER :: icos
	INTEGER, PARAMETER :: version = HMcode2020_feedback
	LOGICAL, PARAMETER :: verbose = .TRUE.
	REAL :: kmin, kmax, zmin, zmax

	status = 0
	call c_f_pointer(config, settings)

	!feedback = settings%feedback
	!TODO if statement: if feedback -> set hmcode 2020 feedback

	! Different HMcode versions
    !version = HMcode2020
	!version = HMcode2020_feedback

	!Fill in the cosmology parameters. We need to convert from CosmoSIS 8-byte reals
	!to HMcode 4-byte reals, hence the extra bit
	status = status + datablock_get(block, cosmo, "omega_m", Om_m)
	status = status + datablock_get(block, cosmo, "omega_lambda", Om_lam)
	status = status + datablock_get(block, cosmo, "omega_b", Om_b)
        status = status + datablock_get_double_default(block, cosmo, "omega_nu", 0.0D0, Om_nu)
	status = status + datablock_get(block, cosmo, "h0", h)
	status = status + datablock_get(block, cosmo, "sigma_8", sig8)
	status = status + datablock_get(block, cosmo, "n_s", n_s)
	status = status + datablock_get_double_default(block, cosmo, "w", -1.0D0, w)


	!status = status + datablock_get_double_default(block, halo, "A", 3.13D0, halo_as)
	!status = status + datablock_get_double_default(block, halo, "eta_0", 0.603D0, halo_eta0)
	!status = status + datablock_get_double_default(block, halo, "log10T_AGN", 7.8, log10T_AGN)


	if (status .ne. 0 ) then
		write(*,*) "Error reading parameters for Mead code"
		return
	endif

	icos = icos_default
	CALL assign_cosmology(icos, cosm, verbose)

    !cosi%om_m=om_m-om_nu !The halo modelling should include only cold matter components (i.e. DM and baryons)
    !I think om_m now includes the neutrinos, but i can't really specify neutrinos as a mass... 
	!i guess i can convert omega_nu to m_nu but annoying
	!also do I need to specify upper or lower case ones?, need to add wa as well
	cosm%iw = iw_wCDM        ! Set to wCDM dark energy
    cosm%Om_v = 0.           ! Force vacuum energy density to zero (note that DE density is non-zero)
	cosm%Om_w=Om_lam
    cosm%Om_b=Om_b
    cosm%h=h
    cosm%w=w
    cosm%sig8=sig8
    cosm%ns=n_s

    !cosi%eta_0 = halo_eta0
    !cosi%As = halo_as
	!cosm%Thead=10**log10T_AGN

    !And get the cosmo power spectrum, again as double precision
    !Also the P is 2D as we get z also
	status = status + datablock_get_double_grid(block, linear_power, &
        "k_h", k_in, "z", z_in, "p_k", p_in)
	a_in = 1/(1+z_in)
	!allocate(k_lin(size(k_in)))
	!allocate(a_lin(size(a_in)))
	!allocate(z_lin(size(z_in)))
	!allocate(pk_lin(size(k_in), size(z_in)))
	!k_lin = k_in
	!a_lin = a_in!(:-1)
	!z_lin = z_in!(:-1)
	!pk_lin = p_in!(:, :-1) !* (k_lin**3.0) * (2.*(pi**2.))
	!DO i=1,size(z_lin)
	!	pk_lin(:, i) = p_in(:, i) * (k_lin**3.0) * (2.*(pi**2.))
	!END DO
	!pk_lin = Pk_Delta(pk_lin, k_lin)

	allocate(k_lin(size(k_in)))
	allocate(a_lin(1))
	allocate(z_lin(1))
	allocate(pk_lin(size(k_in)))
	k_lin = k_in
	a_lin = a_in(1)
	z_lin = z_in(1)
	pk_lin = p_in(:, 1) !* (k_lin**3.0) * (2.*(pi**2.))    
	!Pk_lin = Pk_Delta(Pk_lin, k_lin)
    CALL init_external_linear_power_tables(cosm, k_lin, a_lin, reshape(pk_lin, [nk_lin, 1]))

	if (status .ne. 0 ) then
		write(*,*) "Error reading P(k,z) for Mead code"
		return
	endif

	nk_lin = size(k_lin)
	CALL init_external_linear_power_tables(cosm, k_lin, a_lin, reshape(pk_lin, [nk_lin, 1]))
	!Copy in k
	!allocate(cosi%ktab(size(k_in)))
	!cosi%ktab = k_in

	!Find the index of z where z==0
	!if (z_in(1)==0.0) then
	!	z_index=1
	!elseif (z_in(size(z_in))==0.0) then
	!	z_index=size(z_in)
	!else
	!	write(*,*) "P(k,z=0) not found - please calculate"
	!	status = 1
	!	return
	!endif
	!Copy in P(k) from the right part of P(k,z)
	!allocate(cosi%pktab(size(k_in)))
    !cosi%pktab = p_in(:, z_index) * (cosi%ktab**3.)/(2.*(pi**2.))
    !cosi%itk = 5


	!Set the output ranges in k and z
	!CALL fill_table(real(settings%kmin),real(settings%kmax),k,settings%nk,LOG_SPACING)DELETE
	!CALL fill_table(real(settings%zmin),real(settings%zmax),ztab,settings%nz,LINEAR_SPACING)DELETE

	!create a and k arrays based on parameters
	kmin = settings%kmin
	kmax = settings%kmax
	zmin = settings%zmin
	zmax = settings%zmax
	CALL fill_array_log(kmin, kmax, k, settings%nk)
    CALL fill_array(zmin, zmax, ztab, settings%nz)
	! need to calculate a = 1/(1+z)
	atab = 1 /(1+ztab) 
    !CALL fill_array(amin, amax, a, na) DELTE

	!Fill table for output power
	!ALLOCATE(p_out(settings%nk,settings%nz))

	! Initialise cosmological model
	CALL init_cosmology(cosm)
	CALL print_cosmology(cosm)

	CALL calculate_HMcode(k, atab, pk, settings%nk, settings%nz, cosm, version=version)


	!Loop over redshifts
	!DO j=1,settings%nz

		!Sets the redshift
		!z=ztab(j)

		!Initiliasation for the halomodel calcualtion
		!Also normalises power spectrum (via sigma_8)
		!and fills sigma(R) tables
		
		!CALL halomod_init(z,real(settings%numin),real(settings%numax),lut,cosi)

		!Loop over k values
		!DO i=1,SIZE(k)
		!	plin=p_lin(k(i),cosi)        
		!	CALL halomod(k(i),z,p1h,p2h,pfull,plin,lut,cosi)
			!This outputs k^3 P(k).  We convert back.
		!	p_out(i,j)=pfull / (k(i)**3.0) * (2.*(pi**2.))
		!END DO

		!IF(j==1) THEN
		!	if (settings%feedback) WRITE(*,fmt='(A5,A7)') 'i', 'z'
		!	if (settings%feedback) WRITE(*,fmt='(A13)') '   ============'
		!END IF
		! if (settings%feedback) WRITE(*,fmt='(I5,F8.3)') j, ztab(j)
	!END DO

	!convert to double precision
	allocate(k_out(settings%nk))
	allocate(z_out(settings%nz))
	k_out = k
	z_out = ztab
	allocate(p_out(settings%nk,settings%nz))
	DO i=1,size(z_out)
		!	pk_lin(:, i) = p_in(:, i) * (k_lin**3.0) * (2.*(pi**2.))
		p_out(:, i) = Pk_Delta(pk(:,i), k)
	END DO
	!p_out = pk

	!Convert k to k/h to match other modules
	!Output results to cosmosis
	status = datablock_put_double_grid(block,nl_power, "k_h", k_out, "z", z_out, "p_k", p_out)

	!Free memory
	deallocate(k)
	deallocate(ztab)
	deallocate(p_out)
	deallocate(k_in)
	deallocate(z_in)
	deallocate(p_in)
	deallocate(k_out)
	deallocate(z_out)
	!call deallocate_LUT(lut)
    !IF(ALLOCATED(cosm%rtab)) DEALLOCATE(cosm%rtab)
    !IF(ALLOCATED(cosm%sigtab)) DEALLOCATE(cosm%sigtab)   
    !IF(ALLOCATED(cosm%ktab)) DEALLOCATE(cosm%ktab)
    !IF(ALLOCATED(cosm%tktab)) DEALLOCATE(cosm%tktab)
    !IF(ALLOCATED(cosm%pktab)) DEALLOCATE(cosm%pktab)

end function execute
