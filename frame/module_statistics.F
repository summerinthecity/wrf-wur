#ifdef WRFSTAT
MODULE module_statistics

! Version 1.2.3.2
! Takanobu Yamaguchi
! tak.yamaguchi@noaa.gov
! Since 2011/01
! Implemented in WRF3.3.1 in 2012/03
! - This module outputs the statistical output based on the mass-weighted horizontal mean with the
!   height coordinate.
! - SUBROUTINE wrf_statistics_init should be placed after the WRF initialization and before the
!   integration.
! - SUBROUTINE wrf_statistics shoud be placed where other subroutines that save the WRF history file
!   are called in the SUBROUTINE integrate in /WRF/frame/module_integrate.F.
!
! - This module is created based on the source code of SAM6.8.2.
! - The name of subroutines and variables are based on SAM6.8.2.
! - This module is intended to be used for idealized LES and CSRM simulations with the flat surface.
! - This module is memory intense due to allocating many 3D arrays in order to store local value for
!   horizontal average in the height coordinate.
! - Do not use nesting.
! - Althoguth this module should work with OMP since there is no OMP call in this 
!   module and ohter source codes related to this module, do not compile with OMP for safety reason.
!
! What is changed from v1.2.3.1
! - A bug fixed: ewle_eta(i,k,j) => ewle_eta(iwrf,k,jwrf) where SGS TKE advective tendency.
! - A bug fixed: w_eta(i,k+1,j) in loop k=1,nzm for condavg_mask

	IMPLICIT NONE
	
!!	SAVE   ! BE CAREFUL WHEN REMOVE THIS LINE
	PRIVATE   ! BE CAREFUL WHEN REMOVE THIS LINE
	
	!---------- Subdomain information obtained from WRF & SAM6.8.2/grid.f90 ----------
	INTEGER, SAVE :: i0, j0 ! start indices for WRF array minus 1
	INTEGER, SAVE :: nx     ! number of x direction grid, which does not contain ghost points for MPI
	INTEGER, SAVE :: ny     ! number of y direction grid, which does not contain ghost points for MPI
	INTEGER, SAVE :: nxp    ! nx+1
	INTEGER, SAVE :: nyp    ! ny+1
	INTEGER, SAVE :: nz     ! number of vertical levels for WRF: mass level at nz does not have value
	INTEGER, SAVE :: nzm    ! nz-1: output height level
	INTEGER, SAVE :: nzsh   ! nzsh = (nzm-1)*hrlayers+1: number of levels for zsh
	INTEGER, SAVE :: hrlayers   ! hrlayers-1 = number of layers between zs(k) and zs(k+1) for zsh
	INTEGER, SAVE :: nsubdomains   ! number of sumdomains/tiles/processors
	REAL, SAVE :: dt   ! time step (seconds)
	REAL, SAVE :: dx   ! horizonta grid spacing
	REAL, SAVE :: dy
	REAL, SAVE, ALLOCATABLE, DIMENSION(:) :: zs    ! height of the initial S level (output Z)
!!	REAL, SAVE, ALLOCATABLE, DIMENSION(:) :: zw    ! height of the initial W level
	REAL, SAVE, ALLOCATABLE, DIMENSION(:) :: press ! initial pressure at zs (constant with time)
	REAL, SAVE, ALLOCATABLE, DIMENSION(:) :: zsh   ! high resolution height level
	
	!---------- Variables available from WRF (unit in wrf_statistics) ----------
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: zw_eta    ! heihgt at W
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: ps_eta    ! pressure at S
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: rhod_eta  ! dry air density
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: rhom_eta  ! moist air density
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: u_eta     ! x-wind
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: v_eta     ! y-wind
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: w_eta     ! z-wind
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: tp_eta    ! potential temperature
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: qv_eta    ! water vapor
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: qc_eta    ! liquid water (condensate)
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: qi_eta    ! ice water (condensate)
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: qr_eta    ! rain water (precipitation)
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: qs_eta    ! snow (precipitation)
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: qg_eta    ! graupel (precipitation)
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: eta_eta   ! qv-qsat
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: ssr_eta   ! supersaturation
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: qnc_eta   ! cloud water number concentration
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: qni_eta   ! ice number concentration
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: qnr_eta   ! rain number concentration
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: qns_eta   ! snow number concentration
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: qng_eta   ! graupel number concentration
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: qnn_eta   ! CCN number concentration
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: tke_eta   ! SGS TKE
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: tkm_eta   ! SGS vertical eddy viscosity
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: tkh_eta   ! SGS vertical eddy conductivity
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: radqrlw_eta ! longwave radiative heating rate for tp
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: radqrsw_eta ! shortwave radiative heating rate
	
	!---------- Variables collected elsewhere in WRF (unit in wrf_statistics_init) ----------
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:), PUBLIC :: wadv_eta    ! vertical advective velocity
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:), PUBLIC :: uwle_eta    ! flux from vertical advection
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:), PUBLIC :: vwle_eta    ! v
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:), PUBLIC :: twle_eta    ! theta ==> tabs
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:), PUBLIC :: mwle_eta  ! moist
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:), PUBLIC :: swle_eta  ! scalar
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:), PUBLIC :: ewle_eta    ! SGS TKE
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:), PUBLIC :: uwsb_eta    ! flux from vertical diffuison
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:), PUBLIC :: vwsb_eta    ! v
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:), PUBLIC :: wwsb_eta    ! w
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:), PUBLIC :: twsb_eta    ! theta ==> tabs
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:), PUBLIC :: mwsb_eta  ! moist
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:), PUBLIC :: swsb_eta  ! scalar
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:), PUBLIC :: ewsb_eta    ! SGS TKE
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:), PUBLIC :: flux_eta    ! temp array for tw, mw, sw
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:), PUBLIC :: e_shear_eta ! SGS TKE shear production
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:), PUBLIC :: e_buoy_eta  ! SGS TKE buoyant production
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:), PUBLIC :: e_diss_eta  ! SGS TKE dissipation
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:), PUBLIC :: precipl_eta ! liquid precip rate
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:), PUBLIC :: precipi_eta ! ice precip rate
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:), PUBLIC :: vfmi_eta    ! ice terminal velocity
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:), PUBLIC :: cond_rate_eta   ! condensation rate
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:), PUBLIC :: evap_rate_eta   ! evaporation rate
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:), PUBLIC :: deps_rate_eta   ! deposition rate
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:), PUBLIC :: subl_rate_eta   ! sublimation rate
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:), PUBLIC :: radlwup_eta ! longwave radiative upward flux
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:), PUBLIC :: radlwdn_eta ! downward
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:), PUBLIC :: radswup_eta ! shortwave upward
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:), PUBLIC :: radswdn_eta ! downward
	
	!---------- Variables from the WRF-LES module: arrays are used in module_wrfles.F ----------
	REAL, SAVE, ALLOCATABLE, DIMENSION(:), PUBLIC :: uobs, vobs, wobs, tpobs, qvobs
	REAL, SAVE, ALLOCATABLE, DIMENSION(:), PUBLIC :: unudge, vnudge, tpnudge, qvnudge, qntnudge
	REAL, SAVE, ALLOCATABLE, DIMENSION(:), PUBLIC :: tplshadv, qvlshadv
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:), PUBLIC :: ulsvadv_eta, vlsvadv_eta
	REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:), PUBLIC :: tplsvadv_eta, qvlsvadv_eta
	
	!---------- From SAM6.8.2/grid.f90 & main.f90----------
	INTEGER, SAVE :: nstep      ! current number of performed time steps
	INTEGER, SAVE :: nstat      ! interval in time step to compute statistics
	INTEGER, SAVE :: nstatis    ! interval between substeps to compute statistics = nstat/nstatfq
	INTEGER, SAVE :: nstatfrq   ! frequency of computing statistics
	INTEGER, SAVE :: nstatsteps ! index of samples over nstat. nstatfrq when mode(nstep,nstat)=0
	CHARACTER(LEN=148), SAVE, PUBLIC :: stat_file_name ! output file name, which has to be ###.stat
	INTEGER, SAVE :: rank       ! rank of the current subdomain task (default 0)
	LOGICAL, SAVE :: dompi      ! logical switch to do mutitasking
	LOGICAL, SAVE :: masterproc ! .TRUE. if rank == 0
	INTEGER, SAVE :: mpi_comm   ! communicator
	LOGICAL, SAVE :: docldconditionals     ! liquid & ice cloud
	LOGICAL, SAVE :: docndconditionals     ! liquid & ice condensate
	LOGICAL, SAVE :: docoreconditionals    ! up & down draft core
	LOGICAL, SAVE :: dosatupdnconditionals ! saturated up & down, unsaturated environment
	LOGICAL, SAVE :: doupdnconditionals    ! updraft & downdraft
	LOGICAL, SAVE :: LES    ! use different sampling criteria than non-LES run for conditional stat
	LOGICAL, SAVE :: CSRM   ! true if LES = .FALSE.
	
	!---------- From SAM6.8.2/vars.f90 ----------
	! Quantities sampled for 1D time series output
	REAL, SAVE :: w_max, w_min, u_max, u_min, v_max, v_min
	REAL, SAVE :: s_sst, s_psfc
	REAL, SAVE :: s_acld, s_acldl, s_acldm, s_acldh, s_acldcold   ! shaded cloud fraction
	REAL, SAVE :: s_ar, s_prec, s_precl, s_preci   ! precipitation
	REAL, SAVE :: s_pw, s_cwp, s_iwp, s_rwp, s_swp, s_gwp, s_lwp2, s_iwp2   ! column waters
	REAL, SAVE :: s_lhf, s_shf   ! surface flux
	REAL, SAVE :: s_flns, s_flnt, s_flnsc, s_flntc, s_flds, s_fsns, s_fsnt, s_fsnsc, s_fsntc, &
	              s_fsds, s_solin   ! radiation
	
	! Conditional average statistics
	INTEGER, SAVE :: ncondavg
	INTEGER, SAVE :: icondavg_cld
	INTEGER, SAVE :: icondavg_cldl, icondavg_cldi
	INTEGER, SAVE :: icondavg_cndl, icondavg_cndi
	INTEGER, SAVE :: icondavg_cor, icondavg_cordn
	INTEGER, SAVE :: icondavg_satdn, icondavg_satup, icondavg_env
	INTEGER, SAVE :: icondavg_wup, icondavg_wdn
	REAL, SAVE, ALLOCATABLE :: condavg_factor(:,:)   ! replaces cloud_factor, core_factor
	REAL, SAVE, ALLOCATABLE :: condavg_mask(:,:,:,:) ! indicator array for various conditional means
	CHARACTER(LEN= 8), SAVE, allocatable :: condavgname(:)       ! array of short names
	CHARACTER(LEN=25), SAVE, allocatable :: condavglongname(:)   ! array of long names
	
	!---------- From SAM6.8.2/hbuffer.f90 ----------
	INTEGER, PARAMETER :: hbuf_max_length = 1000
	INTEGER, SAVE :: hbuf_length
	REAL(8), SAVE, ALLOCATABLE, DIMENSION(:) :: hbuf   ! DOUBLE PRECISION
	CHARACTER(LEN=15), SAVE, DIMENSION(hbuf_max_length) :: namelist
	CHARACTER(LEN=80), SAVE, DIMENSION(hbuf_max_length) :: deflist
	CHARACTER(LEN=10), SAVE, DIMENSION(hbuf_max_length) :: unitlist
	INTEGER, SAVE, DIMENSION(hbuf_max_length) :: status
	INTEGER, SAVE, DIMENSION(hbuf_max_length) :: average_type
	
	!---------- Lower thresholds for resetting 0 for microphysics scalars ----------
	LOGICAL, SAVE :: dousednc
	REAL, SAVE :: qc_threshold, qi_threshold, qr_threshold, qs_threshold, qg_threshold
	
	!---------- Misc. ----------
	LOGICAL, SAVE :: flag_skip_after_restart  ! true only the time step at restart
	LOGICAL, SAVE, PUBLIC :: docollect_data   ! true when time to collect data outside of this module
	LOGICAL, SAVE, PUBLIC :: docollect_data_rk! true when docollect_data & the last step of RK step
	LOGICAL, SAVE, PUBLIC :: doavgprecflux   ! true to output time mean over nstatis as a sample
	REAL, SAVE :: factor_kg2mg   ! 1E-6 if the unit of number concentration is #/kg.
	REAL, SAVE :: factor_prec    ! 86400 / nstatis for precipitation flux
	
	! PUBLIC functions and subroutines
	PUBLIC :: wrf_statistics_init
	PUBLIC :: wrf_statistics
	
	CONTAINS
	
	!=====================================================================================
	
	SUBROUTINE wrf_statistics_init( grid, config_flags )
	
	! This subroutine is placed at the beginning of the subroutine wrf_run in /main/module_wrf_top.F
	
		USE module_dm, ONLY: ntasks_x, ntasks_y
		USE module_domain, ONLY: domain, get_ijk_from_grid
		USE module_model_constants, ONLY: g
		USE module_configure, ONLY : grid_config_rec_type
		USE module_state_description
		
		IMPLICIT NONE
		
		! Input data
		TYPE(domain), INTENT(IN) :: grid
		TYPE(grid_config_rec_type), INTENT(IN) :: config_flags   ! contains namelist data
		
		! Grid indices
		INTEGER :: ids, ide, jds, jde, kds, kde   ! Index in terms of the whole domain
		INTEGER :: ims, ime, jms, jme, kms, kme   ! Index including ghost points
		INTEGER :: ips, ipe, jps, jpe, kps, kpe   ! Index inside of the subdomain
		! ids = 1
		! ide = e_we in namelist
		! ips = 1 + e_we / #_subdomain_in_x * index_subdomain
		! ipe = ips + e_we / #_subdomain_in_x - 1 (+1 if e_we is odd)
		! Note: index_subdomain=[0,1,2,...] should be function of rank
		!       ips = ids at the most west subdomain
		!       ipe = ide at the most east subdomain
		! ims = ips - #_ghost_points
		! ime = ipe + #_ghost_points
		!
		! u(ids) = u(ide)
		! s(ids) = s(ide)
		! Because of the above relationship, in order to compute average first in each processor,
		! then MPI, ide is required
		! (1) ide has to be odd number if number of subdomains (processors) in x direction is even
		!     number.
		! (2) ide can be either odd or even number if number of subdomains in x direction is odd
		!     number.
		! Similarly in y direction
		!
		! Example 1: ide = odd (13), nproc = even (4)
		!
		!   | pr1 | pr2 | pr3 | pr4 |
		!   -+-+-+-+-+-+-+-+-+-+-+-+-+ number of grid point in each subdomain is 3 (odd)
		!   us                      us
		!   11                      11
		!
		! Example 2: ide = odd (17), nproc = even (4)
		!
		!   |  pr1  |  pr2  |  pr3  |  pr4  |
		!   -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ number of grid point in each subdomain is 4 (even)
		!   us                              us
		!   11                              11
		!
		! Example 3: ide = even (10), nproc = odd (3)
		!   | pr1 | pr2 | pr3 |
		!   -+-+-+-+-+-+-+-+-+-+ number of grid point in each subdomein is 3 (odd)
		!   us                us
		!   11                11
		!
		! Example 4: ide = odd (13), nproc = odd (3)
		!   |  pr1  |  pr2  |  pr3  |
		!   -+-+-+-+-+-+-+-+-+-+-+-+-+ number of grid point in each subdomain is 4 (even)
		!   us                      us
		!   11                      11
		!
		! - : u grid point
		! + : s grid point
		! 1 : holds same value
		
		! Misc
		REAL :: dzs, dzh
		INTEGER :: k, restart_interval_step, ierr, k1, k2
		
		! Set dompi, masterproc, rank
		CALL task_init
		PRINT*, '========== WRF_STATISTICS_INIT =========='
		PRINT*,'nsubdomains: ', nsubdomains
		
		! Time step
		dt = grid%dt
		PRINT*,'dt: ', dt
		
		! Horizontal grid spacing
		dx = grid%dx
		dy = grid%dy
		PRINT*,'dx, dy: ', dx, dy
		
		! Get grid indices
		CALL get_ijk_from_grid( grid, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, &
		                        ips, ipe, jps, jpe, kps, kpe )
		! Even number subdomain in x (y) direction ==> ide (jde) has to be odd number
		IF ( nsubdomains > 1 .AND. mod((ide-1)/(ipe-ips+1),2) == 0 .AND. mod(ide,2) == 0 ) THEN
			PRINT*,'Set e_we odd number. Abort!'
			CALL task_abort()
		ENDIF
		IF ( nsubdomains > 1 .AND. mod((jde-1)/(jpe-jps+1),2) == 0 .AND. mod(jde,2) == 0 ) THEN
			PRINT*,'Set e_sn odd number. Abort!'
			CALL task_abort()
		ENDIF
		! Check if number of grid points in each direction is divisible without remainder
		IF ( nsubdomains > 1 .AND. MOD(ide-1,ntasks_x) /= 0 ) THEN
			PRINT*,'Set e_we-1 divisible without remainder by ntaskx_x(=nproc_x). Abort!'
			CALL task_abort()
		ENDIF
		IF ( nsubdomains > 1 .AND. MOD(jde-1,ntasks_y) /= 0 ) THEN
			PRINT*,'Set e_sn-1 divisible without remainder by ntasks_y(=nproc_y). Abort!'
			CALL task_abort()
		ENDIF
		
		!==================== VERY IMPORTANT: GRID INDEX INFORMATION ====================
		i0 = ips - 1
		j0 = jps - 1
		nx = ipe - ips + 1
		ny = jpe - jps + 1
		IF ( ipe == ide ) nx = nx - 1   ! at the most east subdomain
		IF ( jpe == jde ) ny = ny - 1   ! at the most north subdomain
		nxp = nx + 1
		nyp = ny + 1
		nz  = kpe
		nzm = kpe - 1
		PRINT*,'rank, nx, ny, nz, nzm: ', rank, nx, ny, nz, nzm
		!==================================================================================
		
		
		! Hieght levels and corresponding initial pressure levels
		! See the subroutine phy_prep in /dyn_em/module_big_step_utilities_em.F
	!!	ALLOCATE( zs(nz), zw(nz), press(nz) )
	!!	DO k = 1, nz
	!!		zw(k) = grid%phb(ips,k,jps) / g
	!!	ENDDO
	!!	DO k = 1, nzm
	!!		zs(k) = 0.5 * ( zw(k) + zw(k+1) )
	!!	ENDDO
		ALLOCATE( zs(nz), press(nz) )
		DO k = 1, nzm
			zs(k) = 0.5 * ( grid%phb(ips,k,jps) + grid%phb(ips,k+1,jps) ) / g
		ENDDO
		DO k = 1, nz
			press(k) = grid%pb(ips,k,jps)
		ENDDO
		press(:) = 0.01 * press(:)   ! Pa ==> mb
	!	PRINT*, 'k, zs, zw, press:'
	!	DO k = nzm, 1, -1
	!		WRITE(*,'(I5,4F11.3,2F8.3)') k, zs(k), zw(k), press(k)
	!	ENDDO
		
		! Setup conditional averages as in setparm.f90
		CALL setparm
		
		! High resolution height level for perturbation calculation
		! zsh(1) = zs(1)
		! ...
		! zsh(1+1*hrlayers) = zs(2)
		! ...
		! zsh(1+2*hrlayers) = zs(3)
		! zsh(nzsh) = zs(nzm)
		nzsh = (nzm-1) * hrlayers + 1
		ALLOCATE( zsh(nzsh) )
		k2 = 0
		DO k = 1, nzm-1
			dzs = zs(k+1) - zs(k)
			dzh = dzs / REAL(hrlayers)
			DO k1 = 1, hrlayers
				k2 = k2 + 1
				zsh(k2) = zs(k) + REAL(k1-1) * dzh
			ENDDO
		ENDDO
		zsh(nzsh) = zs(nzm)
	!	DO k = 1, nzsh
	!		k1 = 1 + (k-1) / hrlayers
	!		WRITE(*,'(I5,2F11.3,I5)'), k, zsh(k), zs(k1), k1
	!	ENDDO
		
		! Check if mod(restart_interval_step,nstat) == 0
		! No restart file is required since output is written when the WRF restart file is written.
		restart_interval_step = REAL( grid%restart_interval * 60 ) / dt
		IF ( mod(restart_interval_step,nstat) /= 0 ) THEN
			PRINT*,'Error: mod(restart_interval_step,nstat) /= 0'
			PRINT*,'restart_interval_step, nstat: ', restart_interval_step, nstat
			CALL task_abort()
		ENDIF
		
		! Set nstep
		nstep = grid%itimestep
		PRINT*,'nstep: ',nstep
		
		! Read list of vertical profile names to store in file
		ALLOCATE( hbuf( hbuf_max_length*nzm+30 ) )   ! extra 30 for 1D time series
		CALL hbuf_init
		
		! Restart run or not
		nstatis = nstat / nstatfrq
		nstatsteps = 0
		IF ( .NOT.grid%restart ) THEN
			flag_skip_after_restart = .FALSE.   ! always false except the time step at restart
			nstep = nstep - 1   ! add 1 at the beginning of wrf_statistics
		ELSE
			flag_skip_after_restart = .TRUE.
		ENDIF
		
		! Allocate arrays (these variables are on the eta level)
		ALLOCATE( &
			! Arrays starting (x,y)=(0,0) for U, V UW, VW points (density calculation)
			! zw_eta has +1 level for zs_eta calculation
			! qn# are included here for adjustment in q#_threshold
			zw_eta(0:nx,nz,0:ny), &
			rhod_eta(0:nx,nzm,0:ny), rhom_eta(0:nx,nzm,0:ny), ps_eta(0:nx,nzm,0:ny), &
			! Normal array
			u_eta(nx,nzm,ny),   v_eta(nx,nzm,ny),    w_eta(nx,nzm,ny), &
			tp_eta(nx,nzm,ny),  qv_eta(nx,nzm,ny),   qc_eta(nx,nzm,ny),  qi_eta(nx,nzm,ny), &
			qr_eta(nx,nzm,ny),  qs_eta(nx,nzm,ny),   qg_eta(nx,nzm,ny),  qnc_eta(nx,nzm,ny), &
			qni_eta(nx,nzm,ny), qnr_eta(nx,nzm,ny),  qns_eta(nx,nzm,ny), qng_eta(nx,nzm,ny), &
			qnn_eta(nx,nzm,ny), eta_eta(nx,nzm,ny),  ssr_eta(nx,nzm,ny), &
			tke_eta(nx,nzm,ny), tkm_eta(nx,nzm,ny),  tkh_eta(nx,nzm,ny), &
			radqrlw_eta(nx,nzm,ny), radqrsw_eta(nx,nzm,ny), &
			! Output collected elsewhere in WRF
			wadv_eta(ims:ime,kms:kme,jms:jme), &
			uwle_eta(ims:ime,kms:kme,jms:jme), uwsb_eta(ims:ime,kms:kme,jms:jme), &
			vwle_eta(ims:ime,kms:kme,jms:jme), vwsb_eta(ims:ime,kms:kme,jms:jme), &
			wwsb_eta(ims:ime,kms:kme,jms:jme), &
			twle_eta(ims:ime,kms:kme,jms:jme), twsb_eta(ims:ime,kms:kme,jms:jme), &
			ewle_eta(ims:ime,kms:kme,jms:jme), ewsb_eta(ims:ime,kms:kme,jms:jme), &
			mwle_eta(ims:ime,kms:kme,jms:jme,num_moist), &
			mwsb_eta(ims:ime,kms:kme,jms:jme,num_moist), &
			swle_eta(ims:ime,kms:kme,jms:jme,num_scalar), &
			swsb_eta(ims:ime,kms:kme,jms:jme,num_scalar), &
			flux_eta(ims:ime,kms:kme,jms:jme), &
			e_shear_eta(ims:ime,kms:kme,jms:jme), e_buoy_eta(ims:ime,kms:kme,jms:jme), &
			e_diss_eta(ims:ime,kms:kme,jms:jme), &
			precipl_eta(ims:ime,kms:kme,jms:jme), precipi_eta(ims:ime,kms:kme,jms:jme), &
			vfmi_eta(ims:ime,kms:kme,jms:jme), &
			cond_rate_eta(ims:ime,kms:kme,jms:jme), evap_rate_eta(ims:ime,kms:kme,jms:jme), &
			deps_rate_eta(ims:ime,kms:kme,jms:jme), subl_rate_eta(ims:ime,kms:kme,jms:jme), &
			radlwup_eta(ims:ime,kms:kme,jms:jme), radlwdn_eta(ims:ime,kms:kme,jms:jme), &
			radswup_eta(ims:ime,kms:kme,jms:jme), radswdn_eta(ims:ime,kms:kme,jms:jme), &
			! Output from the WRF-LES module
			uobs(nzm), vobs(nzm), wobs(nzm), tpobs(nzm), qvobs(nzm), &
			unudge(nzm), vnudge(nzm), tpnudge(nzm), qvnudge(nzm), qntnudge(nzm), &
			tplshadv(nzm), qvlshadv(nzm), &
			ulsvadv_eta(ims:ime,kms:kme,jms:jme), vlsvadv_eta(ims:ime,kms:kme,jms:jme), &
			tplsvadv_eta(ims:ime,kms:kme,jms:jme), qvlsvadv_eta(ims:ime,kms:kme,jms:jme), &
			! Allocated?
			STAT=ierr )
		IF (ierr /= 0) THEN
			PRINT*,'Error: Could not allocate arrays at wrf_statistics_init!'
			CALL task_abort()
		ENDIF
		
		! Initialization of arrays, which are filled not in this module but elsewhere in WRF
		wadv_eta(:,:,:) = 0.   ! Pa/s, later converted to m/s
		uwle_eta(:,:,:) = 0.   ! Pa/s m/s
		vwle_eta(:,:,:) = 0.   ! Pa/s m/s
		twle_eta(:,:,:) = 0.   ! K Pa/s for TP
		mwle_eta(:,:,:,:) = 0.   ! kg/kg Pa/s
		swle_eta(:,:,:,:) = 0.   ! scalar_unit Pa/s
		ewle_eta(:,:,:) = 0.   ! m3/s3
		uwsb_eta(:,:,:) = 0.   ! m2/s2
		vwsb_eta(:,:,:) = 0.   ! m2/s2
		wwsb_eta(:,:,:) = 0.   ! m2/s2
		twsb_eta(:,:,:) = 0.   ! K m/s for TP
		mwsb_eta(:,:,:,:) = 0.   ! kg/kg m/s
		swsb_eta(:,:,:,:) = 0.   ! scalar_unit m/s
		ewsb_eta(:,:,:) = 0.   ! m3/s3
		flux_eta(:,:,:) = 0.
		e_shear_eta(:,:,:) = 0.   ! m2/s3
		e_buoy_eta(:,:,:) = 0.    ! m2/s3
		e_diss_eta(:,:,:) = 0.    ! m2/s3
		! /phys/module_mp_###.F
		precipl_eta(:,:,:) = 0.   ! kg/kg m/s = mm/s
		precipi_eta(:,:,:) = 0.   ! kg/kg m/s = mm/s
		vfmi_eta(:,:,:) = 0.   ! m/s
		cond_rate_eta(:,:,:) = 0.   ! kg/kg/s
		evap_rate_eta(:,:,:) = 0.   ! kg/kg/s
		deps_rate_eta(:,:,:) = 0.   ! kg/kg/s
		subl_rate_eta(:,:,:) = 0.   ! kg/kg/s
		! /phys/module_ra_###.F
		radlwup_eta(:,:,:) = 0.   ! W/m2 for TABS
		radlwdn_eta(:,:,:) = 0.   ! W/m2 for TABS
		radswup_eta(:,:,:) = 0.   ! W/m2 for TABS
		radswdn_eta(:,:,:) = 0.   ! W/m2 for TABS
		! /frame/module_wrfles.F
		uobs(:)  = 0.   ! m/s
		vobs(:)  = 0.   ! m/s
		wobs(:)  = 0.   ! m/s
		tpobs(:) = 0.   ! K for TP
		qvobs(:) = 0.   ! kg/kg
		unudge(:)  = 0.   ! m/s2
		vnudge(:)  = 0.   ! m/s2
		tpnudge(:) = 0.   ! K/s for TP
		qvnudge(:) = 0.   ! kg/kg/s
		qntnudge(:) = 0.   ! #/mg/s
		tplshadv(:) = 0.   ! K/s for TP
		qvlshadv(:) = 0.   ! kg/kg/s
		ulsvadv_eta(:,:,:)  = 0.   ! m/s2
		vlsvadv_eta(:,:,:)  = 0.   ! m/s2
		tplsvadv_eta(:,:,:) = 0.   ! K/s for TP
		qvlsvadv_eta(:,:,:) = 0.   ! kg/kg/s
		
		
		! Set unit conversion factor.
		! 1. factor_kg2mg for number concentration depending on microphysics scheme
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Add other microphysics scheme if factor_kg2mg /= 1.E-6 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		SELECT CASE( config_flags%mp_physics )
#ifdef ESRL
			CASE (FEINGOLD2M)
				factor_kg2mg = 1.
#endif /*ESRL*/
			CASE DEFAULT
				factor_kg2mg = 1.E-6
		END SELECT
		
		! 2. factor_prec
		! precipitation flux is accumulated over nstatis = nstat / nstatfrq
		! 86400 for conversion from mm/s to mm/day
		IF ( doavgprecflux ) THEN
			factor_prec = 86400. / REAL(nstatis)
		ELSE
			factor_prec = 86400.
		ENDIF
		
		PRINT*,'========================================='
		
	END SUBROUTINE wrf_statistics_init
	
	!=====================================================================================
	
	SUBROUTINE wrf_statistics( grid, config_flags )
	
	! This subroutine is placed right after the call of med_before_solve_io and med_last_solve_io.
	
		USE module_domain, ONLY: domain
		USE module_nesting, ONLY: nests_to_open
		USE module_configure, ONLY : grid_config_rec_type
		USE module_state_description
		USE module_model_constants, ONLY: t0, g
		
		IMPLICIT NONE
		
		! Input data
		TYPE(domain), INTENT(IN) :: grid
		TYPE(grid_config_rec_type), INTENT(IN) :: config_flags   ! contains namelist data
		
		! Local
!?????		INTEGER :: nestid, kid !???????????????????????????????????????????????????????????????
		INTEGER :: i, j, k, iwrf, jwrf
		REAL :: dqnn
		LOGICAL :: flag_nesting
		
!?????		IF ( .NOT.nests_to_open( grid, nestid, kid ) ) THEN !??????????????????????????????????
!????? Is this right way to avoid nested domain????????????????????????????????????????????????????
		flag_nesting = .FALSE.   ! temporal
		IF ( (.NOT.flag_nesting) .AND. (.NOT.flag_skip_after_restart) ) THEN
			nstep = nstep + 1
		!	PRINT *,'wrf_statistics: nstep, grid%itimestep = ', nstep, grid%itimestep
			
			! Diagnose some 1D time series output, which are reuired to collect every time step
			! Find w_max, w_min, u_max, u_min, v_max, v_min: code from kurant.f90
			w_max = MAX( w_max, MAXVAL( grid%w_2(i0+1:i0+nx,1:nz, j0+1:j0+ny) ) )   ! m/s
			w_min = MIN( w_max, MINVAL( grid%w_2(i0+1:i0+nx,1:nz, j0+1:j0+ny) ) )   ! m/s
			u_max = MAX( u_max, MAXVAL( grid%u_2(i0+1:i0+nx,1:nzm,j0+1:j0+ny) ) )   ! m/s
			u_min = MIN( u_min, MINVAL( grid%u_2(i0+1:i0+nx,1:nzm,j0+1:j0+ny) ) )   ! m/s
			v_max = MAX( v_max, MAXVAL( grid%v_2(i0+1:i0+nx,1:nzm,j0+1:j0+ny) ) )   ! m/s
			v_min = MIN( v_min, MINVAL( grid%v_2(i0+1:i0+nx,1:nzm,j0+1:j0+ny) ) )   ! m/s
			
			! Time to collect data?
			! The flag docollect_data tells if time to collect data in other source files.
			! Collection has to be done one time step before.
			IF ( mod(nstep+1,nstatis) == 0 ) THEN
				docollect_data = .TRUE.
			!	PRINT*,'Collecting data...'
			ELSE
				docollect_data = .FALSE.
			ENDIF
			docollect_data_rk = .FALSE.   ! Set TRUE in solve_em.F if condition is met.
			
			! Time to compute and save statistics?
			IF ( mod(nstep,nstatis) == 0 ) THEN
				nstatsteps = nstatsteps + 1
				PRINT*,'Collecting statistics...'
				
				! Pass WRF data in "grid" to arrays
				! Height and pressure conversion is in the subroutine phy_prep in
				! /dyn_em/module_big_step_utilities_em.F
				! Height at W: m
				DO j = 0, ny   ! 0 for VW
					jwrf = j0 + j
					DO k = 1, nz
						DO i = 0, nx   ! 0 for UW
							iwrf = i0 + i
							zw_eta(i,k,j) = ( grid%phb(iwrf,k,jwrf) + grid%ph_2(iwrf,k,jwrf) ) / g
						ENDDO
					ENDDO
				ENDDO
				! Pressure at S: Pa => mb (=hPa)
				DO j = 0, ny   ! 0 for V
					jwrf = j0 + j
					DO k = 1, nzm
						DO i = 0, nx   ! 0 for U
							iwrf = i0 + i
							ps_eta(i,k,j) = 0.01 * grid%p(iwrf,k,jwrf) + press(k)
						ENDDO
					ENDDO
				ENDDO 
				! Air density at S: kg/m3
				! alb and al are dry air density
				DO j = 0, ny   ! 0 for V
					jwrf = j0 + j
					DO k = 1, nzm
						DO i = 0, nx   ! 0 for U
							iwrf = i0 + i
							rhod_eta(i,k,j) = 1. / ( grid%alb(iwrf,k,jwrf) + grid%al(iwrf,k,jwrf) )
							rhom_eta(i,k,j) = 1. / ( grid%alb(iwrf,k,jwrf) + grid%al(iwrf,k,jwrf) ) &
							        * ( 1. + grid%moist(iwrf,k,jwrf,P_qv) + grid%moist(iwrf,k,jwrf,P_qc) &
							               + grid%moist(iwrf,k,jwrf,P_qi) + grid%moist(iwrf,k,jwrf,P_qr) &
							               + grid%moist(iwrf,k,jwrf,P_qs) + grid%moist(iwrf,k,jwrf,P_qg) )
						ENDDO
					ENDDO
				ENDDO
				! U, V, W, & Otehr
				DO j = 1, ny
					jwrf = j0 + j
					DO k = 1, nzm
						DO i = 1, nx
							iwrf = i0 + i
							! U, V, W
							u_eta(i,k,j) = grid%u_2(iwrf,k,jwrf)   ! m/s
							v_eta(i,k,j) = grid%v_2(iwrf,k,jwrf)   ! m/s
							w_eta(i,k,j) = grid%w_2(iwrf,k,jwrf)   ! m/s
							! TP
							tp_eta(i,k,j) = grid%t_2(iwrf,k,jwrf) + t0   ! K-t0 => K
							! QV, QC, QI, QR, QS, QG
							qv_eta(i,k,j) = grid%moist(iwrf,k,jwrf,P_qv)   ! kg/kg
							qc_eta(i,k,j) = grid%moist(iwrf,k,jwrf,P_qc)   ! kg/kg
							qi_eta(i,k,j) = grid%moist(iwrf,k,jwrf,P_qi)   ! kg/kg
							qr_eta(i,k,j) = grid%moist(iwrf,k,jwrf,P_qr)   ! kg/kg
							qs_eta(i,k,j) = grid%moist(iwrf,k,jwrf,P_qs)   ! kg/kg
							qg_eta(i,k,j) = grid%moist(iwrf,k,jwrf,P_qg)   ! kg/kg
							! ETA, SSR
#ifdef ESRL
! These are only for the Feingold scheme
							eta_eta(i,k,j) = grid%scalar(iwrf,k,jwrf,P_eta0)  ! kg/kg
							ssr_eta(i,k,j) = grid%scalar(iwrf,k,jwrf,P_ssr)   ! %
#else
							eta_eta(i,k,j) = 0.
							ssr_eta(i,k,j) = 0.
#endif /*ESRL*/
							! QNC, QNI, QNR, QNS, QNG, QNN
							qnc_eta(i,k,j) = grid%scalar(iwrf,k,jwrf,P_qnc) * factor_kg2mg ! #/kg => #/mg
							qni_eta(i,k,j) = grid%scalar(iwrf,k,jwrf,P_qni) * factor_kg2mg ! #/kg => #/mg
							qnr_eta(i,k,j) = grid%scalar(iwrf,k,jwrf,P_qnr) * factor_kg2mg ! #/kg => #/mg
							qns_eta(i,k,j) = grid%scalar(iwrf,k,jwrf,P_qns) * factor_kg2mg ! #/kg => #/mg
							qng_eta(i,k,j) = grid%scalar(iwrf,k,jwrf,P_qng) * factor_kg2mg ! #/kg => #/mg
							qnn_eta(i,k,j) = grid%scalar(iwrf,k,jwrf,P_qnn) * factor_kg2mg ! #/kg => #/mg
							! PRECIPL, PRECIPI
							precipl_eta(iwrf,k,jwrf) = precipl_eta(iwrf,k,jwrf) * factor_prec! mm/s=>mm/day
							precipi_eta(iwrf,k,jwrf) = precipi_eta(iwrf,k,jwrf) * factor_prec! mm/s=>mm/day
							! SGS TKE
							tke_eta(i,k,j) = grid%tke_2(iwrf,k,jwrf)   ! m2/s2
							! TKM & TKH
							tkm_eta(i,k,j) = grid%xkmv(iwrf,k,jwrf)   ! m2/s
							tkh_eta(i,k,j) = grid%xkhv(iwrf,k,jwrf)   ! m2/s
							! RADLWQR, RADSWQR
							radqrlw_eta(i,k,j) = grid%rthratenlw(iwrf,k,jwrf)   ! K/s for TP
							radqrsw_eta(i,k,j) = grid%rthratensw(iwrf,k,jwrf)   ! K/s for TP
						ENDDO
					ENDDO
				ENDDO
				
				IF (LES) THEN
				! Reset qc_eta, qi_eta, qr_eta, qs_eta, qg_eta, qnc_eta, qni_eta, qnr_eta, qns_eta,
				! qng_eta to 0 if the grid value is smaller than the specified threshold.
				! Also add number concentration back to qnn to make total number unchanged.
				! (default = 0 kg/kg or 0 /mg if dousednc)
					IF (.NOT.dousednc) THEN
						DO j = 1, ny
							DO k = 1, nzm
								DO i = 1, nx
									dqnn = 0.
									IF ( qc_eta(i,k,j) < qc_threshold ) THEN
										qc_eta(i,k,j) = 0.
										dqnn = dqnn + qnc_eta(i,k,j)
										qnc_eta(i,k,j) = 0.
									ENDIF
									IF ( qi_eta(i,k,j) < qi_threshold ) THEN
										qi_eta(i,k,j) = 0.
										dqnn = dqnn + qni_eta(i,k,j)
										qni_eta(i,k,j) = 0.
									ENDIF
									IF ( qr_eta(i,k,j) < qr_threshold ) THEN
										qr_eta(i,k,j) = 0.
										dqnn = dqnn + qnr_eta(i,k,j)
										qnr_eta(i,k,j) = 0.
									ENDIF
									IF ( qs_eta(i,k,j) < qs_threshold ) THEN
										qs_eta(i,k,j) = 0.
										dqnn = dqnn + qns_eta(i,k,j)
										qns_eta(i,k,j) = 0.
									ENDIF
									IF ( qg_eta(i,k,j) < qg_threshold ) THEN
										qg_eta(i,k,j) = 0.
										dqnn = dqnn + qng_eta(i,k,j)
										qng_eta(i,k,j) = 0.
									ENDIF
									qnn_eta(i,k,j) = qnn_eta(i,k,j) + dqnn
								ENDDO
							ENDDO
						ENDDO
					ELSE
						DO j = 1, ny
							DO k = 1, nzm
								DO i = 1, nx
									dqnn = 0.
									IF ( qnc_eta(i,k,j) < qc_threshold ) THEN
										qc_eta(i,k,j) = 0.
										dqnn = dqnn + qnc_eta(i,k,j)
										qnc_eta(i,k,j) = 0.
									ENDIF
									IF ( qni_eta(i,k,j) < qi_threshold ) THEN
										qi_eta(i,k,j) = 0.
										dqnn = dqnn + qni_eta(i,k,j)
										qni_eta(i,k,j) = 0.
									ENDIF
									IF ( qnr_eta(i,k,j) < qr_threshold ) THEN
										qr_eta(i,k,j) = 0.
										dqnn = dqnn + qnr_eta(i,k,j)
										qnr_eta(i,k,j) = 0.
									ENDIF
									IF ( qns_eta(i,k,j) < qs_threshold ) THEN
										qs_eta(i,k,j) = 0.
										dqnn = dqnn + qns_eta(i,k,j)
										qns_eta(i,k,j) = 0.
									ENDIF
									IF ( qng_eta(i,k,j) < qg_threshold ) THEN
										qg_eta(i,k,j) = 0.
										dqnn = dqnn + qng_eta(i,k,j)
										qng_eta(i,k,j) = 0.
									ENDIF
									qnn_eta(i,k,j) = qnn_eta(i,k,j) + dqnn
								ENDDO
							ENDDO
						ENDDO
					ENDIF   ! .NOT.dousednc
				ENDIF   ! LES
				
				! 1D time series output
				DO j = 1, ny
					jwrf = j0 + j
					DO i = 1, nx
						iwrf = i0 + i
						! SST
						s_sst = s_sst + grid%sst(iwrf,jwrf)   ! K
						! Surface pressure
						s_psfc = s_psfc + grid%psfc(iwrf,jwrf) * 0.01   ! Pa => hPa (mb)
						! Surface sensible and latent heat flux
						s_shf = s_shf + grid%hfx(iwrf,jwrf)   ! W/m2
						s_lhf = s_lhf + grid%lh(iwrf,jwrf)    ! W/m2
						! Surface precipitation
						s_prec = s_prec + precipl_eta(iwrf,1,jwrf) + precipi_eta(iwrf,1,jwrf)   ! mm/day
						s_precl = s_precl + precipl_eta(iwrf,1,jwrf)   ! mm/day
						s_preci = s_preci + precipi_eta(iwrf,1,jwrf)   ! mm/day
						! Surface precipitation fraction
						IF ( precipl_eta(iwrf,1,jwrf)+precipi_eta(iwrf,1,jwrf) > 0.1 ) s_ar = s_ar + 1.
					ENDDO
				ENDDO
				
				! Compute statistics
				CALL statistics
				
				! Save statistics
				IF ( mod(nstep,nstat) == 0 ) THEN
					PRINT *,'Writing statistics:nstatstep=',nstatsteps
					CALL hbuf_average()
					CALL hbuf_write()
					CALL hbuf_flush()
					nstatsteps = 0
					PRINT *,'Finish writing statistics'
				ENDIF
				
				! Initialize precipl_eta, precipi_eta for next sample period.
				! Important if doavgprecflux=.TRUE.
				precipl_eta(:,:,:) = 0.
				precipi_eta(:,:,:) = 0.
				
			ENDIF   ! if compute statistics
		ENDIF   ! if not nesting & not after restart
		
		! Make sure the restart flag is false, now
		flag_skip_after_restart = .FALSE.
	
	END SUBROUTINE wrf_statistics
	
	!==============================================================================================
	!************************* SUBROUTINES BELOW DO NOT INTERACT WITH WRF *************************
	!==============================================================================================
	
	SUBROUTINE setparm
	
	! Initialize parameters
	
		IMPLICIT NONE
		
#ifdef WRFLES
		CHARACTER (LEN=50) :: casename
#endif /*WRFLES*/
		INTEGER :: icondavg, ierr
		
		NAMELIST /STAT/ &
			stat_file_name, nstat, nstatfrq, LES, doavgprecflux, &
			docldconditionals, docndconditionals, docoreconditionals, dosatupdnconditionals, &
			doupdnconditionals, &
			hrlayers, dousednc, qc_threshold, qi_threshold, qr_threshold, qs_threshold, qg_threshold
		
		! Set defaults
		stat_file_name = 'wrf.stat'
		nstat    = 1
		nstatfrq = 1
		CSRM = .TRUE.
		LES = .FALSE.
		doavgprecflux = .FALSE.
		docldconditionals  = .FALSE.
		docndconditionals  = .FALSE.
		docoreconditionals = .FALSE.
		dosatupdnconditionals = .FALSE.
		doupdnconditionals = .FALSE.
		hrlayers = 10
		dousednc = .FALSE.
		qc_threshold = 0.   ! g/kg or #/mg
		qi_threshold = 0.
		qr_threshold = 0.
		qs_threshold = 0.
		qg_threshold = 0.
		
		! Read namelist variables from prm
#ifdef WRFLES
		OPEN(66, FILE='./CaseName', STATUS='old', FORM='formatted')
		READ(66,*) casename
		CLOSE(66)
		OPEN(55, FILE='./'//TRIM(casename)//'/prm', STATUS='old', FORM='formatted')
#else
		OPEN(55, FILE='./stat/prm', STATUS='old', FORM='formatted')
#endif /*WRFLES*/
		READ(55, STAT)
		CLOSE(55)
#ifdef WRFLES
		stat_file_name = TRIM(casename)//'_'//TRIM(stat_file_name)
#endif /*WRFLES*/
		IF (LES) CSRM = .FALSE.
		PRINT*,'Output file name: ', TRIM(stat_file_name)
		PRINT*,'nstat: ', nstat
		PRINT*,'nstatfrq: ', nstatfrq
		PRINT*,'CSRM: ', CSRM
		PRINT*,'LES: ', LES
		PRINT*,'doavgprecflux:', doavgprecflux
		PRINT*,'docldconditionals:', docldconditionals
		PRINT*,'docndconditionals:', docndconditionals
		PRINT*,'docoreconditionals: ', docoreconditionals
		PRINT*,'dosatupdnconditionals: ', dosatupdnconditionals
		PRINT*,'doupdnconditionals: ', doupdnconditionals
		PRINT*,'hrlayers: ', hrlayers
		PRINT*,'qc_threshold: ', qc_threshold
		PRINT*,'qi_threshold: ', qi_threshold
		PRINT*,'qr_threshold: ', qr_threshold
		PRINT*,'qs_threshold: ', qs_threshold
		PRINT*,'qg_threshold: ', qg_threshold
		PRINT*,'dousednc: ', dousednc
		
		! hrlayers >= 1 and even number
		IF ( hrlayers == 0 ) THEN
			WRITE(*,*) 'Error: hrlayers has to be >= 1!'
			CALL task_abort()
		ENDIF
		
		! Convert unit of thresholds from g/kg to kg/kg
		IF (.NOT.dousednc) THEN
			qc_threshold = qc_threshold * 1.E-3
			qi_threshold = qi_threshold * 1.E-3
			qr_threshold = qr_threshold * 1.E-3
			qs_threshold = qs_threshold * 1.E-3
			qg_threshold = qg_threshold * 1.E-3
		ENDIF
		
		! bloss: set up conditional averages
		ncondavg = 1 ! always output CLD conditional average
		IF (docldconditionals) ncondavg = ncondavg + 2
		IF (docndconditionals) ncondavg = ncondavg + 2
		IF (docoreconditionals) ncondavg = ncondavg + 2
		IF (dosatupdnconditionals) ncondavg = ncondavg + 3
		IF (doupdnconditionals) ncondavg = ncondavg + 2
		ALLOCATE(condavg_factor(nzm,ncondavg), &   ! replaces old cloud_factor, core_factor
		         condavg_mask(nx,nzm,ny,ncondavg), &   ! nx x ny x nzm indicator arrays
		         condavgname(ncondavg), &   ! short names (e.g. CLD, COR, SATUP)
		         condavglongname(ncondavg), &   ! long names (e.g. cloud, core, saturated updraft)
		         STAT=ierr)
		IF (ierr /= 0) THEN
			WRITE(*,*) 'Error: Could not allocate arrays for conditional statistics!'
			CALL task_abort()
		ENDIF
		
		! indicators that can be used to tell whether a particular average is present.
		! If >0, these give the index into the condavg arrays where this particular conditional
		! average appears.
		icondavg_cld = -1
		icondavg_cldl = -1
		icondavg_cldi = -1
		icondavg_cndl = -1
		icondavg_cndi = -1
		icondavg_cor = -1
		icondavg_cordn = -1
		icondavg_satup = -1
		icondavg_satdn= -1
		icondavg_env = -1
		
		icondavg = 0
		icondavg = icondavg + 1
		condavgname(icondavg) = 'CLD'
		condavglongname(icondavg) = 'cloud'
		icondavg_cld = icondavg
		
		IF (docldconditionals) THEN
			icondavg = icondavg + 1
			condavgname(icondavg) = 'CLDL'
			condavglongname(icondavg) = 'liquid cloud'
			icondavg_cldl = icondavg
			
			icondavg = icondavg + 1
			condavgname(icondavg) = 'CLDI'
			condavglongname(icondavg) = 'ice cloud'
			icondavg_cldi = icondavg
		ENDIF
		
		IF (docndconditionals) THEN
			icondavg = icondavg + 1
			condavgname(icondavg) = 'CNDL'
			condavglongname(icondavg) = 'liquid water'
			icondavg_cndl = icondavg
			
			icondavg = icondavg + 1
			condavgname(icondavg) = 'CNDI'
			condavglongname(icondavg) = 'ice water'
			icondavg_cndi = icondavg
		ENDIF
		
		IF (docoreconditionals) THEN
			icondavg = icondavg + 1
			condavgname(icondavg) = 'COR'
			condavglongname(icondavg) = 'core'
			icondavg_cor = icondavg
			
			icondavg = icondavg + 1
			condavgname(icondavg) = 'CDN'
			condavglongname(icondavg) = 'downdraft core'
			icondavg_cordn = icondavg
		ENDIF
		
		IF (dosatupdnconditionals) THEN
			icondavg = icondavg + 1
			condavgname(icondavg) = 'SUP'
			condavglongname(icondavg) = 'saturated updrafts'
			icondavg_satup = icondavg
			
			icondavg = icondavg + 1
			condavgname(icondavg) = 'SDN'
			condavglongname(icondavg) = 'saturated downdrafts'
			icondavg_satdn = icondavg
			
			icondavg = icondavg + 1
			condavgname(icondavg) = 'ENV'
			condavglongname(icondavg) = 'unsaturated environment'
			icondavg_env = icondavg
		ENDIF
		
		IF (doupdnconditionals) THEN
			icondavg = icondavg + 1
			condavgname(icondavg) = 'WUP'
			condavglongname(icondavg) = 'updrafts'
			icondavg_wup = icondavg
			
			icondavg = icondavg + 1
			condavgname(icondavg) = 'WDN'
			condavglongname(icondavg) = 'downdrafts'
			icondavg_wdn = icondavg
		ENDIF
		
	END SUBROUTINE setparm
	
	!=====================================================================================
	
	SUBROUTINE statistics
	
	! This subroutine is based on statistics.f90, but significantly modified.
	!
	! WRF variables are on the eta level. The height on the constant eta level changes with time.
	! In order to compute horizontal mean profiles with a height level, a value at each output
	! level has to be computed with vertical interpolation from the local eta level.
	!
	! Since some variables are defined non-linearly, interpolation is taken as late as possible,
	! that is, these variables will be diagnosed on the local WRF grid before interpolated. The
	! resulting profiles is not necessary to satisfy the exact relationship like qn = qc + qi.
	!
	! A monotone cubic interpolation (Steffen, M., 1990 Astron. Astrophys.) is used.
	!
	! This subroutine is memory intense. Using 1D arrays in the innermost loop could be done...
	
	!---------- NOTE: PLACE WHERE TENDENCIES AND FLUXES CAN BE COLLECTED ----------
	! solve_em
	!   IF ( docollect_data & rk_step == rk_order ) docollect_data_rk = .TRUE., otherwise FALSE.
	!   CALL rk_step_pres
	!   <=== vertical velocity for advection can be collected here.
	!   IF ( rk_step == 1) THEN
	!      CALL first_rk_step_part1
	!      CALL first_rk_step_part2 <== diffusive tendency
	!         CALL vertical_diffusion_2 <== SGS vertical flux
	!            CALL vertical_diffusion_2_u <== SGS vertical flux
	!            CALL vertical_diffusion_2_v <== SGS vertical flux
	!            CALL vertical_diffusion_s <== SGS vertical flux, flux_eta is used
	!   ENDIF
	!   CALL rk_tendency <== momentum advective tendency if docollect_data_rk
	!      CALL advect_u <== vertical flux if docollect_data_rk
	!      CALL advect_v <== vertical flux if docollect_data_rk
	!      CALL advect_scalar <== vertical flux if docollect_data_rk, flux_eta is used
	!   CALL rk_scalar_tend <== scalar advective tendency
	!                       <== vertical flux if docollect_data_rk, flux_eta is used
	!
	! 1. ADVECTION
	! namelist: v_mom_adv_order=3, v_sca_adv_order=3
	!           moist_adv_opt=2, scalar_adv_opt=2, tke_adv_opt=2
	! Note: Other option will not produce output.
	!
	! 1.1 VERTICAL FLUX OF ADVECTION (resolved scale flux, TKE tendency)
	! /dyn_em/module_advect_em.F/advect_u, advect_v, advect_w, advect_scalar, advect_scalar_mono
	! (1) docollect_data_rk == .TRUE.?
	! (2) Grab flux only vertical advection.
	! (3) Unit conersion and coordinate transformation from eta to z
	!     The diagnostic relation for the dry density is given as
	!        d( g z )/deta = - mu / rho ===> dz/deta = - mu / ( rho g )
	!     where mu and rho represent for dry air. Usually mud=mut, rhod=1/alt in WRF
	!     Since w = dz/dt & omega = deta/dt
	!        dz/dt * dt/deta = - mu / ( rho g ) ===> dz/dt = - mu / ( rho g ) deta/dt
	!        ===> w = - mu / (rho g) omega ===> omega = - rho g w / mu
	!     Vertical flux is computed as
	!        flux = mu * omega * var = mu * ( - rho g w / mu ) * var
	!             = - rho g w var
	!        ===> w var = - flux / ( rho g )
	!     rho is dry air density.
	!     This conversion is performed in this subroutine.
	! (5) Theta flux from advect_scalar is stored from flux_eta to twle_eta after
	!     CALL advect_scalar in /dyn_em/module_em.F/rk_tendency.
	! (6) Flux from advect_scalar_mono is stored from flux_eta to the designated array after
	!     CALL rk_scalar_tend in /dyn_em/solve_em.F
	! Fluxes are defined at the w level except w flux, which is at the s level
	!
	! 1.2 ADVECTIVE TENDENCY (resolved scale turbulence advective tendency of TKE, SGS TKE)
	! Momentum & theta @ /dyn_em/module_em.F/rk_tendency
	! (1) docollect_data_rk == .TRUE.?
	! (2) ru_tend, rv_tend, rw_tend, t_tend after CALL advect_u, advect_v, advect_w, advect_scalar
	! (3) Unit conversion: ru_tend/muu, rv_tend/muv, rw_tend/mut to m/s2, & t_tend/mut to K/s
	! Moist, scalar, SGS TKE @ /dyn_em/solve_em/F
	! (1) docollect_data_rk == .TRUE.?
	! (2) advect_tend right after CALL rk_scalar_tend
	! (3) Unit conversion: advect_tend/muts, where muts is updated mut in RK steps.
	! These advective tendency requires the coordinate transformation for horizontal tendency,
	! and also correction to full density from dry density. Since horizontal mean cancels
	! horizontal advective tendency, the vertical advective tendency is computed in this
	! subroutine.
	!
	! 2. DIFFUSION
	! namelist: diff_opt=2, km_opt=2 to use SGS TKE scheme
	!           diff_opt=2, bl_pbl_physics=0 for diffusion
	!           isfflx=1 to apply WRFLSE surface fluxes
	! NOTE: km_opt=3 (Smagorinsky) works but no SGS TKE output by design.
	!
	! 2.1 VERTICAL FLUX OF DIFFUSION (SGS flux)
	! /dyn_em/module_diffusion_em.F/vertical_diffusion_2, vertical_diffusion_2_u,
	! vertical_diffusion_2_v, vertical_diffusion_2_w, vertical_dissusion_s
	! (1) Grab flux only vertical diffusion in #_2_u, #_2_v, #_2_w, #_s, #=vertical_diffusion
	! (2) Grab surface flux in vertical_diffusion_2
	! (3) Flux from vertical_diffusion_s is stored from flux_eta to the designated array after
	!     CALL verrtical_diffusion_s.
	! Unit conversion is not required since flux with correct unit can be obtained.
	! Coordinate transformation is not required since flux is computed with height coordinate.
	! Fluxes are defined at the w level except w flux, which is at the s level.
	!
	! 2.2 DIFFUSIVE TENDENCY (SGS turbulence transport of resolved scale & SGS TKE)
	! /dyn_em/module_first_rk_step_part2.F
	! (1) Store ru_tendf, rv_tendf, rw_tendf, t_tendf, tke_tend, moist_tend, scalar_tend right
	!     before CALL vertical_diffusion_2
	! (2) Subtract the stored tendency from tendency right after CALL horizontal_diffusion_2.
	! (3) Unit conversion by dividing muu, muv, mut (w and s)
	! Coordinate transformation is not required since tendency is computed with height coordinate.
	! However correction to full density from dry density has to be performed. Since horizontal
	! mean cancels horizontal diffusive tendency, the vertical diffusive tendency is computed in
	! this subroutine.
	!
	! 3. SGS TKE BUDGET
	! namelist: diff_opt=2, km_opt=2 to use SGS TKE scheme
	! /dyn_em/module_diffusion_em.F/tke_rhs
	! (1) Store tendency right before CALL tke_shear, CALL tke_buoyancy, CALL tke_dissip
	! (2) Subtract the stored tendency from tendency right after each CALL.
	! (3) Unit conversion to m2/s3 by tendency/mut
	! Coordinate transformation for SGS shear production is not required since it uses the stress
	! tensor computed for diffusion with height coordinate.
	!
	! 4. BUOYANCY (resolved scale buoyant production)
	!!!!! This need conversion from mu' to rho' and I do not know how to do it now.
	!!!!! So. buoyant production is diagnosed in this subrutine.
	! /dyn_em/module_big_step_utilities_em.F/pg_buoy_w
	! (1) docollect_data_rk == .TRUE.?
	! (2) buoyancy (m/s2, positive upward) = - g * mu / (mu + mub)
	!
	! NOTE: All map-scale factors are assumed to 1.
	!---------------------------------------------------------------------------
		
	!	USE isccp, only: isccp_get <== future option?
		USE module_state_description 
		USE module_model_constants, ONLY: &
			cp, ggr => g, lcond => xlv, lfus => xlf, lsub => xls, rv => r_v, rgas => r_d
		
		IMPLICIT NONE
		
		! Physical constants from params.f90
	!	REAL, PARAMETER :: cp = 1004.5          ! Specific heat of air, J/kg/K
	!	REAL, PARAMETER :: ggr = 9.81           ! Gravity acceleration, m/s2
	!	REAL, PARAMETER :: lcond = 2.5104e+06   ! Latent heat of condensation, J/kg
	!	REAL, PARAMETER :: lfus = 0.3336e+06    ! Latent heat of fusion, J/kg
	!	REAL, PARAMETER :: lsub = 2.8440e+06    ! Latent heat of sublimation, J/kg
	!	REAL, PARAMETER :: rv = 461.6           ! Gas constant for water vapor, J/kg/K
	!	REAL, PARAMETER :: rgas = 287.          ! Gas constant for dry air, J/kg/K
	!Unused! REAL, PARAMETER :: diffelq = 2.21e-05   ! Diffusivity of water vapor, m2/s
	!Unused!	REAL, PARAMETER :: therco = 2.40e-02    ! Thermal conductivity of air, J/m/s/K
	!Unused!	REAL, PARAMETER :: muelq = 1.717e-05    ! Dynamic viscosity of air
		REAL, PARAMETER :: fac_ggr = ggr / cp
		REAL, PARAMETER :: fac_cond = lcond/cp
		REAL, PARAMETER :: fac_fus = lfus/cp
		REAL, PARAMETER :: fac_sub = lsub/cp
		REAL, PARAMETER :: epsv = (rv-rgas)/rgas   ! ~0.61
		
		! 1D array for summation over subdomain for output defined at ZS
		! 1st moment
		REAL, DIMENSION(nzm) :: pres
		REAL, DIMENSION(nzm) :: rho   ! summation ==> mean
		REAL, DIMENSION(nzm) :: u, v
		REAL, DIMENSION(nzm) :: dse, lwse, liwse, mse, vdse, smse
		REAL, DIMENSION(nzm) :: tp, tpl, tpli, tpe, tpv, tabs
		REAL, DIMENSION(nzm) :: qt, qv, qc, qi, qr, qs, qg, qn, qp, qcl, qci, qct, qsat
		REAL, DIMENSION(nzm) :: relh, eta, ssr
		REAL, DIMENSION(nzm) :: qnc, qni, qnr, qns, qng, qnn, qnt
		REAL, DIMENSION(nzm) :: cond_rate, evap_rate, ncond_rate, deps_rate, subl_rate, &
		                        ndeps_rate, tcond_rate
		! Flux
		REAL, DIMENSION(nzm) :: uw, uwle, uwsb, vw, vwle, vwsb
		REAL, DIMENSION(nzm) :: lwsef, lwsefr, lwsefs, liwsef, liwsefr, liwsefs
		REAL, DIMENSION(nzm) :: tplf, tplfr, tplfs, tpvf, tpvfr, tpvfs
		REAL, DIMENSION(nzm) :: qtf, qtfr, qtfs, qpf, qpfr, qpfs, qcf, qcfr, qcfs, qif, qifr, qifs, &
		                        qrf, qrfr, qrfs, qsf, qsfr, qsfs, qgf, qgfr, qgfs
		REAL, DIMENSION(nzm) :: prec, precl, preci, vfmi
		! Higher-order moment
		REAL, DIMENSION(nzm) :: u2, v2, w2, w3, w4, wskw
		REAL, DIMENSION(nzm) :: liwse2, liwse3, mse2, mse3, vdse2, vdse3
		REAL, DIMENSION(nzm) :: qt2, qt3, qc2, qi2, qcl2, qci2, qm2
		REAL, DIMENSION(nzm) :: slqt, aup
		! TKE
		REAL, DIMENSION(nzm) :: tke, tkeshear, tkebuoy, tkeadv, tkepress, tkediss, tkediff
		REAL, DIMENSION(nzm) :: tkesb, tkesbshear, tkesbbuoy, tkesbadv, tkesbdiss
		REAL, DIMENSION(nzm) :: tkm, tkh
		! Radiation
		REAL, DIMENSION(nzm) :: radlwup, radlwdn, radswup, radswdn, radqrlw, radqrsw, radqr
		REAL, DIMENSION(nzm) :: radqrs, radqrc
		! WRF-LES module
		REAL, DIMENSION(nzm) :: ulsvadv, vlsvadv, tplsvadv, qvlsvadv
		! Conditional average
		REAL, DIMENSION(nzm) :: cld, rhocl, ucl, vcl, wcl
		REAL, DIMENSION(nzm) :: dsecl, msecl, liwsecl
		REAL, DIMENSION(nzm) :: tpvcl, tpvcla, tabscl
		REAL, DIMENSION(nzm) :: qtcl, qccl, qicl, qncl, qpcl
		REAL, DIMENSION(nzm) :: etacl, ssrcl
		REAL, DIMENSION(nzm) :: qnccl, qnicl, qnrcl, qnscl, qngcl, qnncl
		REAL, DIMENSION(nzm) :: wacl, tvwcl, twcl, qwcl, qcwcl, qiwcl
		
		! 1D array for high resolution mean profiles defined at ZSH
		REAL, DIMENSION(nzsh) :: u0, v0, rho0, pres0, liwse0, mse0, vdse0, tpv0
		REAL, DIMENSION(nzsh) :: qt0, qc0, qi0, qr0, qp0
		
		! 3D arrays
		! Height
		REAL, DIMENSION(0:nx,nzm,0:ny) :: zs_eta   ! S
		REAL, DIMENSION(1:nx,nzm,1:ny) :: zu_eta   ! U
		REAL, DIMENSION(1:nx,nzm,1:ny) :: zv_eta   ! V
		REAL, DIMENSION(1:nx,nzm,1:ny) :: zuw_eta  ! UW
		REAL, DIMENSION(1:nx,nzm,1:ny) :: zvw_eta  ! VW
		! Indices for interpolation
		INTEGER, DIMENSION(0:nx,nzm,0:ny) :: km1s, k00s, kp1s, kp2s ! S => ZS
		INTEGER, DIMENSION(1:nx,nzm,1:ny) :: km1w, k00w, kp1w, kp2w ! W => ZS
		INTEGER, DIMENSION(1:nx,nzm,1:ny) :: km1u, k00u, kp1u, kp2u ! U => ZS
		INTEGER, DIMENSION(1:nx,nzm,1:ny) :: km1v, k00v, kp1v, kp2v ! V => ZS
		INTEGER, DIMENSION(1:nx,nzm,1:ny) :: km1uw, k00uw, kp1uw, kp2uw ! UW => ZS
		INTEGER, DIMENSION(1:nx,nzm,1:ny) :: km1vw, k00vw, kp1vw, kp2vw ! VW => ZS
		INTEGER, DIMENSION(0:nx,nzm,0:ny) :: km1s2w, k00s2w, kp1s2w, kp2s2w ! S => W
		INTEGER, DIMENSION(1:nx,nzm,1:ny) :: km1s02s, k00s02s, kp1s02s, kp2s02s ! ZSH => S
		INTEGER, DIMENSION(1:nx,nzsh,1:ny) :: km1s2s0, k00s2s0, kp1s2s0, kp2s2s0 ! S => ZSH
		! Diagnostic variables
		! - Variables used for flux calculation: unit conversion, SGS flux diagnosis
		REAL, DIMENSION(0:nx,nzm,0:ny) :: rhow_eta    ! dry air density at W level ==> moist air
		REAL, DIMENSION(1:nx,nzm,1:ny) :: exnerw_eta  ! exner function at W level
		REAL, DIMENSION(1:nx,nzm,1:ny) :: tkhw_eta    ! tkh at W level
		! - Variable used for higher-order moments and conditional average
		REAL, DIMENSION(1:nx,nzm,1:ny) :: dse_eta   ! dry static energy
		REAL, DIMENSION(1:nx,nzm,1:ny) :: liwse_eta ! liquid-ice water static energy
		REAL, DIMENSION(1:nx,nzm,1:ny) :: mse_eta   ! moist static energy
		REAL, DIMENSION(1:nx,nzm,1:ny) :: vdse_eta  ! virtual dry static energy
		REAL, DIMENSION(1:nx,nzm,1:ny) :: tpv_eta   ! virtual potential temperature
		REAL, DIMENSION(1:nx,nzm,1:ny) :: tabs_eta  ! temperature
		REAL, DIMENSION(1:nx,nzm,1:ny) :: qt_eta    ! qv+qc+qi
		REAL, DIMENSION(1:nx,nzm,1:ny) :: qp_eta    ! qr+qs+qg
		REAL, DIMENSION(1:nx,nzm,1:ny) :: qn_eta    ! qc+qi
		REAL, DIMENSION(1:nx,nzm,1:ny) :: qcl_eta   ! qc+qr
		REAL, DIMENSION(1:nx,nzm,1:ny) :: qci_eta   ! qi+qs+qg
		! Mean value interpolated from ZSH
		! - Used for higher-order moments, TKE tendency terms, conditional average
		REAL, DIMENSION(1:nx,nzm,1:ny) :: u0_eta, v0_eta
		REAL, DIMENSION(1:nx,nzm,1:ny) :: pres0_eta
		REAL, DIMENSION(1:nx,nzm,1:ny) :: tpv0_eta
		! Interpolatd from eta level to ZS
		REAL, DIMENSION(0:nx,nzm,0:ny) :: mweight   ! = rho3d/mean(rho), mass weight
		REAL, DIMENSION(0:nx,nzm,0:ny) :: rho3d
		REAL, DIMENSION(1:nx,nzm,1:ny) :: u3d, v3d, w3d
		REAL, DIMENSION(1:nx,nzm,1:ny) :: dse3d, liwse3d, mse3d, vdse3d
		REAL, DIMENSION(1:nx,nzm,1:ny) :: tpv3d, tabs3d
		REAL, DIMENSION(1:nx,nzm,1:ny) :: qt3d, qc3d, qi3d, qn3d, qp3d
		REAL, DIMENSION(1:nx,nzm,1:ny) :: eta3d, ssr3d
		REAL, DIMENSION(1:nx,nzm,1:ny) :: qnc3d, qni3d, qnr3d, qns3d, qng3d, qnn3d
		REAL, DIMENSION(1:nx,nzm,1:ny) :: radqr3d
		! Mass weight for ZSH level
		REAL, DIMENSION(0:nx,nzsh,0:ny) :: mweight0
		
		! 2D arrays
		! Diagnostic variables
		REAL, DIMENSION(nx,nzm) :: lwse_eta  ! liquid water static energy
		REAL, DIMENSION(nx,nzm) :: smse_eta  ! saturation mse
		REAL, DIMENSION(nx,nzm) :: tpl_eta   ! liquid water potential temperature
		REAL, DIMENSION(nx,nzm) :: tpli_eta  ! Liquid-ice water potential temperature
		REAL, DIMENSION(nx,nzm) :: tpe_eta   ! Equivalent potential temperature
		REAL, DIMENSION(nx,nzm) :: qct_eta   ! qc+qr+qi+qs+qg
		REAL, DIMENSION(nx,nzm) :: qsat_eta  ! saturation mixing ratio
		REAL, DIMENSION(nx,nzm) :: relh_eta  ! ralative humidity
		REAL, DIMENSION(nx,nzm) :: qnt_eta   ! total number concentration (currently qnn+qnc)
		REAL, DIMENSION(nx,nzm) :: ncond_rate_eta   ! net condensation rate
		REAL, DIMENSION(nx,nzm) :: ndeps_rate_eta   ! net deposition rate
		REAL, DIMENSION(nx,nzm) :: tcond_rate_eta   ! total conensation rate (ncond+ndeps)
		REAL, DIMENSION(nx,nzm) :: uw_eta, vw_eta
		REAL, DIMENSION(nx,nzm) :: lwsef_eta, lwsefr_eta, lwsefs_eta
		REAL, DIMENSION(nx,nzm) :: liwsef_eta, liwsefr_eta, liwsefs_eta
		REAL, DIMENSION(nx,nzm) :: tplf_eta, tplfr_eta, tplfs_eta
		REAL, DIMENSION(nx,nzm) :: tpvf_eta, tpvfr_eta, tpvfs_eta
		REAL, DIMENSION(nx,nzm) :: qtf_eta, qtfr_eta, qtfs_eta
		REAL, DIMENSION(nx,nzm) :: qpf_eta, qpfr_eta, qpfs_eta
		REAL, DIMENSION(nx,nzm) :: qcf_eta
		REAL, DIMENSION(nx,nzm) :: qif_eta
		REAL, DIMENSION(nx,nzm) :: qrf_eta
		REAL, DIMENSION(nx,nzm) :: qsf_eta
		REAL, DIMENSION(nx,nzm) :: qgf_eta
		REAL, DIMENSION(nx,nzm) :: precip_eta
		REAL, DIMENSION(nx,nzm) :: u2_eta, v2_eta
		REAL, DIMENSION(nx,nzm) :: w2_eta, w3_eta, w4_eta
		REAL, DIMENSION(nx,nzm) :: liwse2_eta, liwse3_eta
		REAL, DIMENSION(nx,nzm) :: mse2_eta, mse3_eta
		REAL, DIMENSION(nx,nzm) :: vdse2_eta, vdse3_eta
		REAL, DIMENSION(nx,nzm) :: qt2_eta, qt3_eta
		REAL, DIMENSION(nx,nzm) :: qc2_eta
		REAL, DIMENSION(nx,nzm) :: qi2_eta
		REAL, DIMENSION(nx,nzm) :: qcl2_eta
		REAL, DIMENSION(nx,nzm) :: qci2_eta
		REAL, DIMENSION(nx,nzm) :: qm2_eta
		REAL, DIMENSION(nx,nzm) :: slqt_eta
		REAL, DIMENSION(nx,nzm) :: shear_uweta, shear_vweta
		REAL, DIMENSION(nx,nzm) :: adv_ueta, adv_veta, adv_seta, adv_uweta, adv_vweta
		REAL, DIMENSION(nx,nzm) :: dif_ueta, dif_veta, dif_weta
		REAL, DIMENSION(nx,nzm) :: prs_weta, prs_ueta, prs_veta
		REAL, DIMENSION(nx,nzm) :: rhop_eta
		
		! Other miscellaneous variables
		INTEGER :: i, j, k, iwrf, jwrf   ! loop index
		INTEGER :: khunt, khunt1, khunt2, khunt3, khunt4, khunt5   ! used for find_indices
		INTEGER :: km1, k00, kp1, kp2   ! interpolation indices
		REAL :: exner, iexner, presw, tmp, f0, fprime, num, ppm1, pp00
	!!	REAL :: fm1, f00, fp1, fp2   ! NOT USED
		REAL :: factor_xy
		REAL(8) :: factor_n
		! MPI buffer 
		INTEGER, PARAMETER :: nbufzsh = 13   ! number of high resolution mean profile
		INTEGER, PARAMETER :: nbufzs  = 1    ! number of mean profile at zs (RHO)
		INTEGER, PARAMETER :: nbufs   = 2    ! number of scalar (LWP, IWP)
		INTEGER :: nbuffer   ! = nzsh*nbufzsh+nzm*nbufzs+nbufs
		REAL(8), DIMENSION(nzsh*nbufzsh+nzm*nbufzs+nbufs) :: buffer1, buffer2
		! Conditional sample
		INTEGER :: ncond
		LOGICAL :: condition, condition_cl
		REAL :: coef
		! For 1D time series
		REAL :: lwp0, iwp0
		INTEGER, DIMENSION(nx,ny) :: topind
		REAL, DIMENSION(nx,ny) :: cwp, cwpl, cwpm, cwph, lwp, iwp
		REAL :: cwpmax
		
		nbuffer = nzsh*nbufzsh + nzm*nbufzs + nbufs
		factor_xy = 1. / REAL(nx*ny)
		factor_n = 1. / DBLE(nsubdomains)
		
		
		!========== Local grid height for S, U, V, UW, VW levels ==========
		! Heihgt conversion is in the subroutine phy_prep in /dyn_em/module_big_step_utilities_em.F
		! S
		DO j = 0, ny   ! 0 for V
			DO k = 1, nzm
				DO i = 0, nx   ! 0 for U
					zs_eta(i,k,j) = 0.5 * ( zw_eta(i,k,j) + zw_eta(i,k+1,j) )
				ENDDO
			ENDDO
		ENDDO
		! U, V, UW, VW
		DO j = 1, ny 
			DO k = 1, nzm
				DO i = 1, nx
					! U
					zu_eta(i,k,j) = 0.5 * ( zs_eta(i-1,k,j) + zs_eta(i,k,j) )
					! V
					zv_eta(i,k,j) = 0.5 * ( zs_eta(i,k,j-1) + zs_eta(i,k,j) )
					! UW
					zuw_eta(i,k,j) = 0.5 * ( zw_eta(i-1,k,j) + zw_eta(i,k,j) )
					! VW
					zvw_eta(i,k,j) = 0.5 * ( zw_eta(i,k,j-1) + zw_eta(i,k,j) )
				ENDDO
			ENDDO
		ENDDO
		
		
		!========== Prepare for interpolation ==========
		! Interpolation information from all grid points to zs is stored.
		! The k loop has to be the most inner loop for find_indices.
		DO j = 0, ny   ! 0 for V
			DO i = 0, nx   ! 0 for U
				! eta ==> zs
				khunt = 0
				DO k = 1, nzm
					! S
					CALL find_indices( zs(k), zs_eta(i,1:nzm,j), nzm, &
						km1s(i,k,j), k00s(i,k,j), kp1s(i,k,j), kp2s(i,k,j), khunt )
				ENDDO
			ENDDO
		ENDDO
		DO j = 1, ny
			DO i = 1, nx
				! eta ==> zs
				khunt1 = 0
				khunt2 = 0
				khunt3 = 0
				khunt4 = 0
				khunt5 = 0
				DO k = 1, nzm
					! W
					CALL find_indices( zs(k), zw_eta(i,1:nzm,j), nzm, &
							km1w(i,k,j), k00w(i,k,j), kp1w(i,k,j), kp2w(i,k,j), khunt1 )
					! U
					CALL find_indices( zs(k), zu_eta(i,1:nzm,j), nzm, &
							km1u(i,k,j), k00u(i,k,j), kp1u(i,k,j), kp2u(i,k,j), khunt2 )
					! V
					CALL find_indices( zs(k), zv_eta(i,1:nzm,j), nzm, &
							km1v(i,k,j), k00v(i,k,j), kp1v(i,k,j), kp2v(i,k,j), khunt3 )
					! UW
					CALL find_indices( zs(k), zuw_eta(i,1:nzm,j), nzm, &
							km1uw(i,k,j), k00uw(i,k,j), kp1uw(i,k,j), kp2uw(i,k,j), khunt4 )
					! VW
					CALL find_indices( zs(k), zvw_eta(i,1:nzm,j), nzm, &
							km1vw(i,k,j), k00vw(i,k,j), kp1vw(i,k,j), kp2vw(i,k,j), khunt5 )
				ENDDO
			ENDDO
		ENDDO
		
		
		!========== Diagnose variables at eta level / unit conversion ==========
		DO j = 1, ny
			DO k = 1, nzm
				DO i = 1, nx
					! QT, QN, QP, QCL, QCI
					qt_eta(i,k,j) = qv_eta(i,k,j) + qc_eta(i,k,j) + qi_eta(i,k,j)
					qn_eta(i,k,j) = qc_eta(i,k,j) + qi_eta(i,k,j)
					qp_eta(i,k,j) = qr_eta(i,k,j) + qs_eta(i,k,j) + qg_eta(i,k,j)
					qcl_eta(i,k,j) = qc_eta(i,k,j) + qr_eta(i,k,j)
					qci_eta(i,k,j) = qi_eta(i,k,j) + qs_eta(i,k,j) + qg_eta(i,k,j)
					! TPV
					tpv_eta(i,k,j) = tp_eta(i,k,j) * (1.+epsv*qv_eta(i,k,j)-qn_eta(i,k,j)-qp_eta(i,k,j))
					! Exner function
					exner = ( ps_eta(i,k,j) / 1000. ) ** (rgas/cp)
					! TABS
					tabs_eta(i,k,j) = tp_eta(i,k,j) * exner
					! DSE
					dse_eta(i,k,j) = tabs_eta(i,k,j) + fac_ggr * zs_eta(i,k,j)
					! LIWSE
					liwse_eta(i,k,j) = dse_eta(i,k,j) - fac_cond*qcl_eta(i,k,j) - fac_sub*qci_eta(i,k,j)
					! MSE
					mse_eta(i,k,j) = dse_eta(i,k,j) + fac_cond * qv_eta(i,k,j)
					! VDSE
					vdse_eta(i,k,j) = tpv_eta(i,k,j) * exner + fac_ggr * zs_eta(i,k,j)
					
					! RADQRLW, RADQRSW: K/s for TP ==> K/s for TABS
					radqrlw_eta(i,k,j) = radqrlw_eta(i,k,j) * exner
					radqrsw_eta(i,k,j) = radqrsw_eta(i,k,j) * exner
				ENDDO
			ENDDO
		ENDDO
		! RHOW, EXNERW, TKHW, coordinate transformation for vertical flux
		DO j = 0, ny   ! VW
			DO i = 0, nx   ! UW
				khunt = 0
				DO k = 1, nzm   ! The k loop has to be the most inner loop for find_indices
					! S ==> W
					CALL find_indices( zw_eta(i,k,j), zs_eta(i,1:nzm,j), nzm, &
							km1s2w(i,k,j), k00s2w(i,k,j), kp1s2w(i,k,j), kp2s2w(i,k,j), khunt )
				ENDDO
			ENDDO
		ENDDO
		! RHOW
		DO j = 0, ny   ! VW
			DO k = 1, nzm
				DO i = 0, nx   ! UW
					km1 = km1s2w(i,k,j)
					k00 = k00s2w(i,k,j)
					kp1 = kp1s2w(i,k,j)
					kp2 = kp2s2w(i,k,j)
					rhow_eta(i,k,j) =  mono_cubic_interp( zw_eta(i,k,j), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j),   zs_eta(i,k00,j),   zs_eta(i,kp1,j),   zs_eta(i,kp2,j), &
							rhod_eta(i,km1,j), rhod_eta(i,k00,j), rhod_eta(i,kp1,j), rhod_eta(i,kp2,j) )
				ENDDO
			ENDDO
		ENDDO
		! EXNERW, TKHW
		DO j = 1, ny
			DO k = 1, nzm
				DO i = 1, nx
					km1 = km1s2w(i,k,j)
					k00 = k00s2w(i,k,j)
					kp1 = kp1s2w(i,k,j)
					kp2 = kp2s2w(i,k,j)
					! EXNERW
					presw = mono_cubic_interp( zw_eta(i,k,j), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							ps_eta(i,km1,j), ps_eta(i,k00,j), ps_eta(i,kp1,j), ps_eta(i,kp2,j) )
					exnerw_eta(i,k,j) = ( presw / 1000. ) ** (rgas/cp)
					! TKHW
					tkhw_eta(i,k,j) = mono_cubic_interp( zw_eta(i,k,j), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j),  zs_eta(i,k00,j),  zs_eta(i,kp1,j),  zs_eta(i,kp2,j), &
							tkh_eta(i,km1,j), tkh_eta(i,k00,j), tkh_eta(i,kp1,j), tkh_eta(i,kp2,j) )
				ENDDO
			ENDDO
		ENDDO
		! Transformation from omega to w for vertical flux [unit m/s]
		DO j = 1, ny
			jwrf = j0 + j
			DO k = 1, nzm
				DO i = 1, nx
					iwrf = i0 + i
					wadv_eta(iwrf,k,jwrf) = - wadv_eta(iwrf,k,jwrf) / rhow_eta(i,k,j) / ggr
					uwle_eta(iwrf,k,jwrf) = - uwle_eta(iwrf,k,jwrf) &
					                      / ( 0.5 * ( rhow_eta(i-1,k,j) + rhow_eta(i,k,j) ) * ggr )
					vwle_eta(iwrf,k,jwrf) = - vwle_eta(iwrf,k,jwrf) &
					                      / ( 0.5 * ( rhow_eta(i,k,j-1) + rhow_eta(i,k,j) ) * ggr )
					twle_eta(iwrf,k,jwrf) = - twle_eta(iwrf,k,jwrf) / rhow_eta(i,k,j) / ggr
					mwle_eta(iwrf,k,jwrf,:) = - mwle_eta(iwrf,k,jwrf,:) / rhow_eta(i,k,j) / ggr
					swle_eta(iwrf,k,jwrf,:) = - swle_eta(iwrf,k,jwrf,:) / rhow_eta(i,k,j) / ggr
					ewle_eta(iwrf,k,jwrf) = - ewle_eta(iwrf,k,jwrf) / rhow_eta(i,k,j) / ggr
				ENDDO
			ENDDO
		ENDDO
		! RHOW: Convert from dry air to moist air
		DO j = 0, ny   ! VW
			DO k = 1, nzm
				DO i = 0, nx   ! UW
					! S ==> W
					km1 = km1s2w(i,k,j)
					k00 = k00s2w(i,k,j)
					kp1 = kp1s2w(i,k,j)
					kp2 = kp2s2w(i,k,j)
					rhow_eta(i,k,j) = mono_cubic_interp( zw_eta(i,k,j), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j),   zs_eta(i,k00,j),   zs_eta(i,kp1,j),   zs_eta(i,kp2,j), &
							rhom_eta(i,km1,j), rhom_eta(i,k00,j), rhom_eta(i,kp1,j), rhom_eta(i,kp2,j) )
				ENDDO
			ENDDO
		ENDDO
		
		
		!========== Inerpolate from eta to ZS for futur use ==========
		DO j = 0, ny   ! 0 for V
			DO k = 1, nzm 
				DO i = 0, nx  ! 0 for U
					! RHO
					km1 = km1s(i,k,j)
					k00 = k00s(i,k,j)
					kp1 = kp1s(i,k,j)
					kp2 = kp2s(i,k,j)
					rho3d(i,k,j) = mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j),   zs_eta(i,k00,j),   zs_eta(i,kp1,j),   zs_eta(i,kp2,j), &
							rhom_eta(i,km1,j), rhom_eta(i,k00,j), rhom_eta(i,kp1,j), rhom_eta(i,kp2,j) )
				ENDDO
			ENDDO
		ENDDO
		DO j = 1, ny
			DO k = 1, nzm
				DO i = 1, nx
					! U
					km1 = km1u(i,k,j)
					k00 = k00u(i,k,j)
					kp1 = kp1u(i,k,j)
					kp2 = kp2u(i,k,j)
					u3d(i,k,j) = mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zu_eta(i,km1,j), zu_eta(i,k00,j), zu_eta(i,kp1,j), zu_eta(i,kp2,j), &
							u_eta(i,km1,j),  u_eta(i,k00,j),  u_eta(i,kp1,j),  u_eta(i,kp2,j) )
					! V
					km1 = km1v(i,k,j)
					k00 = k00v(i,k,j)
					kp1 = kp1v(i,k,j)
					kp2 = kp2v(i,k,j)
					v3d(i,k,j) = mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zv_eta(i,km1,j), zv_eta(i,k00,j), zv_eta(i,kp1,j), zv_eta(i,kp2,j), &
							v_eta(i,km1,j),  v_eta(i,k00,j),  v_eta(i,kp1,j),  v_eta(i,kp2,j) )
					! W
					km1 = km1w(i,k,j)
					k00 = k00w(i,k,j)
					kp1 = kp1w(i,k,j)
					kp2 = kp2w(i,k,j)
					w3d(i,k,j) = mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zw_eta(i,km1,j), zw_eta(i,k00,j), zw_eta(i,kp1,j), zw_eta(i,kp2,j), &
							w_eta(i,km1,j),  w_eta(i,k00,j),  w_eta(i,kp1,j),  w_eta(i,kp2,j) )
					! SCALR
					km1 = km1s(i,k,j)
					k00 = k00s(i,k,j)
					kp1 = kp1s(i,k,j)
					kp2 = kp2s(i,k,j)
					! DSE
					dse3d(i,k,j) = mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j),  zs_eta(i,k00,j),  zs_eta(i,kp1,j),  zs_eta(i,kp2,j), &
							dse_eta(i,km1,j), dse_eta(i,k00,j), dse_eta(i,kp1,j), dse_eta(i,kp2,j) )
					! LIWSE
					liwse3d(i,k,j) = mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j),   zs_eta(i,k00,j),   zs_eta(i,kp1,j),   zs_eta(i,kp2,j), &
							liwse_eta(i,km1,j),liwse_eta(i,k00,j),liwse_eta(i,kp1,j),liwse_eta(i,kp2,j) )
					! MSE
					mse3d(i,k,j) = mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j),  zs_eta(i,k00,j),  zs_eta(i,kp1,j),  zs_eta(i,kp2,j), &
							mse_eta(i,km1,j), mse_eta(i,k00,j), mse_eta(i,kp1,j), mse_eta(i,kp2,j) )
					! VDSE
					vdse3d(i,k,j) = mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j),   zs_eta(i,k00,j),   zs_eta(i,kp1,j),   zs_eta(i,kp2,j), &
							vdse_eta(i,km1,j), vdse_eta(i,k00,j), vdse_eta(i,kp1,j), vdse_eta(i,kp2,j) )
					! TPV
					tpv3d(i,k,j) = mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j),  zs_eta(i,k00,j),  zs_eta(i,kp1,j),  zs_eta(i,kp2,j), &
							tpv_eta(i,km1,j), tpv_eta(i,k00,j), tpv_eta(i,kp1,j), tpv_eta(i,kp2,j) )
					! TABS
					tabs3d(i,k,j) = mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j),   zs_eta(i,k00,j),   zs_eta(i,kp1,j),   zs_eta(i,kp2,j), &
							tabs_eta(i,km1,j), tabs_eta(i,k00,j), tabs_eta(i,kp1,j), tabs_eta(i,kp2,j) )
					! QT
					qt3d(i,k,j) = mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							qt_eta(i,km1,j), qt_eta(i,k00,j), qt_eta(i,kp1,j), qt_eta(i,kp2,j) )
					! QC
					qc3d(i,k,j) = mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							qc_eta(i,km1,j), qc_eta(i,k00,j), qc_eta(i,kp1,j), qc_eta(i,kp2,j) )
					! QI
					qi3d(i,k,j) = mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							qi_eta(i,km1,j), qi_eta(i,k00,j), qi_eta(i,kp1,j), qi_eta(i,kp2,j) )
					! QN
					qn3d(i,k,j) = mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							qn_eta(i,km1,j), qn_eta(i,k00,j), qn_eta(i,kp1,j), qn_eta(i,kp2,j) )
					! QP
					qp3d(i,k,j) = mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							qp_eta(i,km1,j), qp_eta(i,k00,j), qp_eta(i,kp1,j), qp_eta(i,kp2,j) )
					! ETA
					eta3d(i,k,j) = mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j),  zs_eta(i,k00,j),  zs_eta(i,kp1,j),  zs_eta(i,kp2,j), &
							eta_eta(i,km1,j), eta_eta(i,k00,j), eta_eta(i,kp1,j), eta_eta(i,kp2,j) )
					! SSR
					ssr3d(i,k,j) = mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j),  zs_eta(i,k00,j),  zs_eta(i,kp1,j),  zs_eta(i,kp2,j), &
							ssr_eta(i,km1,j), ssr_eta(i,k00,j), ssr_eta(i,kp1,j), ssr_eta(i,kp2,j) )
					! QNC
					qnc3d(i,k,j) = mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j),  zs_eta(i,k00,j),  zs_eta(i,kp1,j),  zs_eta(i,kp2,j), &
							qnc_eta(i,km1,j), qnc_eta(i,k00,j), qnc_eta(i,kp1,j), qnc_eta(i,kp2,j) )
					! QNI
					qni3d(i,k,j) = mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j),  zs_eta(i,k00,j),  zs_eta(i,kp1,j),  zs_eta(i,kp2,j), &
							qni_eta(i,km1,j), qni_eta(i,k00,j), qni_eta(i,kp1,j), qni_eta(i,kp2,j) )
					! QNR
					qnr3d(i,k,j) = mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j),  zs_eta(i,k00,j),  zs_eta(i,kp1,j),  zs_eta(i,kp2,j), &
							qnr_eta(i,km1,j), qnr_eta(i,k00,j), qnr_eta(i,kp1,j), qnr_eta(i,kp2,j) )
					! QNS
					qns3d(i,k,j) = mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j),  zs_eta(i,k00,j),  zs_eta(i,kp1,j),  zs_eta(i,kp2,j), &
							qns_eta(i,km1,j), qns_eta(i,k00,j), qns_eta(i,kp1,j), qns_eta(i,kp2,j) )
					! QNG
					qng3d(i,k,j) = mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j),  zs_eta(i,k00,j),  zs_eta(i,kp1,j),  zs_eta(i,kp2,j), &
							qng_eta(i,km1,j), qng_eta(i,k00,j), qng_eta(i,kp1,j), qng_eta(i,kp2,j) )
					! QNN
					qnn3d(i,k,j) = mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j),  zs_eta(i,k00,j),  zs_eta(i,kp1,j),  zs_eta(i,kp2,j), &
							qnn_eta(i,km1,j), qnn_eta(i,k00,j), qnn_eta(i,kp1,j), qnn_eta(i,kp2,j) )
				ENDDO ! i
			ENDDO ! k
		ENDDO ! j
		
		
		!========== RHO ==========
		! The mean profile of density at zs is required for the mass-weighted horizontal mean,
		! which is an appropriate way for the full compressible system.
		! Summation here, averaging and MPI exchange later
		rho(:) = 0.
		DO j = 1, ny
			DO k = 1, nzm
				DO i = 1, nx
					rho(k) = rho(k) + rho3d(i,k,j)
				ENDDO
			ENDDO
		ENDDO
		! Store rho for output
		CALL hbuf_put( 'RHO', rho, factor_xy )   ! kg/m3
		
		
		!========== High resolution mean profile for higher order moments ==========
		! The high resolution mean profiles are a mean profile computed to the zsh levels.
		! Number of zsh level is (nzm-1)*hrlayers+1. Higher order moments, e.g., (u_eta-u0)**2,
		! are computed on the eta level, then interpolated to zs level for summation. This requires
		! interpolation of the mean profile to the eta level. In order to have the mean value on
		! the eta level accurate enough, the mean profile should be constructed with a high
		! resolution. Imagine all horizontal grid points have a same value, so each local value is
		! exactly equal to the mean. If the mean profile is constructed with low resolution, then
		! the resulting mean value at the local eta height, which is interpolated from the
		! constructed mean profile with a coarse resolution, will be largely under/over estimated
		! compared with the true mean value.
		!
		! RHO
		! Store interpolated value to mweight0
		DO j = 0, ny   ! 0 for V
			DO i = 0, nx   ! 0 for U
				khunt = 0
				DO k = 1, nzsh   ! The k loop has to be the most inner loop for find_indices
					CALL find_indices( zsh(k), zs_eta(i,1:nzm,j), nzm, km1, k00, kp1, kp2, khunt )
					mweight0(i,k,j) = mono_cubic_interp( zsh(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j),   zs_eta(i,k00,j),   zs_eta(i,kp1,j),   zs_eta(i,kp2,j), &
							rhom_eta(i,km1,j), rhom_eta(i,k00,j), rhom_eta(i,kp1,j), rhom_eta(i,kp2,j) )
					IF ( i > 0 .AND. j > 0 ) THEN   ! Not good way but...
						km1s2s0(i,k,j) = km1
						k00s2s0(i,k,j) = k00
						kp1s2s0(i,k,j) = kp1
						kp2s2s0(i,k,j) = kp2
					ENDIF
				ENDDO
			ENDDO
		ENDDO
		! Summation (averaging and MPI exchange later)
		! Interpolated quantity is multiplied by mweight0 (density).
		u0(:) = 0.
		v0(:) = 0.
		rho0(:) = 0.
		pres0(:)  = 0.
		liwse0(:) = 0.
		mse0(:)  = 0.
		vdse0(:) = 0.
		tpv0(:)  = 0.
		qt0(:) = 0.
		qc0(:) = 0.
		qi0(:) = 0.
		qr0(:) = 0.
		qp0(:) = 0.
		DO j = 1, ny
			DO i = 1, nx
				khunt1 = 0
				khunt2 = 0
				DO k = 1, nzsh   ! The k loop has to be the most inner loop for find_indices
					! U
					CALL find_indices( zsh(k), zu_eta(i,1:nzm,j), nzm, km1, k00, kp1, kp2, khunt1 )
					u0(k) = u0(k) + mono_cubic_interp( zsh(k), km1, k00, kp1, kp2, &
							zu_eta(i,km1,j), zu_eta(i,k00,j), zu_eta(i,kp1,j), zu_eta(i,kp2,j), &
							u_eta(i,km1,j),  u_eta(i,k00,j),  u_eta(i,kp1,j),  u_eta(i,kp2,j) ) &
							* 0.5 * ( mweight0(i-1,k,j) + mweight0(i,k,j) )
					! V
					CALL find_indices( zsh(k), zv_eta(i,1:nzm,j), nzm, km1, k00, kp1, kp2, khunt2 )
					v0(k) = v0(k) + mono_cubic_interp( zsh(k), km1, k00, kp1, kp2, &
							zv_eta(i,km1,j), zv_eta(i,k00,j), zv_eta(i,kp1,j), zv_eta(i,kp2,j), &
							v_eta(i,km1,j),  v_eta(i,k00,j),  v_eta(i,kp1,j),  v_eta(i,kp2,j) ) &
							* 0.5 * ( mweight0(i,k,j-1) + mweight0(i,k,j) )
				ENDDO
			ENDDO
		ENDDO
		DO j = 1, ny
			DO k = 1, nzsh
				DO i = 1, nx
					! S
					km1 = km1s2s0(i,k,j)
					k00 = k00s2s0(i,k,j)
					kp1 = kp1s2s0(i,k,j)
					kp2 = kp2s2s0(i,k,j)
					! RHO
					rho0(k) = rho0(k) + mweight0(i,k,j)
					! PRES: no need to multiply by density
					pres0(k) = pres0(k) + mono_cubic_interp( zsh(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							ps_eta(i,km1,j), ps_eta(i,k00,j), ps_eta(i,kp1,j), ps_eta(i,kp2,j) )
					! LIWSE
					liwse0(k) = liwse0(k) + mono_cubic_interp( zsh(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j),   zs_eta(i,k00,j),   zs_eta(i,kp1,j),   zs_eta(i,kp2,j), &
							liwse_eta(i,km1,j),liwse_eta(i,k00,j),liwse_eta(i,kp1,j),liwse_eta(i,kp2,j) ) &
							* mweight0(i,k,j)
					! VDSE
					vdse0(k) = vdse0(k) + mono_cubic_interp( zsh(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j),   zs_eta(i,k00,j),   zs_eta(i,kp1,j),   zs_eta(i,kp2,j), &
							vdse_eta(i,km1,j), vdse_eta(i,k00,j), vdse_eta(i,kp1,j), vdse_eta(i,kp2,j) ) &
							* mweight0(i,k,j)
					! MSE
					mse0(k) = mse0(k) + mono_cubic_interp( zsh(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j),  zs_eta(i,k00,j),  zs_eta(i,kp1,j),  zs_eta(i,kp2,j), &
							mse_eta(i,km1,j), mse_eta(i,k00,j), mse_eta(i,kp1,j), mse_eta(i,kp2,j) ) &
							* mweight0(i,k,j)
					! TPV
					tpv0(k) = tpv0(k) + mono_cubic_interp( zsh(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j),  zs_eta(i,k00,j),  zs_eta(i,kp1,j),  zs_eta(i,kp2,j), &
							tpv_eta(i,km1,j), tpv_eta(i,k00,j), tpv_eta(i,kp1,j), tpv_eta(i,kp2,j) ) &
							* mweight0(i,k,j)
					! QT
					qt0(k) = qt0(k) + mono_cubic_interp( zsh(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							qt_eta(i,km1,j), qt_eta(i,k00,j), qt_eta(i,kp1,j), qt_eta(i,kp2,j) ) &
							* mweight0(i,k,j)
					! QC
					qc0(k) = qc0(k) + mono_cubic_interp( zsh(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							qc_eta(i,km1,j), qc_eta(i,k00,j), qc_eta(i,kp1,j), qc_eta(i,kp2,j) ) &
							* mweight0(i,k,j)
					! QI
					qi0(k) = qi0(k) + mono_cubic_interp( zsh(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							qi_eta(i,km1,j), qi_eta(i,k00,j), qi_eta(i,kp1,j), qi_eta(i,kp2,j) ) &
							* mweight0(i,k,j)
					! QR
					qr0(k) = qr0(k) + mono_cubic_interp( zsh(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							qr_eta(i,km1,j), qr_eta(i,k00,j), qr_eta(i,kp1,j), qr_eta(i,kp2,j) ) &
							* mweight0(i,k,j)
					! QP
					qp0(k) = qp0(k) + mono_cubic_interp( zsh(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							qp_eta(i,km1,j), qp_eta(i,k00,j), qp_eta(i,kp1,j), qp_eta(i,kp2,j) ) &
							* mweight0(i,k,j)
				ENDDO ! k
			ENDDO ! i
		ENDDO ! j
		
		
		!========== LWP2, IWP2, PW, CWP, IWP, RWP, SWP, GWP, cloud fractions ==========
		! The subdomain-sum of liquid wter path & ice water path are computed for the variance output.
		! Also, other column waters and cloud fractions are prepard for 1D output.
		cwp(:,:) = 0.   ! cloud water+ice path doe cloud fraction diagnosis
		lwp(:,:) = 0.   ! liquid water path = qc + qr
		iwp(:,:) = 0.   ! ice water path = qi + qs + qg
		cwpl(:,:) = 0.
		cwpm(:,:) = 0.
		cwph(:,:) = 0.
		topind(:,:) = 1
		IF (CSRM) THEN
			cwpmax = 0.02
		ELSE
			cwpmax = 0.0
		ENDIF
		DO j = 1, ny
			DO k = nzm, 1, -1
				DO i = 1, nx
					tmp = rhom_eta(i,k,j) * ( zw_eta(i,k+1,j) - zw_eta(i,k,j) )
					! Shaded cloud fraction
					cwp(i,j) = cwp(i,j) + tmp * qn_eta(i,k,j)   ! kg/m2
					IF (ps_eta(i,k,j) >= 700.) THEN
						cwpl(i,j) = cwpl(i,j) + tmp * qn_eta(i,k,j)
					ELSE IF (ps_eta(i,k,j) <= 400.) THEN
						cwph(i,j) = cwph(i,j) + tmp * qn_eta(i,k,j)
					ELSE
						cwpm(i,j) = cwpm(i,j) + tmp * qn_eta(i,k,j)
					ENDIF
					IF(cwp(i,j) > cwpmax .AND. topind(i,j) == 1) topind(i,j) = k
					! 1D output
					s_pw  = s_pw  + tmp * qv_eta(i,k,j)   ! mm <== (kg/m2)/(1000kg/m3 rho of H2O)=m/1000
					s_cwp = s_cwp + tmp * qc_eta(i,k,j) * 1.E3   ! g/m2 <==(kg/m2)1000
					s_iwp = s_iwp + tmp * qi_eta(i,k,j) * 1.E3   ! g/m2
					s_rwp = s_rwp + tmp * qr_eta(i,k,j) * 1.E3   ! g/m2
					s_swp = s_swp + tmp * qs_eta(i,k,j) * 1.E3   ! g/m2
					s_gwp = s_gwp + tmp * qg_eta(i,k,j) * 1.E3   ! g/m2
					! Variance of LWP and IWP
					lwp(i,j) = lwp(i,j) + tmp * qcl_eta(i,k,j)   ! kg/m2
					iwp(i,j) = iwp(i,j) + tmp * qci_eta(i,k,j)   ! kg/m2
				ENDDO
			ENDDO
		ENDDO
		! Shaded cloud fraction
		DO j = 1, ny
			DO i = 1, nx
				IF (cwp(i,j)  > cwpmax) s_acld  = s_acld  + 1.
				IF (cwpl(i,j) > cwpmax) s_acldl = s_acldl + 1.
				IF (cwpm(i,j) > cwpmax) s_acldm = s_acldm + 1.
				IF (cwph(i,j) > cwpmax) s_acldh = s_acldh + 1.
				IF (tabs_eta(i,topind(i,j),j) < 245.) s_acldcold = s_acldcold + 1.
			ENDDO
		ENDDO
		! Summation of LWP and IWP (averaging and MPI exchange next)
		lwp0 = 0.   ! mean liquid water path
		iwp0 = 0.   ! mean ice water path
		DO j = 1, ny
			DO i = 1, nx
				lwp0 = lwp0 + lwp(i,j)
				iwp0 = iwp0 + iwp(i,j)
			ENDDO
		ENDDO
		
		
		!========== Averaging RHO, high resolution mean profile, LWP, IWP ==========
		! ZSH
		u0(:)     = u0(:) * factor_xy
		v0(:)     = v0(:) * factor_xy
		rho0(:)   = rho0(:) * factor_xy
		pres0(:)  = pres0(:) * factor_xy
		liwse0(:) = liwse0(:) * factor_xy
		mse0(:)   = mse0(:) * factor_xy
		vdse0(:)  = vdse0(:) * factor_xy
		tpv0(:)   = tpv0(:) * factor_xy
		qt0(:)    = qt0(:) * factor_xy
		qc0(:)    = qc0(:) * factor_xy
		qi0(:)    = qi0(:) * factor_xy
		qr0(:)    = qr0(:) * factor_xy
		qp0(:)    = qp0(:) * factor_xy
		! ZS
		rho(:)  = rho(:) * factor_xy
		! 1 element
		lwp0 = lwp0 * factor_xy
		iwp0 = iwp0 * factor_xy
		IF (dompi) THEN
			buffer1(        1:   nzsh) = u0(:)
			buffer1(   nzsh+1: 2*nzsh) = v0(:)
			buffer1( 2*nzsh+1: 3*nzsh) = rho0(:)
			buffer1( 3*nzsh+1: 4*nzsh) = pres0(:)
			buffer1( 4*nzsh+1: 5*nzsh) = liwse0(:)
			buffer1( 5*nzsh+1: 6*nzsh) = mse0(:)
			buffer1( 6*nzsh+1: 7*nzsh) = vdse0(:)
			buffer1( 7*nzsh+1: 8*nzsh) = tpv0(:)
			buffer1( 8*nzsh+1: 9*nzsh) = qt0(:)
			buffer1( 9*nzsh+1:10*nzsh) = qc0(:)
			buffer1(10*nzsh+1:11*nzsh) = qi0(:)
			buffer1(11*nzsh+1:12*nzsh) = qr0(:)
			buffer1(12*nzsh+1:13*nzsh) = qp0(:)
			buffer1(13*nzsh+1:13*nzsh+nzm) = rho(:)
			buffer1(13*nzsh+nzm+1) = lwp0
			buffer1(13*nzsh+nzm+2) = iwp0
			CALL task_sum_real8( buffer1, buffer2, nbuffer )
			buffer2(:) = buffer2(:) * factor_n
			u0(:)     = buffer2(        1:   nzsh)
			v0(:)     = buffer2(   nzsh+1: 2*nzsh)
			rho0(:)   = buffer2( 2*nzsh+1: 3*nzsh)
			pres0(:)  = buffer2( 3*nzsh+1: 4*nzsh)
			liwse0(:) = buffer2( 4*nzsh+1: 5*nzsh)
			mse0(:)   = buffer2( 5*nzsh+1: 6*nzsh)
			vdse0(:)  = buffer2( 6*nzsh+1: 7*nzsh)
			tpv0(:)   = buffer2( 7*nzsh+1: 8*nzsh)
			qt0(:)    = buffer2( 8*nzsh+1: 9*nzsh)
			qc0(:)    = buffer2( 9*nzsh+1:10*nzsh)
			qi0(:)    = buffer2(10*nzsh+1:11*nzsh)
			qr0(:)    = buffer2(11*nzsh+1:12*nzsh)
			qp0(:)    = buffer2(12*nzsh+1:13*nzsh)
			rho(:)    = buffer2(13*nzsh+1:13*nzsh+nzm)
			lwp0      = buffer2(13*nzsh+nzm+1)
			iwp0      = buffer2(13*nzsh+nzm+2)
		ENDIF ! dompi
		! Mass-weighted average for high resolution mean profile
		u0(:) = u0(:) / rho0(:)
		v0(:) = v0(:) / rho0(:)
		liwse0(:) = liwse0(:) / rho0(:)
		mse0(:)   = mse0(:) / rho0(:)
		vdse0(:)  = vdse0(:) / rho0(:)
		tpv0(:)   = tpv0(:) / rho0(:)
		qt0(:) = qt0(:) / rho0(:)
		qc0(:) = qc0(:) / rho0(:)
		qi0(:) = qi0(:) / rho0(:)
		qr0(:) = qr0(:) / rho0(:)
		qp0(:) = qp0(:) / rho0(:)
		! Mass weight for output
		DO j = 0, ny   ! 0 for V
			DO k = 1, nzm
				DO i = 0, nx   ! 0 for U
					mweight(i,k,j) = rho3d(i,k,j) / rho(k)
				ENDDO
			ENDDO
		ENDDO
		! LWP2 & IWP2 for 1D output
		DO j = 1, ny
			DO i = 1, nx
				s_lwp2 = s_lwp2 + ( ( lwp(i,j) - lwp0 ) ** 2 ) * 1.E6   ! g2/m4
				s_iwp2 = s_iwp2 + ( ( iwp(i,j) - iwp0 ) ** 2 ) * 1.E6   ! g2/m4
			ENDDO
		ENDDO
		
		
		!========== Interpolate high resolution mean profile from ZSH to eta ==========
		DO j = 1, ny
			DO i = 1, nx
				khunt1 = 0
				khunt2 = 0
				khunt3 = 0
				DO k = 1, nzm   ! The k loop has to be the most inner loop for find_indices
					! U
					CALL find_indices( zu_eta(i,k,j), zsh, nzsh, km1, k00, kp1, kp2, khunt1 )
					u0_eta(i,k,j) = mono_cubic_interp( zu_eta(i,k,j), km1, k00, kp1, kp2, &
							zsh(km1), zsh(k00), zsh(kp1), zsh(kp2), u0(km1), u0(k00), u0(kp1), u0(kp2) )
					! V
					CALL find_indices( zv_eta(i,k,j), zsh, nzsh, km1, k00, kp1, kp2, khunt2 )
					v0_eta(i,k,j) = mono_cubic_interp( zv_eta(i,k,j), km1, k00, kp1, kp2, &
							zsh(km1), zsh(k00), zsh(kp1), zsh(kp2), v0(km1), v0(k00), v0(kp1), v0(kp2) )
					! S
					CALL find_indices( zs_eta(i,k,j), zsh, nzsh, km1, k00, kp1, kp2, khunt3 )
					! PRES
					pres0_eta(i,k,j) = mono_cubic_interp( zs_eta(i,k,j), km1, k00, kp1, kp2, &
							zsh(km1),   zsh(k00),   zsh(kp1),   zsh(kp2), &
							pres0(km1), pres0(k00), pres0(kp1), pres0(kp2) )
					! TPV
					tpv0_eta(i,k,j) = mono_cubic_interp( zs_eta(i,k,j), km1, k00, kp1, kp2, &
							zsh(km1),  zsh(k00 ), zsh(kp1),  zsh(kp2), &
							tpv0(km1), tpv0(k00), tpv0(kp1), tpv0(kp2) )
					! Store interpolation indices
					km1s02s(i,k,j) = km1
					k00s02s(i,k,j) = k00
					kp1s02s(i,k,j) = kp1
					kp2s02s(i,k,j) = kp2
				ENDDO ! k
			ENDDO ! i
		ENDDO ! j
		
		
		!========== RHO for conditional average, clear & cloudy sky column ==========
		! This is a good place to prepare RHO for standard conditional average with condavg_mask and
		! other conditional statistics, e.g. clear & cloudy sky column.
!!!!!!!!!! FOR FUTURE VERSION !!!!!!!!!!
		
		
		
		
		
		!========== PRES, U, V, and thermodynamic variables ==========
		pres(:) = 0.
		u(:) = 0.
		v(:) = 0.
		dse(:)   = 0.
		lwse(:)  = 0.
		liwse(:) = 0.
		mse(:)  = 0.
		vdse(:) = 0.
		smse(:) = 0.
		tp(:)   = 0.
		tpl(:)  = 0.
		tpli(:) = 0.
		tpe(:) = 0.
		tpv(:) = 0.
		tabs(:) = 0.
		qt(:) = 0.
		qv(:) = 0.
		qc(:) = 0.
		qi(:) = 0.
		qr(:) = 0.
		qs(:) = 0.
		qg(:) = 0.
		qn(:) = 0.
		qp(:) = 0.
		qcl(:) = 0.
		qci(:) = 0.
		qct(:) = 0.
		qsat(:) = 0.
		relh(:) = 0.
		eta(:) = 0.
		ssr(:) = 0.
		qnc(:) = 0.
		qni(:) = 0.
		qnr(:) = 0.
		qns(:) = 0.
		qng(:) = 0.
		qnn(:) = 0.
		qnt(:) = 0.
		cond_rate(:)  = 0.
		evap_rate(:)  = 0.
		ncond_rate(:) = 0.
		deps_rate(:)  = 0.
		subl_rate(:)  = 0.
		ndeps_rate(:) = 0.
		tcond_rate(:) = 0.
		DO j = 1, ny
			jwrf = j0 + j
			! Diagnose variables and store them into 2D array
			DO k = 1, nzm
				DO i = 1, nx
					iwrf = i0 + i
					! LWSE, SMSE
					lwse_eta(i,k) = dse_eta(i,k,j) - fac_cond * qcl_eta(i,k,j)
					smse_eta(i,k) = dse_eta(i,k,j) + fac_cond * qsatw( tabs_eta(i,k,j), ps_eta(i,k,j) )
					! Inverse of exner funciton
					iexner = ( 1000. / ps_eta(i,k,j) ) ** (rgas/cp)
					! TPL, TPLI, TPE
					tpl_eta(i,k)  = tp_eta(i,k,j) - fac_cond * qcl_eta(i,k,j) * iexner
					tpli_eta(i,k) = tpl_eta(i,k)  - fac_sub * qci_eta(i,k,j) * iexner
					tpe_eta(i,k)  = tp_eta(i,k,j) + fac_cond * qv_eta(i,k,j) * iexner
					! QCT
					qct_eta(i,k)  = qcl_eta(i,k,j) + qci_eta(i,k,j)
					! QSAT, RELH
					qsat_eta(i,k) = qsatw( tabs_eta(i,k,j), ps_eta(i,k,j) )
					relh_eta(i,k) = qv_eta(i,k,j) / qsat_eta(i,k)
					! QNT
					qnt_eta(i,k)  = qnn_eta(i,k,j) + qnc_eta(i,k,j)
					! NCOND, NDEPS, TCOND
					ncond_rate_eta(i,k) = cond_rate_eta(iwrf,k,jwrf)-evap_rate_eta(iwrf,k,jwrf)
					ndeps_rate_eta(i,k) = deps_rate_eta(iwrf,k,jwrf)-subl_rate_eta(iwrf,k,jwrf)
					tcond_rate_eta(i,k) = ncond_rate_eta(i,k) + ndeps_rate_eta(i,k)
				ENDDO
			ENDDO
			
			! Interpolation
			DO k = 1, nzm
				DO i = 1, nx
					iwrf = i0 + i
					
					! U, V
					u(k) = u(k) + u3d(i,k,j) * 0.5 * ( mweight(i-1,k,j) + mweight(i,k,j) )
					v(k) = v(k) + v3d(i,k,j) * 0.5 * ( mweight(i,k,j-1) + mweight(i,k,j) )
					
					! Scalars
					km1 = km1s(i,k,j)
					k00 = k00s(i,k,j)
					kp1 = kp1s(i,k,j)
					kp2 = kp2s(i,k,j)
					
					! PRES
					! Pressure: mass-weighting is not required (See Canuto 1997 Atrophys. J.)
					pres(k) = pres(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							ps_eta(i,km1,j), ps_eta(i,k00,j), ps_eta(i,kp1,j), ps_eta(i,kp2,j) )
					
					! DSE, LWSE, LIWSE, MSE, VDSE
					dse(k) = dse(k)  + dse3d(i,k,j)  * mweight(i,k,j)
					
					lwse(k) = lwse(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							lwse_eta(i,km1), lwse_eta(i,k00), lwse_eta(i,kp1), lwse_eta(i,kp2) ) &
							* mweight(i,k,j)
					
					liwse(k) = liwse(k) + liwse3d(i,k,j) * mweight(i,k,j)
					mse(k)   = mse(k)   + mse3d(i,k,j)   * mweight(i,k,j)
					vdse(k)  = vdse(k)  + vdse3d(i,k,j)  * mweight(i,k,j)
					
					! SMSE
					smse(k) = smse(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							smse_eta(i,km1), smse_eta(i,k00), smse_eta(i,kp1), smse_eta(i,kp2) ) &
							* mweight(i,k,j)
					
					! TP, TPL, TPLI, TPE, TPV, TABS
					tp(k) = tp(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							tp_eta(i,km1,j), tp_eta(i,k00,j), tp_eta(i,kp1,j), tp_eta(i,kp2,j) ) &
							* mweight(i,k,j)
					
					tpl(k)  = tpl(k)  + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							tpl_eta(i,km1),  tpl_eta(i,k00),  tpl_eta(i,kp1),  tpl_eta(i,kp2) ) &
							* mweight(i,k,j)
					
					tpli(k) = tpli(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							tpli_eta(i,km1), tpli_eta(i,k00), tpli_eta(i,kp1), tpli_eta(i,kp2) ) &
							* mweight(i,k,j)
					
					tpe(k)  = tpe(k)  + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							tpe_eta(i,km1),  tpe_eta(i,k00),  tpe_eta(i,kp1),  tpe_eta(i,kp2) ) &
							* mweight(i,k,j)
					
					tpv(k)  = tpv(k)  + tpv3d(i,k,j)  * mweight(i,k,j)
					tabs(k) = tabs(k) + tabs3d(i,k,j) * mweight(i,k,j)
					
					! QT, QV, QC, QI, QR, QS, QG, QN, QP, QCL, QCI, QCT, QSAT, RELH
					qt(k) = qt(k) + qt3d(i,k,j) * mweight(i,k,j)
					
					qv(k) = qv(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							qv_eta(i,km1,j), qv_eta(i,k00,j), qv_eta(i,kp1,j), qv_eta(i,kp2,j) ) &
							* mweight(i,k,j)
					
					qc(k) = qc(k) + qc3d(i,k,j) * mweight(i,k,j)
					qi(k) = qi(k) + qi3d(i,k,j) * mweight(i,k,j)
					
					qr(k) = qr(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							qr_eta(i,km1,j), qr_eta(i,k00,j), qr_eta(i,kp1,j), qr_eta(i,kp2,j) ) &
							* mweight(i,k,j)
					
					qs(k) = qs(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							qs_eta(i,km1,j), qs_eta(i,k00,j), qs_eta(i,kp1,j), qs_eta(i,kp2,j) ) &
							* mweight(i,k,j)
					
					qg(k) = qg(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							qg_eta(i,km1,j), qg_eta(i,k00,j), qg_eta(i,kp1,j), qg_eta(i,kp2,j) ) &
							* mweight(i,k,j)
					
					qn(k)  = qn(k)  + qn3d(i,k,j)  * mweight(i,k,j)
					qp(k)  = qp(k)  + qp3d(i,k,j)  * mweight(i,k,j)
					
					qcl(k) = qcl(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j),  zs_eta(i,k00,j),  zs_eta(i,kp1,j),  zs_eta(i,kp2,j), &
							qcl_eta(i,km1,j), qcl_eta(i,k00,j), qcl_eta(i,kp1,j), qcl_eta(i,kp2,j) ) &
							* mweight(i,k,j)
					
					qci(k) = qci(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j),  zs_eta(i,k00,j),  zs_eta(i,kp1,j),  zs_eta(i,kp2,j), &
							qci_eta(i,km1,j), qci_eta(i,k00,j), qci_eta(i,kp1,j), qci_eta(i,kp2,j) ) &
							* mweight(i,k,j)
					
					qct(k) = qct(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							qct_eta(i,km1),  qct_eta(i,k00),  qct_eta(i,kp1),  qct_eta(i,kp2) ) &
							* mweight(i,k,j)
					
					qsat(k) = qsat(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							qsat_eta(i,km1), qsat_eta(i,k00), qsat_eta(i,kp1), qsat_eta(i,kp2) ) &
							* mweight(i,k,j)
					
					relh(k) = relh(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							relh_eta(i,km1), relh_eta(i,k00), relh_eta(i,kp1), relh_eta(i,kp2) ) &
							* mweight(i,k,j)
					
					! ETA
					eta(k) = eta(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j),  zs_eta(i,k00,j),  zs_eta(i,kp1,j),  zs_eta(i,kp2,j), &
							eta_eta(i,km1,j), eta_eta(i,k00,j), eta_eta(i,kp1,j), eta_eta(i,kp2,j) )
					
					! SSR, QNC, QNI, QNR, QNS, QNG, QNN, QNT
					eta(k) = eta(k) + eta3d(i,k,j) * mweight(i,k,j)
					ssr(k) = ssr(k) + ssr3d(i,k,j) * mweight(i,k,j)
					qnc(k) = qnc(k) + qnc3d(i,k,j) * mweight(i,k,j)
					qni(k) = qni(k) + qni3d(i,k,j) * mweight(i,k,j)
					qnr(k) = qnr(k) + qnr3d(i,k,j) * mweight(i,k,j)
					qns(k) = qns(k) + qns3d(i,k,j) * mweight(i,k,j)
					qng(k) = qng(k) + qng3d(i,k,j) * mweight(i,k,j)
					qnn(k) = qnn(k) + qnn3d(i,k,j) * mweight(i,k,j)
					
					qnt(k) = qnt(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							qnt_eta(i,km1),  qnt_eta(i,k00),  qnt_eta(i,kp1),  qnt_eta(i,kp2) ) &
							* mweight(i,k,j)
					
					! COND, EVAP, NCOND, DEPS, SUBL, NDEPS, TCOND
					! The unit of ####_rate_eta is kg/kg/s
					cond_rate(k) = cond_rate(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							cond_rate_eta(iwrf,km1,jwrf), cond_rate_eta(iwrf,k00,jwrf), &
							cond_rate_eta(iwrf,kp1,jwrf), cond_rate_eta(iwrf,kp2,jwrf) ) * mweight(i,k,j)
					
					evap_rate(k) = evap_rate(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							evap_rate_eta(iwrf,km1,jwrf), evap_rate_eta(iwrf,k00,jwrf), &
							evap_rate_eta(iwrf,kp1,jwrf), evap_rate_eta(iwrf,kp2,jwrf) ) * mweight(i,k,j)
					
					ncond_rate(k) = ncond_rate(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							ncond_rate_eta(i,km1), ncond_rate_eta(i,k00), &
							ncond_rate_eta(i,kp1), ncond_rate_eta(i,kp2) ) * mweight(i,k,j)
					
					deps_rate(k) = deps_rate(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							deps_rate_eta(iwrf,km1,jwrf), deps_rate_eta(iwrf,k00,jwrf), &
							deps_rate_eta(iwrf,kp1,jwrf), deps_rate_eta(iwrf,kp2,jwrf) ) * mweight(i,k,j)
					
					subl_rate(k) = subl_rate(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							subl_rate_eta(iwrf,km1,jwrf), subl_rate_eta(iwrf,k00,jwrf), &
							subl_rate_eta(iwrf,kp1,jwrf), subl_rate_eta(iwrf,kp2,jwrf) ) * mweight(i,k,j)
					
					ndeps_rate(k) = ndeps_rate(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							ndeps_rate_eta(i,km1), ndeps_rate_eta(i,k00), &
							ndeps_rate_eta(i,kp1), ndeps_rate_eta(i,kp2) ) * mweight(i,k,j)
					
					tcond_rate(k) = tcond_rate(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							tcond_rate_eta(i,km1), tcond_rate_eta(i,k00), &
							tcond_rate_eta(i,kp1), tcond_rate_eta(i,kp2) ) * mweight(i,k,j)
				ENDDO   ! i
			ENDDO   ! k
		ENDDO   ! j
		
		CALL hbuf_put( 'PRES', pres,factor_xy )   ! mb
		CALL hbuf_put( 'U', u, factor_xy )   ! m/s
		CALL hbuf_put( 'V', v, factor_xy )   ! m/s
		CALL hbuf_put( 'DSE',  dse,  factor_xy)    ! K
		CALL hbuf_put( 'LWSE', lwse, factor_xy )   ! K
		CALL hbuf_put( 'LIWSE', liwse, factor_xy )   ! K
		CALL hbuf_put( 'MSE',  mse,  factor_xy )   ! K
		CALL hbuf_put( 'VDSE', vdse, factor_xy )   ! K
		CALL hbuf_put( 'SMSE', smse, factor_xy )   ! K
		CALL hbuf_put( 'TP',   tp,  factor_xy )   ! K
		CALL hbuf_put( 'TPL',  tpl, factor_xy )   ! K
		CALL hbuf_put( 'TPLI', tpli, factor_xy )  ! K
		CALL hbuf_put( 'TPE',  tpe, factor_xy )   ! K
		CALL hbuf_put( 'TPV',  tpv, factor_xy )   ! K
		CALL hbuf_put( 'TABS', tabs, factor_xy )  ! K
		CALL hbuf_put( 'QT', qt, 1.E3*factor_xy )   ! kg/kg ==> g/kg
		CALL hbuf_put( 'QV', qv, 1.E3*factor_xy )   ! kg/kg ==> g/kg
		CALL hbuf_put( 'QC', qc, 1.E3*factor_xy )   ! kg/kg ==> g/kg
		CALL hbuf_put( 'QI', qi, 1.E3*factor_xy )   ! kg/kg ==> g/kg
		CALL hbuf_put( 'QR', qr, 1.E3*factor_xy )   ! kg/kg ==> g/kg
		CALL hbuf_put( 'QS', qs, 1.E3*factor_xy )   ! kg/kg ==> g/kg
		CALL hbuf_put( 'QG', qg, 1.E3*factor_xy )   ! kg/kg ==> g/kg
		CALL hbuf_put( 'QN', qn, 1.E3*factor_xy )   ! kg/kg ==> g/kg
		CALL hbuf_put( 'QP', qp, 1.E3*factor_xy )   ! kg/kg ==> g/kg
		CALL hbuf_put( 'QCL', qcl, 1.E3*factor_xy )   ! kg/kg ==> g/kg
		CALL hbuf_put( 'QCI', qci, 1.E3*factor_xy )   ! kg/kg ==> g/kg
		CALL hbuf_put( 'QCT', qct, 1.E3*factor_xy )   ! kg/kg ==> g/kg
		CALL hbuf_put( 'QSAT', qsat, 1.E3*factor_xy )   ! kg/kg ==> g/kg
		CALL hbuf_put( 'RELH', relh, 1.E2*factor_xy )   ! fraction ==> %
		CALL hbuf_put( 'ETA', eta, 1.E3*factor_xy )   ! kg/kg ==> g/kg
		CALL hbuf_put( 'SSR', ssr, factor_xy )   ! %
		CALL hbuf_put( 'QNC', qnc, factor_xy )   ! #/mg
		CALL hbuf_put( 'QNI', qni, factor_xy )   ! #/mg
		CALL hbuf_put( 'QNR', qnr, factor_xy )   ! #/mg
		CALL hbuf_put( 'QNS', qns, factor_xy )   ! #/mg
		CALL hbuf_put( 'QNG', qng, factor_xy )   ! #/mg
		CALL hbuf_put( 'QNN', qnn, factor_xy )   ! #/mg
		CALL hbuf_put( 'QNT', qnt, factor_xy )   ! #/mg
		CALL hbuf_put( 'COND',  cond_rate,  factor_xy*3600000. )   ! kg/kg/s ==> g/kg/hour
		CALL hbuf_put( 'EVAP',  evap_rate,  factor_xy*3600000. )   ! kg/kg/s ==> g/kg/hour
		CALL hbuf_put( 'NCOND', ncond_rate, factor_xy*3600000. )   ! kg/kg/s ==> g/kg/hour
		CALL hbuf_put( 'DEPS',  deps_rate,  factor_xy*3600000. )   ! kg/kg/s ==> g/kg/hour
		CALL hbuf_put( 'SUBL',  subl_rate,  factor_xy*3600000. )   ! kg/kg/s ==> g/kg/hour
		CALL hbuf_put( 'NDEPS', ndeps_rate, factor_xy*3600000. )   ! kg/kg/s ==> g/kg/hour
		CALL hbuf_put( 'TCOND', tcond_rate, factor_xy*3600000. )   ! kg/kg/s ==> g/kg/hour
		
		
		!========== Flux ==========
		! Resolved scale flux is the horizontal mean of flux computed for vertical advection.
		! SGS flux is the vertical diffusive flux computed for diffusion.
		! Both fluxes are defined at the W level.
		! twle_eta and twsb_eta have been converted from TP to TABS.
		uw(:) = 0.
		uwle(:) = 0.
		uwsb(:) = 0.
		vw(:) = 0.
		vwle(:) = 0.
		vwsb(:) = 0.
		lwsef(:) = 0.
		lwsefr(:) = 0.
		lwsefs(:) = 0.
		tplf(:)  = 0.
		tplfr(:) = 0.
		tplfs(:) = 0.
		tpvf(:)  = 0.
		tpvfr(:) = 0.
		tpvfs(:) = 0.
		qtf(:)  = 0.
		qtfr(:) = 0.
		qtfs(:) = 0.
		qpf(:)  = 0.
		qpfr(:) = 0.
		qpfs(:) = 0.
		qcf(:)  = 0.
		qcfr(:) = 0.
		qcfs(:) = 0.
		qif(:)  = 0.
		qifr(:) = 0.
		qifs(:) = 0.
		qrf(:)  = 0.
		qrfr(:) = 0.
		qrfs(:) = 0.
		qsf(:)  = 0.
		qsfr(:) = 0.
		qsfs(:) = 0.
		qgf(:)  = 0.
		qgfr(:) = 0.
		qgfs(:) = 0.
		prec(:) = 0.
		precl(:) = 0.
		preci(:) = 0.
		vfmi(:)  = 0.
		DO j = 1, ny
			jwrf = j0 + j
			! Diagnose flux into 2D array
			DO k = 1, nzm
				DO i = 1, nx
					iwrf = i0 + i
					! U
					uw_eta(i,k) = uwle_eta(iwrf,k,jwrf) + uwsb_eta(iwrf,k,jwrf)
					! V
					vw_eta(i,k) = vwle_eta(iwrf,k,jwrf) + vwsb_eta(iwrf,k,jwrf)
					! LWSE
					lwsefr_eta(i,k) = cp * exnerw_eta(i,k,j) * twle_eta(iwrf,k,jwrf) &
					     - lcond * ( mwle_eta(iwrf,k,jwrf,p_qc) + mwle_eta(iwrf,k,jwrf,p_qr) )
					lwsefs_eta(i,k) = cp * exnerw_eta(i,k,j) * twsb_eta(iwrf,k,jwrf) &
					     - lcond * ( mwsb_eta(iwrf,k,jwrf,p_qc) + mwsb_eta(iwrf,k,jwrf,p_qr) )
					lwsef_eta(i,k) = lwsefr_eta(i,k) + lwsefs_eta(i,k)
					! LIWSE
					liwsefr_eta(i,k) = lwsefr_eta(i,k) - lsub * ( mwle_eta(iwrf,k,jwrf,p_qi) &
					                 + mwle_eta(iwrf,k,jwrf,p_qs) + mwle_eta(iwrf,k,jwrf,p_qg) )
					liwsefs_eta(i,k) = lwsefs_eta(i,k) - lsub * ( mwsb_eta(iwrf,k,jwrf,p_qi) &
					                 + mwsb_eta(iwrf,k,jwrf,p_qs) + mwsb_eta(iwrf,k,jwrf,p_qg) )
					liwsef_eta(i,k) = liwsefr_eta(i,k) + liwsefs_eta(i,k)
					! TPL
					tplfr_eta(i,k) = cp * twle_eta(iwrf,k,jwrf) - lcond / exnerw_eta(i,k,j) &
					           * ( mwle_eta(iwrf,k,jwrf,p_qc) + mwle_eta(iwrf,k,jwrf,p_qr) )
					tplfs_eta(i,k) = cp * twsb_eta(iwrf,k,jwrf) - lcond / exnerw_eta(i,k,j) &
					           * ( mwsb_eta(iwrf,k,jwrf,p_qc) + mwsb_eta(iwrf,k,jwrf,p_qr) )
					tplf_eta(i,k) = tplfr_eta(i,k) + tplfs_eta(i,k)
					! QT
					qtfr_eta(i,k) = lcond * ( mwle_eta(iwrf,k,jwrf,p_qv) + mwle_eta(iwrf,k,jwrf,p_qc) ) &
					              + lsub * mwle_eta(iwrf,k,jwrf,p_qi)
					qtfs_eta(i,k) = lcond * ( mwsb_eta(iwrf,k,jwrf,p_qv) + mwsb_eta(iwrf,k,jwrf,p_qc) ) &
					              + lsub * mwsb_eta(iwrf,k,jwrf,p_qi)
					qtf_eta(i,k) = qtfr_eta(i,k) + qtfs_eta(i,k)
					! QP
					qpfr_eta(i,k) = lcond * mwle_eta(iwrf,k,jwrf,p_qr) &
					              + lsub * ( mwle_eta(iwrf,k,jwrf,p_qs) + mwle_eta(iwrf,k,jwrf,p_qg) )
					qpfs_eta(i,k) = lcond * mwsb_eta(iwrf,k,jwrf,p_qr) &
					              + lsub * ( mwsb_eta(iwrf,k,jwrf,p_qs) + mwsb_eta(iwrf,k,jwrf,p_qg) )
					qpf_eta(i,k) = qpfr_eta(i,k) + qpfs_eta(i,k)
					! QC, QI, QR, QS, QG
					qcf_eta(i,k) = mwle_eta(iwrf,k,jwrf,p_qc) + mwsb_eta(iwrf,k,jwrf,p_qc)
					qif_eta(i,k) = mwle_eta(iwrf,k,jwrf,p_qi) + mwsb_eta(iwrf,k,jwrf,p_qi)
					qrf_eta(i,k) = mwle_eta(iwrf,k,jwrf,p_qr) + mwsb_eta(iwrf,k,jwrf,p_qr)
					qsf_eta(i,k) = mwle_eta(iwrf,k,jwrf,p_qs) + mwsb_eta(iwrf,k,jwrf,p_qs)
					qgf_eta(i,k) = mwle_eta(iwrf,k,jwrf,p_qg) + mwsb_eta(iwrf,k,jwrf,p_qg)
					! PRECIP
					precip_eta(i,k) = precipl_eta(iwrf,k,jwrf) + precipi_eta(iwrf,k,jwrf)
				ENDDO
			ENDDO
			! Flux related with virtual temperature
			! TPV
			DO k = 2, nzm
				DO i = 1, nx
					iwrf = i0 + i
					tpvfr_eta(i,k) = cp * twle_eta(iwrf,k,jwrf) &
					           * ( 1. + ( epsv * mwle_eta(iwrf,k,jwrf,p_qv) &
					           - mwle_eta(iwrf,k,jwrf,p_qc) - mwle_eta(iwrf,k,jwrf,p_qr) &
					           - mwle_eta(iwrf,k,jwrf,p_qi) - mwle_eta(iwrf,k,jwrf,p_qs) &
					           - mwle_eta(iwrf,k,jwrf,p_qg) ) &
					       / ( SIGN(1.,wadv_eta(iwrf,k,jwrf))*MAX(ABS(wadv_eta(iwrf,k,jwrf)),1.E-15) ) )
					tpvfs_eta(i,k) = cp * twsb_eta(iwrf,k,jwrf) - tkhw_eta(i,k,j) &
					     * ( tp_eta(i,k,j)*( epsv*qv_eta(i,k,j)-qn_eta(i,k,j)-qp_eta(i,k,j) ) &
					       - tp_eta(i,k-1,j)*( epsv*qv_eta(i,k-1,j)-qn_eta(i,k-1,j)-qp_eta(i,k-1,j) ) ) &
					     / ( zs_eta(i,k,j) - zs_eta(i,k-1,j) )
					tpvf_eta(i,k) = tpvfr_eta(i,k) + tpvfs_eta(i,k)
				ENDDO
			ENDDO
			! Add surface flux to tpvfs. Assume SGS [[w''tp''qv'']]=0
			k = 1
			DO i = 1, nx
				iwrf = i0 + i
				tpvfr_eta(i,k) = 0.
				tpvfs_eta(i,k) = cp * ( ( 1. + epsv * qv_eta(i,k,j) ) * twsb_eta(iwrf,k,jwrf) &
				               + epsv * tp_eta(i,k,j) * mwsb_eta(iwrf,k,jwrf,p_qv) )
				tpvf_eta(i,k)  = tpvfr_eta(i,k) + tpvfs_eta(i,k)
			ENDDO
			
			! Interpolation
			DO k = 1, nzm
				DO i = 1, nx
					iwrf = i0 + i
					
					! U flux (m2/s2) @ UW
					km1 = km1uw(i,k,j)
					k00 = k00uw(i,k,j)
					kp1 = kp1uw(i,k,j)
					kp2 = kp2uw(i,k,j)
					! Total
					uw(k) = uw(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zuw_eta(i,km1,j), zuw_eta(i,k00,j), zuw_eta(i,kp1,j), zuw_eta(i,kp2,j), &
							uw_eta(i,km1),    uw_eta(i,k00),    uw_eta(i,kp1),    uw_eta(i,kp2) ) &
							* mweight(i,k,j)
					! Resolved
					uwle(k) = uwle(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zuw_eta(i,km1,j), zuw_eta(i,k00,j), zuw_eta(i,kp1,j), zuw_eta(i,kp2,j), &
							uwle_eta(iwrf,km1,jwrf), uwle_eta(iwrf,k00,jwrf), &
							uwle_eta(iwrf,kp1,jwrf), uwle_eta(iwrf,kp2,jwrf) ) * mweight(i,k,j)
					! SGS
					uwsb(k) = uwsb(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zuw_eta(i,km1,j), zuw_eta(i,k00,j), zuw_eta(i,kp1,j), zuw_eta(i,kp2,j), &
							uwsb_eta(iwrf,km1,jwrf), uwsb_eta(iwrf,k00,jwrf), &
							uwsb_eta(iwrf,kp1,jwrf), uwsb_eta(iwrf,kp2,jwrf) ) * mweight(i,k,j)
					
					! V flux (m2/s2) @ VW
					km1 = km1vw(i,k,j)
					k00 = k00vw(i,k,j)
					kp1 = kp1vw(i,k,j)
					kp2 = kp2vw(i,k,j)
					! Total
					vw(k) = vw(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zvw_eta(i,km1,j), zvw_eta(i,k00,j), zvw_eta(i,kp1,j), zvw_eta(i,kp2,j), &
							vw_eta(i,km1),    vw_eta(i,k00),    vw_eta(i,kp1),    vw_eta(i,kp2) ) &
							* mweight(i,k,j)
					! Resolve
					vwle(k) = vwle(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zvw_eta(i,km1,j), zvw_eta(i,k00,j), zvw_eta(i,kp1,j), zvw_eta(i,kp2,j), &
							vwle_eta(iwrf,km1,jwrf), vwle_eta(iwrf,k00,jwrf), &
							vwle_eta(iwrf,kp1,jwrf), vwle_eta(iwrf,kp2,jwrf) ) * mweight(i,k,j)
					! SGS
					vwsb(k) = vwsb(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zvw_eta(i,km1,j),  zvw_eta(i,k00,j),  zvw_eta(i,kp1,j),  zvw_eta(i,kp2,j), &
							vwsb_eta(iwrf,km1,jwrf), vwsb_eta(iwrf,k00,jwrf), &
							vwsb_eta(iwrf,kp1,jwrf), vwsb_eta(iwrf,kp2,jwrf) ) * mweight(i,k,j)
					
			!!		! TPV flux (W/m2)
			!!		! Perturbation of TPV has to be computed here.
			!!		km1 = km1s(i,k,j)
			!!		k00 = k00s(i,k,j)
			!!		kp1 = kp1s(i,k,j)
			!!		kp2 = kp2s(i,k,j)
			!!		fprime = mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
			!!				zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
			!!				tpv_eta(i,km1,j)-tpv0_eta(i,km1,j), tpv_eta(i,k00,j)-tpv0_eta(i,k00,j), &
			!!				tpv_eta(i,kp1,j)-tpv0_eta(i,kp1,j), tpv_eta(i,kp2,j)-tpv0_eta(i,kp2,j) )
			!!		tpvfr(k) = tpvfr(k) + cp * rho3d(i,k,j) * w3d(i,k,j) * fprime
					
					! S flux (W/m2)
					km1 = km1w(i,k,j)
					k00 = k00w(i,k,j)
					kp1 = kp1w(i,k,j)
					kp2 = kp2w(i,k,j)
					
					! LWSE
					! Total
					lwsef(k) = lwsef(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zw_eta(i,km1,j),  zw_eta(i,k00,j),  zw_eta(i,kp1,j),  zw_eta(i,kp2,j), &
							lwsef_eta(i,km1), lwsef_eta(i,k00), lwsef_eta(i,kp1), lwsef_eta(i,kp2) ) &
							* rho3d(i,k,j)
					! Resolved
					lwsefr(k) = lwsefr(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zw_eta(i,km1,j), zw_eta(i,k00,j), zw_eta(i,kp1,j), zw_eta(i,kp2,j), &
							lwsefr_eta(i,km1), lwsefr_eta(i,k00), lwsefr_eta(i,kp1), lwsefr_eta(i,kp2) ) &
							* rho3d(i,k,j)
					! SGS
					lwsefs(k) = lwsefs(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zw_eta(i,km1,j), zw_eta(i,k00,j), zw_eta(i,kp1,j), zw_eta(i,kp2,j), &
							lwsefs_eta(i,km1), lwsefs_eta(i,k00), lwsefs_eta(i,kp1), lwsefs_eta(i,kp2) ) &
							* rho3d(i,k,j)
					
					! LIWSE
					! Total
					liwsef(k) = liwsef(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zw_eta(i,km1,j),   zw_eta(i,k00,j),   zw_eta(i,kp1,j),   zw_eta(i,kp2,j), &
							liwsef_eta(i,km1), liwsef_eta(i,k00), liwsef_eta(i,kp1), liwsef_eta(i,kp2) ) &
							* rho3d(i,k,j)
					! Resolved
					liwsefr(k) = liwsefr(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zw_eta(i,km1,j),   zw_eta(i,k00,j),   zw_eta(i,kp1,j),   zw_eta(i,kp2,j), &
							liwsefr_eta(i,km1),liwsefr_eta(i,k00),liwsefr_eta(i,kp1),liwsefr_eta(i,kp2) ) &
							* rho3d(i,k,j)
					! SGS
					liwsefs(k) = liwsefs(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zw_eta(i,km1,j),   zw_eta(i,k00,j),   zw_eta(i,kp1,j),   zw_eta(i,kp2,j), &
							liwsefs_eta(i,km1),liwsefs_eta(i,k00),liwsefs_eta(i,kp1),liwsefs_eta(i,kp2) ) &
							* rho3d(i,k,j)
					
					! TPL
					! Total
					tplf(k) = tplf(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zw_eta(i,km1,j), zw_eta(i,k00,j), zw_eta(i,kp1,j), zw_eta(i,kp2,j), &
							tplf_eta(i,km1), tplf_eta(i,k00), tplf_eta(i,kp1), tplf_eta(i,kp2) ) &
							* rho3d(i,k,j)
					! Resolved
					tplfr(k) = tplfr(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zw_eta(i,km1,j), zw_eta(i,k00,j), zw_eta(i,kp1,j), zw_eta(i,kp2,j), &
							tplfr_eta(i,km1), tplfr_eta(i,k00), tplfr_eta(i,kp1), tplfr_eta(i,kp2) ) &
							* rho3d(i,k,j)
					! SGS
					tplfs(k) = tplfs(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zw_eta(i,km1,j), zw_eta(i,k00,j), zw_eta(i,kp1,j), zw_eta(i,kp2,j), &
							tplfs_eta(i,km1), tplfs_eta(i,k00), tplfs_eta(i,kp1), tplfs_eta(i,kp2) ) &
							* rho3d(i,k,j)
					
					! TPV
					! Total
					tpvf(k) = tpvf(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zw_eta(i,km1,j), zw_eta(i,k00,j), zw_eta(i,kp1,j), zw_eta(i,kp2,j), &
							tpvf_eta(i,km1), tpvf_eta(i,k00), tpvf_eta(i,kp1), tpvf_eta(i,kp2) ) &
							* rho3d(i,k,j)
					! Resolved
					tpvfr(k) = tpvfr(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zw_eta(i,km1,j),  zw_eta(i,k00,j),  zw_eta(i,kp1,j),  zw_eta(i,kp2,j), &
							tpvfr_eta(i,km1), tpvfr_eta(i,k00), tpvfr_eta(i,kp1), tpvfr_eta(i,kp2) ) &
							* rho3d(i,k,j)
					! SGS
					tpvfs(k) = tpvfs(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zw_eta(i,km1,j),  zw_eta(i,k00,j),  zw_eta(i,kp1,j),  zw_eta(i,kp2,j), &
							tpvfs_eta(i,km1), tpvfs_eta(i,k00), tpvfs_eta(i,kp1), tpvfs_eta(i,kp2) ) &
							* rho3d(i,k,j)
					
					! QT
					! Total
					qtf(k) = qtf(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zw_eta(i,km1,j), zw_eta(i,k00,j), zw_eta(i,kp1,j), zw_eta(i,kp2,j), &
							qtf_eta(i,km1),  qtf_eta(i,k00),  qtf_eta(i,kp1),  qtf_eta(i,kp2) ) &
							* rho3d(i,k,j)
					! Resolved
					qtfr(k) = qtfr(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zw_eta(i,km1,j), zw_eta(i,k00,j), zw_eta(i,kp1,j), zw_eta(i,kp2,j), &
							qtfr_eta(i,km1), qtfr_eta(i,k00), qtfr_eta(i,kp1), qtfr_eta(i,kp2) ) &
							* rho3d(i,k,j)
					! SGS
					qtfs(k) = qtfs(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zw_eta(i,km1,j), zw_eta(i,k00,j), zw_eta(i,kp1,j), zw_eta(i,kp2,j), &
							qtfs_eta(i,km1), qtfs_eta(i,k00), qtfs_eta(i,kp1), qtfs_eta(i,kp2) ) &
							* rho3d(i,k,j)
					
					! QP
					! Total
					qpf(k) = qpf(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zw_eta(i,km1,j), zw_eta(i,k00,j), zw_eta(i,kp1,j), zw_eta(i,kp2,j), &
							qpf_eta(i,km1),  qpf_eta(i,k00),  qpf_eta(i,kp1),  qpf_eta(i,kp2) ) &
							* rho3d(i,k,j)
					! Resolved
					qpfr(k) = qpfr(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zw_eta(i,km1,j), zw_eta(i,k00,j), zw_eta(i,kp1,j), zw_eta(i,kp2,j), &
							qpfr_eta(i,km1), qpfr_eta(i,k00), qpfr_eta(i,kp1), qpfr_eta(i,kp2) ) &
							* rho3d(i,k,j)
					! SGS
					qpfs(k) = qpfs(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zw_eta(i,km1,j), zw_eta(i,k00,j), zw_eta(i,kp1,j), zw_eta(i,kp2,j), &
							qpfs_eta(i,km1), qpfs_eta(i,k00), qpfs_eta(i,kp1), qpfs_eta(i,kp2) ) &
							* rho3d(i,k,j)
					
					! QC
					! Total
					qcf(k) = qcf(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zw_eta(i,km1,j), zw_eta(i,k00,j), zw_eta(i,kp1,j), zw_eta(i,kp2,j), &
							qcf_eta(i,km1),  qcf_eta(i,k00),  qcf_eta(i,kp1),  qcf_eta(i,kp2) ) &
							* lcond * rho3d(i,k,j)
					! Resolved
					qcfr(k) = qcfr(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zw_eta(i,km1,j), zw_eta(i,k00,j), zw_eta(i,kp1,j), zw_eta(i,kp2,j), &
							mwle_eta(iwrf,km1,jwrf,p_qc), mwle_eta(iwrf,k00,jwrf,p_qc), &
							mwle_eta(iwrf,kp1,jwrf,p_qc), mwle_eta(iwrf,kp2,jwrf,p_qc) ) &
							* lcond * rho3d(i,k,j)
					! SGS
					qcfs(k) = qcfs(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zw_eta(i,km1,j), zw_eta(i,k00,j), zw_eta(i,kp1,j), zw_eta(i,kp2,j), &
							mwsb_eta(iwrf,km1,jwrf,p_qc), mwsb_eta(iwrf,k00,jwrf,p_qc), &
							mwsb_eta(iwrf,kp1,jwrf,p_qc), mwsb_eta(iwrf,kp2,jwrf,p_qc) ) &
							* lcond * rho3d(i,k,j)
					
					! QI
					! Total
					qif(k) = qif(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zw_eta(i,km1,j), zw_eta(i,k00,j), zw_eta(i,kp1,j), zw_eta(i,kp2,j), &
							qif_eta(i,km1),  qif_eta(i,k00),  qif_eta(i,kp1),  qif_eta(i,kp2) ) &
							* lsub * rho3d(i,k,j)
					! Resolved
					qifr(k) = qifr(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zw_eta(i,km1,j), zw_eta(i,k00,j), zw_eta(i,kp1,j), zw_eta(i,kp2,j), &
							mwle_eta(iwrf,km1,jwrf,p_qi), mwle_eta(iwrf,k00,jwrf,p_qi), &
							mwle_eta(iwrf,kp1,jwrf,p_qi), mwle_eta(iwrf,kp2,jwrf,p_qi) ) &
							* lsub * rho3d(i,k,j)
					! SGS
					qifs(k) = qifs(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zw_eta(i,km1,j), zw_eta(i,k00,j), zw_eta(i,kp1,j), zw_eta(i,kp2,j), &
							mwsb_eta(iwrf,km1,jwrf,p_qi), mwsb_eta(iwrf,k00,jwrf,p_qi), &
							mwsb_eta(iwrf,kp1,jwrf,p_qi), mwsb_eta(iwrf,kp2,jwrf,p_qi) ) &
							* lsub * rho3d(i,k,j)
					
					! QR
					! Total
					qrf(k) = qrf(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zw_eta(i,km1,j), zw_eta(i,k00,j), zw_eta(i,kp1,j), zw_eta(i,kp2,j), &
							qrf_eta(i,km1),  qrf_eta(i,k00),  qrf_eta(i,kp1),  qrf_eta(i,kp2) ) &
							* lcond * rho3d(i,k,j)
					! Resolved
					qrfr(k) = qrfr(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zw_eta(i,km1,j), zw_eta(i,k00,j), zw_eta(i,kp1,j), zw_eta(i,kp2,j), &
							mwle_eta(iwrf,km1,jwrf,p_qr), mwle_eta(iwrf,k00,jwrf,p_qr), &
							mwle_eta(iwrf,kp1,jwrf,p_qr), mwle_eta(iwrf,kp2,jwrf,p_qr) ) &
							* lcond * rho3d(i,k,j)
					! SGS
					qrfs(k) = qrfs(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zw_eta(i,km1,j), zw_eta(i,k00,j), zw_eta(i,kp1,j), zw_eta(i,kp2,j), &
							mwsb_eta(iwrf,km1,jwrf,p_qr), mwsb_eta(iwrf,k00,jwrf,p_qr), &
							mwsb_eta(iwrf,kp1,jwrf,p_qr), mwsb_eta(iwrf,kp2,jwrf,p_qr)) &
							* lcond * rho3d(i,k,j)
					
					! QS
					! Total
					qsf(k) = qsf(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zw_eta(i,km1,j), zw_eta(i,k00,j), zw_eta(i,kp1,j), zw_eta(i,kp2,j), &
							qsf_eta(i,km1),  qsf_eta(i,k00),  qsf_eta(i,kp1),  qsf_eta(i,kp2) ) &
							* lsub * rho3d(i,k,j)
					! Resolved
					qsfr(k) = qsfr(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zw_eta(i,km1,j), zw_eta(i,k00,j), zw_eta(i,kp1,j), zw_eta(i,kp2,j), &
							mwle_eta(iwrf,km1,jwrf,p_qs), mwle_eta(iwrf,k00,jwrf,p_qs), &
							mwle_eta(iwrf,kp1,jwrf,p_qs), mwle_eta(iwrf,kp2,jwrf,p_qs) ) &
							* lsub * rho3d(i,k,j)
					! SGS
					qsfs(k) = qsfs(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zw_eta(i,km1,j), zw_eta(i,k00,j), zw_eta(i,kp1,j), zw_eta(i,kp2,j), &
							mwsb_eta(iwrf,km1,jwrf,p_qs), mwsb_eta(iwrf,k00,jwrf,p_qs), &
							mwsb_eta(iwrf,kp1,jwrf,p_qs), mwsb_eta(iwrf,kp2,jwrf,p_qs) ) &
							* lsub * rho3d(i,k,j)
					
					! QG
					! Total
					qgf(k) = qgf(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zw_eta(i,km1,j), zw_eta(i,k00,j), zw_eta(i,kp1,j), zw_eta(i,kp2,j), &
							qgf_eta(i,km1),  qgf_eta(i,k00),  qgf_eta(i,kp1),  qgf_eta(i,kp2) ) &
							* lsub * rho3d(i,k,j)
					! Resolved
					qgfr(k) = qgfr(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zw_eta(i,km1,j), zw_eta(i,k00,j), zw_eta(i,kp1,j), zw_eta(i,kp2,j), &
							mwle_eta(iwrf,km1,jwrf,p_qg), mwle_eta(iwrf,k00,jwrf,p_qg), &
							mwle_eta(iwrf,kp1,jwrf,p_qg), mwle_eta(iwrf,kp2,jwrf,p_qg) ) &
							* lsub * rho3d(i,k,j)
					! SGS
					qgfs(k) = qgfs(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zw_eta(i,km1,j), zw_eta(i,k00,j), zw_eta(i,kp1,j), zw_eta(i,kp2,j), &
							mwsb_eta(iwrf,km1,jwrf,p_qg), mwsb_eta(iwrf,k00,jwrf,p_qg), &
							mwsb_eta(iwrf,kp1,jwrf,p_qg), mwsb_eta(iwrf,kp2,jwrf,p_qg) ) &
							* lsub * rho3d(i,k,j)
					
					! PRECIP
					prec(k) = prec(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zw_eta(i,km1,j),   zw_eta(i,k00,j),   zw_eta(i,kp1,j),  zw_eta(i,kp2,j), &
							precip_eta(i,km1), precip_eta(i,k00), precip_eta(i,kp1), precip_eta(i,kp2) )
					
					! PRECIPL
					precl(k) = precl(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zw_eta(i,km1,j), zw_eta(i,k00,j), zw_eta(i,kp1,j), zw_eta(i,kp2,j), &
							precipl_eta(iwrf,km1,jwrf), precipl_eta(iwrf,k00,jwrf), &
							precipl_eta(iwrf,kp1,jwrf), precipl_eta(iwrf,kp2,jwrf) )
					
					! PRECIPI
					preci(k) = preci(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zw_eta(i,km1,j), zw_eta(i,k00,j), zw_eta(i,kp1,j), zw_eta(i,kp2,j), &
							precipi_eta(iwrf,km1,jwrf), precipi_eta(iwrf,k00,jwrf), &
							precipi_eta(iwrf,kp1,jwrf), precipi_eta(iwrf,kp2,jwrf) )
					
					! VFMI
					vfmi(k) = vfmi(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zw_eta(i,km1,j), zw_eta(i,k00,j), zw_eta(i,kp1,j), zw_eta(i,kp2,j), &
							vfmi_eta(iwrf,km1,jwrf), vfmi_eta(iwrf,k00,jwrf), &
							vfmi_eta(iwrf,kp1,jwrf), vfmi_eta(iwrf,kp2,jwrf) )
				ENDDO   ! k
			ENDDO   ! i
		ENDDO   ! j
		
		CALL hbuf_put( 'UW',  uw,   factor_xy )   ! m2/s2
		CALL hbuf_put( 'UWR', uwle, factor_xy )   ! m2/s2
		CALL hbuf_put( 'UWS', uwsb, factor_xy )   ! m2/s2
		CALL hbuf_put( 'VW',  vw,   factor_xy )   ! m2/s2
		CALL hbuf_put( 'VWR', vwle, factor_xy )   ! m2/s2
		CALL hbuf_put( 'VWS', vwsb, factor_xy )   ! m2/s2
		CALL hbuf_put( 'LWSEFLUX',  lwsef,  factor_xy )   ! W/m2
		CALL hbuf_put( 'LWSEFLUXR', lwsefr, factor_xy )   ! W/m2
		CALL hbuf_put( 'LWSEFLUXS', lwsefs, factor_xy )   ! W/m2
		CALL hbuf_put( 'LIWSEFLUX',  lwsef,  factor_xy )   ! W/m2
		CALL hbuf_put( 'LIWSEFLUXR', lwsefr, factor_xy )   ! W/m2
		CALL hbuf_put( 'LIWSEFLUXS', lwsefs, factor_xy )   ! W/m2
		CALL hbuf_put( 'TPLFLUX',  tplf,  factor_xy )   ! W/m2
		CALL hbuf_put( 'TPLFLUXR', tplfr, factor_xy )   ! W/m2
		CALL hbuf_put( 'TPLFLUXS', tplfs, factor_xy )   ! W/m2
		CALL hbuf_put( 'TPVFLUX',  tpvf,  factor_xy )   ! W/m2
		CALL hbuf_put( 'TPVFLUXR', tpvfr, factor_xy )   ! W/m2
		CALL hbuf_put( 'TPVFLUXS', tpvfs, factor_xy )   ! W/m2
		CALL hbuf_put( 'QTFLUX',  qtf,  factor_xy )   ! W/m2
		CALL hbuf_put( 'QTFLUXR', qtfr, factor_xy )   ! W/m2
		CALL hbuf_put( 'QTFLUXS', qtfs, factor_xy )   ! W/m2
		CALL hbuf_put( 'QPFLUX',  qpf,  factor_xy )   ! W/m2
		CALL hbuf_put( 'QPFLUXR', qpfr, factor_xy )   ! W/m2
		CALL hbuf_put( 'QPFLUXS', qpfs, factor_xy )   ! W/m2
		CALL hbuf_put( 'QCFLUX',  qcf,  factor_xy )   ! W/m2
		CALL hbuf_put( 'QCFLUXR', qcfr, factor_xy )   ! W/m2
		CALL hbuf_put( 'QCFLUXS', qcfs, factor_xy )   ! W/m2
		CALL hbuf_put( 'QIFLUX',  qif,  factor_xy )   ! W/m2
		CALL hbuf_put( 'QIFLUXR', qifr, factor_xy )   ! W/m2
		CALL hbuf_put( 'QIFLUXS', qifs, factor_xy )   ! W/m2
		CALL hbuf_put( 'QRFLUX',  qrf,  factor_xy )   ! W/m2
		CALL hbuf_put( 'QRFLUXR', qrfr, factor_xy )   ! W/m2
		CALL hbuf_put( 'QRFLUXS', qrfs, factor_xy )   ! W/m2
		CALL hbuf_put( 'QSFLUX',  qsf,  factor_xy )   ! W/m2
		CALL hbuf_put( 'QSFLUXR', qsfr, factor_xy )   ! W/m2
		CALL hbuf_put( 'QSFLUXS', qsfs, factor_xy )   ! W/m2
		CALL hbuf_put( 'QGFLUX',  qgf,  factor_xy )   ! W/m2
		CALL hbuf_put( 'QGFLUXR', qgfr, factor_xy )   ! W/m2
		CALL hbuf_put( 'QGFLUXS', qgfs, factor_xy )   ! W/m2
		CALL hbuf_put( 'PRECIP',  prec,  factor_xy )   ! mm/day
		CALL hbuf_put( 'PRECIPL', precl, factor_xy )   ! mm/day
		CALL hbuf_put( 'PRECIPI', preci, factor_xy )   ! mm/day
		CALL hbuf_put( 'VFMI', vfmi, factor_xy )   ! m/s
		
		
		!========== Higher order moments ==========
		! Resolved higher order moments
		! 1. At eta level, compute perturbation with high resolution mean profile
		! 2. Compute element of moment, e.g., u'u'
		! 3. Interpolate to zs
		! 4. Summation
		! This order most minimizes loss of moments due to interpolation.
		! w' is not computed since [s'w] = [s'w'] + [s'[w]] = [s'w'] + [s'][w] = [s'w'] 
		u2(:)   = 0.
		v2(:)   = 0.
		w2(:)   = 0.
		w3(:)   = 0.
		w4(:)   = 0.
		wskw(:) = 0.
		liwse2(:) = 0.
		liwse3(:) = 0.
		mse2(:)  = 0.
		mse3(:)  = 0.
		vdse2(:) = 0.
		vdse3(:) = 0.
		qt2(:)  = 0.
		qt3(:)  = 0.
		qc2(:)  = 0.
		qi2(:)  = 0.
		qcl2(:) = 0.
		qci2(:) = 0.
		qm2(:)  = 0.
		slqt(:) = 0.
		aup(:)  = 0.
		tke(:)  = 0.
		DO j = 1, ny
			! Diagnose higher-order moment into 2D array
			DO k = 1, nzm
				DO i = 1, nx
					! U, V, W
					u2_eta(i,k) = ( u_eta(i,k,j) - u0_eta(i,k,j) ) ** 2
					v2_eta(i,k) = ( v_eta(i,k,j) - v0_eta(i,k,j) ) ** 2
					w2_eta(i,k) = w_eta(i,k,j) ** 2
					w3_eta(i,k) = w_eta(i,k,j) ** 3
					w4_eta(i,k) = w_eta(i,k,j) ** 4
					
					! S
					km1 = km1s02s(i,k,j)
					k00 = k00s02s(i,k,j)
					kp1 = kp1s02s(i,k,j)
					kp2 = kp2s02s(i,k,j)
					! LIWSE
					f0 = mono_cubic_interp( zs_eta(i,k,j), km1, k00, kp1, kp2, &
							zsh(km1),    zsh(k00),    zsh(kp1),    zsh(kp2), &
							liwse0(km1), liwse0(k00), liwse0(kp1), liwse0(kp2) )
					liwse2_eta(i,k) = ( liwse_eta(i,k,j) - f0 ) ** 2
					liwse3_eta(i,k) = ( liwse_eta(i,k,j) - f0 ) ** 3
					fprime = liwse_eta(i,k,j) - f0   ! for correlation between LIWSE and QT
					! MSE
					f0 = mono_cubic_interp( zs_eta(i,k,j), km1, k00, kp1, kp2, &
							zsh(km1),  zsh(k00),  zsh(kp1),  zsh(kp2), &
							mse0(km1), mse0(k00), mse0(kp1), mse0(kp2) )
					mse2_eta(i,k) = ( mse_eta(i,k,j) - f0 ) ** 2
					mse3_eta(i,k) = ( mse_eta(i,k,j) - f0 ) ** 3
					! VDSE
					f0 = mono_cubic_interp( zs_eta(i,k,j), km1, k00, kp1, kp2, &
							zsh(km1),   zsh(k00),   zsh(kp1),   zsh(kp2), &
							vdse0(km1), vdse0(k00), vdse0(kp1), vdse0(kp2) )
					vdse2_eta(i,k) = ( vdse_eta(i,k,j) - f0 ) ** 2
					vdse3_eta(i,k) = ( vdse_eta(i,k,j) - f0 ) ** 3
					! QT
					f0 = mono_cubic_interp( zs_eta(i,k,j), km1, k00, kp1, kp2, &
						zsh(km1), zsh(k00), zsh(kp1), zsh(kp2), qt0(km1), qt0(k00), qt0(kp1), qt0(kp2) )
					qt2_eta(i,k) = ( qt_eta(i,k,j) - f0 ) ** 2
					qt3_eta(i,k) = ( qt_eta(i,k,j) - f0 ) ** 3
					! LIWSEQT
					slqt_eta(i,k) = fprime * ( qt_eta(i,k,j) - f0 )
					! QC
					f0 = mono_cubic_interp( zs_eta(i,k,j), km1, k00, kp1, kp2, &
						zsh(km1), zsh(k00), zsh(kp1), zsh(kp2), qc0(km1), qc0(k00), qc0(kp1), qc0(kp2) )
					qc2_eta(i,k) = ( qc_eta(i,k,j) - f0 ) ** 2
					! QI
					f0 = mono_cubic_interp( zs_eta(i,k,j), km1, k00, kp1, kp2, &
						zsh(km1), zsh(k00), zsh(kp1), zsh(kp2), qi0(km1), qi0(k00), qi0(kp1), qi0(kp2) )
					qi2_eta(i,k) = ( qi_eta(i,k,j) - f0 ) ** 2
					! QCL
					f0 = mono_cubic_interp( zs_eta(i,k,j), km1, k00, kp1, kp2, &
							zsh(km1),          zsh(k00),          zsh(kp1),          zsh(kp2), &
							qc0(km1)+qr0(km1), qc0(k00)+qr0(k00), qc0(kp1)+qr0(kp1), qc0(kp2)+qr0(kp2) )
					qcl2_eta(i,k) = ( qcl_eta(i,k,j) - f0 ) ** 2
					! QCI
					f0 = mono_cubic_interp( zs_eta(i,k,j), km1, k00, kp1, kp2, &
							zsh(km1), zsh(k00), zsh(kp1), zsh(kp2), &
							qi0(km1)+qp0(km1)-qr0(km1), qi0(k00)+qp0(k00)-qr0(k00), &
							qi0(kp1)+qp0(kp1)-qr0(kp1), qi0(kp2)+qp0(k00)-qr0(kp2) )
					qci2_eta(i,k) = ( qci_eta(i,k,j) - f0 ) ** 2
					! QM=QT+QP
					f0 = mono_cubic_interp( zs_eta(i,k,j), km1, k00, kp1, kp2, &
							zsh(km1),          zsh(k00),          zsh(kp1),          zsh(kp2), &
							qt0(km1)+qp0(km1), qt0(k00)+qp0(k00), qt0(kp1)+qp0(kp1), qt0(kp2)+qp0(kp2) )
					qm2_eta(i,k) = ( qt_eta(i,k,j) + qp_eta(i,k,j) - f0 ) ** 2
				ENDDO
			ENDDO
			
			! Interpolation
			DO k = 1, nzm
				DO i = 1, nx
					! U2
					! Interpolation indices from eta to zs
					km1 = km1u(i,k,j)
					k00 = k00u(i,k,j)
					kp1 = kp1u(i,k,j)
					kp2 = kp2u(i,k,j)
					u2(k) = u2(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zu_eta(i,km1,j), zu_eta(i,k00,j), zu_eta(i,kp1,j), zu_eta(i,kp2,j), &
							u2_eta(i,km1),   u2_eta(i,k00),   u2_eta(i,kp1),   u2_eta(i,kp2) ) &
							* 0.5 * ( mweight(i-1,k,j) + mweight(i,k,j) )
					
					! V2
					km1 = km1v(i,k,j)
					k00 = k00v(i,k,j)
					kp1 = kp1v(i,k,j)
					kp2 = kp2v(i,k,j)
					v2(k) = v2(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zv_eta(i,km1,j), zv_eta(i,k00,j), zv_eta(i,kp1,j), zv_eta(i,kp2,j), &
							v2_eta(i,km1),   v2_eta(i,k00),   v2_eta(i,kp1),   v2_eta(i,kp2) ) &
							* 0.5 * ( mweight(i,k,j-1) + mweight(i,k,j) )
					
					! W2, W3, W4, AUP
					km1 = km1w(i,k,j)
					k00 = k00w(i,k,j)
					kp1 = kp1w(i,k,j)
					kp2 = kp2w(i,k,j)
					w2(k) = w2(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zw_eta(i,km1,j), zw_eta(i,k00,j), zw_eta(i,kp1,j), zw_eta(i,kp2,j), &
							w2_eta(i,km1),   w2_eta(i,k00),   w2_eta(i,kp1),   w2_eta(i,kp2) ) &
							* mweight(i,k,j)
					w3(k) = w3(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zw_eta(i,km1,j), zw_eta(i,k00,j), zw_eta(i,kp1,j), zw_eta(i,kp2,j), &
							w3_eta(i,km1),   w3_eta(i,k00),   w3_eta(i,kp1),   w3_eta(i,kp2) ) &
							* mweight(i,k,j)
					w4(k) = w4(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zw_eta(i,km1,j), zw_eta(i,k00,j), zw_eta(i,kp1,j), zw_eta(i,kp2,j), &
							w4_eta(i,km1),   w4_eta(i,k00),   w4_eta(i,kp1),   w4_eta(i,kp2) ) &
							* mweight(i,k,j)
					IF ( w3d(i,k,j) > 0. ) aup(k) = aup(k) + 1.
					
					! S moments
					km1 = km1s(i,k,j)
					k00 = k00s(i,k,j)
					kp1 = kp1s(i,k,j)
					kp2 = kp2s(i,k,j)
					
					! LWSE2, LWSE3
					liwse2(k) = liwse2(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j),   zs_eta(i,k00,j),   zs_eta(i,kp1,j),   zs_eta(i,kp2,j), &
							liwse2_eta(i,km1), liwse2_eta(i,k00), liwse2_eta(i,kp1), liwse2_eta(i,kp2) ) &
							* mweight(i,k,j)
					liwse3(k) = liwse3(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							liwse3_eta(i,km1), liwse3_eta(i,k00), liwse3_eta(i,kp1), liwse3_eta(i,kp2) ) &
							* mweight(i,k,j)
					
					! MSE2, MSE3
					mse2(k) = mse2(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							mse2_eta(i,km1), mse2_eta(i,k00), mse2_eta(i,kp1), mse2_eta(i,kp2) ) &
							* mweight(i,k,j)
					mse3(k) = mse3(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							mse3_eta(i,km1), mse3_eta(i,k00), mse3_eta(i,kp1), mse3_eta(i,kp2) ) &
							* mweight(i,k,j)
					
					! VDSE2, VDSE3
					vdse2(k) = vdse2(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j),  zs_eta(i,k00,j),  zs_eta(i,kp1,j),  zs_eta(i,kp2,j), &
							vdse2_eta(i,km1), vdse2_eta(i,k00), vdse2_eta(i,kp1), vdse2_eta(i,kp2) ) &
							* mweight(i,k,j)
					vdse3(k) = vdse3(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							vdse3_eta(i,km1), vdse3_eta(i,k00), vdse3_eta(i,kp1), vdse3_eta(i,kp2) ) &
							* mweight(i,k,j)
					
					! QT2, QT3
					qt2(k) = qt2(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							qt2_eta(i,km1),  qt2_eta(i,k00),  qt2_eta(i,kp1),  qt2_eta(i,kp2) ) &
							* mweight(i,k,j)
					qt3(k) = qt3(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							qt3_eta(i,km1),  qt3_eta(i,k00),  qt3_eta(i,kp1),  qt3_eta(i,kp2) ) &
							* mweight(i,k,j)
					
					! QC2
					qc2(k) = qc2(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							qc2_eta(i,km1),  qc2_eta(i,k00),  qc2_eta(i,kp1),  qc2_eta(i,kp2) ) &
							* mweight(i,k,j)
					
					! QI2
					qi2(k) = qi2(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							qi2_eta(i,km1),  qi2_eta(i,k00),  qi2_eta(i,kp1),  qi2_eta(i,kp2) ) &
							* mweight(i,k,j)
					
					! QCL2
					qcl2(k) = qcl2(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							qcl2_eta(i,km1), qcl2_eta(i,k00), qcl2_eta(i,kp1), qcl2_eta(i,kp2) ) &
							* mweight(i,k,j)
					
					! QCI2
					qci2(k) = qci2(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							qci2_eta(i,km1), qci2_eta(i,k00), qci2_eta(i,kp1), qci2_eta(i,kp2) ) &
							* mweight(i,k,j)
					
					! QM2
					qm2(k) = qm2(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							qm2_eta(i,km1),  qm2_eta(i,k00),  qm2_eta(i,kp1),  qm2_eta(i,kp2) ) &
							* mweight(i,k,j)
					
					! SLQT
					slqt(k) = slqt(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							slqt_eta(i,km1), slqt_eta(i,k00), slqt_eta(i,kp1), slqt_eta(i,kp2) ) &
							* mweight(i,k,j)
				ENDDO   ! k
			ENDDO   ! i
		ENDDO   ! j
		! WSKEW, TKE
		DO k = 1, nzm
			wskw(k) = w3(k) / ( w2(k)*factor_xy+1.E-5 )**1.5
			tke(k) = 0.5 * ( u2(k) + v2(k) + w2(k) )
		ENDDO
		
		CALL hbuf_put( 'U2', u2, factor_xy )   ! m2/s2
		CALL hbuf_put( 'V2', v2, factor_xy )   ! m2/s2
		CALL hbuf_put( 'W2', w2, factor_xy )   ! m2/s2
		CALL hbuf_put( 'W3', w3, factor_xy )   ! m3/s3
		CALL hbuf_put( 'W4', w4, factor_xy )   ! m4/s4
		CALL hbuf_put( 'WSKEW', wskw,  factor_xy )
		CALL hbuf_put( 'LIWSE2', liwse2, factor_xy )   ! K2
		CALL hbuf_put( 'LIWSE3', liwse3, factor_xy )   ! K3
		CALL hbuf_put( 'MSE2',  mse2,  factor_xy )   ! K2
		CALL hbuf_put( 'MSE3',  mse3,  factor_xy )   ! K3
		CALL hbuf_put( 'VDSE2', vdse2, factor_xy )   ! K2
		CALL hbuf_put( 'VDSE3', vdse3, factor_xy )   ! K3
		CALL hbuf_put( 'QT2',  qt2,  1.E6*factor_xy )   ! kg2/kg2 ==> g2/kg2
		CALL hbuf_put( 'QT3',  qt3,  1.E9*factor_xy )   ! kg2/kg2 ==> g2/kg2
		CALL hbuf_put( 'QC2',  qc2,  1.E6*factor_xy )   ! kg2/kg2 ==> g2/kg2
		CALL hbuf_put( 'QI2',  qi2,  1.E6*factor_xy )   ! kg2/kg2 ==> g2/kg2
		CALL hbuf_put( 'QCL2', qcl2, 1.E6*factor_xy )   ! kg2/kg2 ==> g2/kg2
		CALL hbuf_put( 'QCI2', qci2, 1.E6*factor_xy )   ! kg2/kg2 ==> g2/kg2
		CALL hbuf_put( 'QM2',  qm2,  1.E6*factor_xy )   ! kg2/kg2 ==> g2/kg2
		CALL hbuf_put( 'LIWSEQT', slqt, factor_xy*fac_cond )   ! K kg/kg ==> K2
		CALL hbuf_put( 'AUP', aup, factor_xy )
		CALL hbuf_put( 'TKE', tke, factor_xy )   ! m2/s2
		
		
		!========== TKE budget ==========
		! Pressure transport + redstiribution is not available yet.
		! See technical note for detail.
		tkeshear(:) = 0.
		tkebuoy(:)  = 0.
		tkeadv(:)   = 0.
		tkepress(:) = 0.
		tkediss(:) = 0.
		tkediff(:) = 0.
		DO j = 1, ny
			jwrf = j0 + j
			! Diagnose tendency
			shear_uweta(:,:) = 0.
			shear_vweta(:,:) = 0.
			adv_ueta(:,:)  = 0.
			adv_veta(:,:)  = 0.
			adv_seta(:,:)  = 0.
			adv_uweta(:,:) = 0.
			adv_vweta(:,:) = 0.
			dif_ueta(:,:) = 0.
			dif_veta(:,:) = 0.
			dif_weta(:,:) = 0.
			prs_ueta(:,:) = 0.
			prs_veta(:,:) = 0.
			prs_weta(:,:) = 0.
			rhop_eta(:,:) = 0.
			DO k = 2, nzm   ! surface value is zero
				DO i = 1, nx
					iwrf = i0 + i
					! Shear production
					! - uwle_eta * du0_eta/dz
					shear_uweta(i,k) = - uwle_eta(iwrf,k,jwrf) * ( u0_eta(i,k,j) - u0_eta(i,k-1,j) ) &
					                 / ( zu_eta(i,k,j) - zu_eta(i,k-1,j) )
					! - vwle_eta * dv0_eta/dz
					shear_vweta(i,k) = - vwle_eta(iwrf,k,jwrf) * ( v0_eta(i,k,j) - v0_eta(i,k-1,j) ) &
					                 / ( zv_eta(i,k,j) - zv_eta(i,k-1,j) )
					
					! Turbulence transport
					! - 1/2 uwle_eta du_eta/dz + uwle_eta * du0_eta/dz
					adv_uweta(i,k) = - uwle_eta(iwrf,k,jwrf) * ( u_eta(i,k,j) - u_eta(i,k-1,j) ) &
					               * 0.5 / ( zu_eta(i,k,j) - zu_eta(i,k-1,j) ) - shear_uweta(i,k)
					! - 1/2 vwle_eta dv_eta/dz + vwle_eta * dv0_eta/dz
					adv_vweta(i,k) = - vwle_eta(iwrf,k,jwrf) * ( v_eta(i,k,j) - v_eta(i,k-1,j) ) &
					               * 0.5 / ( zv_eta(i,k,j) - zv_eta(i,k-1,j) ) - shear_vweta(i,k)
					
					! Diffusion
					! w_eta d[ rhow_eta wwsb_eta ]/dz
					dif_weta(i,k) = w_eta(i,k,j) * ( rhow_eta(i,k,j) * wwsb_eta(iwrf,k,jwrf) &
					              - rhow_eta(i,k-1,j) * wwsb_eta(iwrf,k-1,jwrf) ) &
					              / ( zs_eta(i,k,j) - zs_eta(i,k-1,j) )
					
					! Pressure transport
					! - w_eta d[ps_eta - p0_eta]/dz
					! Unit conversion from hPa to Pa later
					prs_weta(i,k) = - w_eta(i,k,j) * ( ps_eta(i,k,j) - pres0_eta(i,k,j) &
					              - ps_eta(i,k-1,j) + pres0_eta(i,k-1,j) ) &
					              / ( zs_eta(i,k,j) - zs_eta(i,k-1,j) )
				ENDDO
			ENDDO
			DO k = 1, nzm-1   ! top value is zero
				DO i = 1, nx
					iwrf = i0 + i
					! Turbulence transport
					! - ( 1/2 u_eta - u0_eta ) d[ rhow_eta uwle_eta ]/dz
					adv_ueta(i,k) = - ( 0.5 * u_eta(i,k,j) - u0_eta(i,k,j) ) &
					    * 0.5 * ( ( rhow_eta(i-1,k+1,j) + rhow_eta(i,k+1,j) ) * uwle_eta(iwrf,k+1,jwrf) &
					            - ( rhow_eta(i-1,k,j) + rhow_eta(i,k,j) ) * uwle_eta(iwrf,k,jwrf) ) &
					    / ( zuw_eta(i,k+1,j) - zuw_eta(i,k,j) )
					! - ( 1/2 v_eta - v0_eta ) d[ rhow_eta vwle_eta ]/dz
					adv_veta(i,k) = - ( 0.5 * v_eta(i,k,j) - v0_eta(i,k,j) ) &
					    * 0.5 * ( ( rhow_eta(i,k+1,j-1) + rhow_eta(i,k+1,j) ) * vwle_eta(iwrf,k+1,jwrf) &
					            - ( rhow_eta(i,k,j-1) + rhow_eta(i,k,j) ) * vwle_eta(iwrf,k,jwrf) ) &
					    / ( zv_eta(i,k+1,j) - zv_eta(i,k,j) )
					! - 1/2 d[ rhow_eta w w w ]/dz
					adv_seta(i,k) = - 0.5 * ( rhow_eta(i,k+1,j) * ( w_eta(i,k+1,j) ** 3 ) &
					              - rhow_eta(i,k,j) * ( w_eta(i,k,j) ** 3 ) ) &
					              / ( zw_eta(i,k+1,j)-zw_eta(i,k,j) )
					! Diffusion
					! ( u_eta - u0_eta ) d[ rhow_eta uwsb_eta ]/dz
					dif_ueta(i,k) = ( u_eta(i,k,j) - u0_eta(i,k,j) ) &
					    * 0.5 * ( ( rhow_eta(i-1,k+1,j) + rhow_eta(i,k+1,j) ) * uwsb_eta(iwrf,k+1,jwrf) &
					            - ( rhow_eta(i-1,k,j) + rhow_eta(i,k,j) ) * uwsb_eta(iwrf,k,jwrf) ) &
					    / ( zuw_eta(i,k+1,j) - zuw_eta(i,k,j) )
					! ( v_eta - v0_eta ) d[ rhow_eta vwsb_eta ]/dz
					dif_veta(i,k) = ( 0.5 * v_eta(i,k,j) - v0_eta(i,k,j) ) &
					    * 0.5 * ( ( rhow_eta(i,k+1,j-1) + rhow_eta(i,k+1,j) ) * vwsb_eta(iwrf,k+1,jwrf) &
					            - ( rhow_eta(i,k,j-1) + rhow_eta(i,k,j) ) * vwsb_eta(iwrf,k,jwrf) ) &
					    / ( zv_eta(i,k+1,j) - zv_eta(i,k,j) )
				ENDDO
			ENDDO
			! Pressure transport
			! Horizontal gradient is computed with interpolated ps_eta to the U and V levels
			DO i = 1, nx
				khunt1 = 0
				khunt2 = 0
				khunt3 = 0
				khunt4 = 0
				DO k = 1, nzm   ! k loop inside for find_indices
					! - u_eta dps_eta/dx
					! Unit conversion from hPa to Pa later
					CALL find_indices(zu_eta(i,k,j), zs_eta(i-1,1:nzm,j), nzm, km1, k00, kp1, kp2,khunt1)
					ppm1 = mono_cubic_interp( zu_eta(i,k,j), km1, k00, kp1, kp2, &
							zs_eta(i-1,km1,j), zs_eta(i-1,k00,j), zs_eta(i-1,kp1,j), zs_eta(i-1,kp2,j), &
							ps_eta(i-1,km1,j), ps_eta(i-1,k00,j), ps_eta(i-1,kp1,j), ps_eta(i-1,kp2,j) )
					CALL find_indices(zu_eta(i,k,j), zs_eta(i,1:nzm,j), nzm, km1, k00, kp1, kp2, khunt2)
					pp00 = mono_cubic_interp( zu_eta(i,k,j), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							ps_eta(i,km1,j), ps_eta(i,k00,j), ps_eta(i,kp1,j), ps_eta(i,kp2,j) )
					prs_ueta(i,k) = - u_eta(i,k,j) * ( pp00 - ppm1 ) / dx
					
					! - v_eta dps_eta/dy
					! Unit conversion from hPa to Pa later
					CALL find_indices(zv_eta(i,k,j), zs_eta(i,1:nzm,j-1), nzm, km1, k00, kp1, kp2,khunt3)
					ppm1 = mono_cubic_interp( zv_eta(i,k,j), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j-1), zs_eta(i,k00,j-1), zs_eta(i,kp1,j-1), zs_eta(i,kp2,j-1), &
							ps_eta(i,km1,j-1), ps_eta(i,k00,j-1), ps_eta(i,kp1,j-1), ps_eta(i,kp2,j-1) )
					CALL find_indices(zv_eta(i,k,j), zs_eta(i,1:nzm,j), nzm, km1, k00, kp1, kp2, khunt4)
					pp00 = mono_cubic_interp( zv_eta(i,k,j), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							ps_eta(i,km1,j), ps_eta(i,k00,j), ps_eta(i,kp1,j), ps_eta(i,kp2,j) )
					prs_veta(i,k) = - v_eta(i,k,j) * ( pp00 - ppm1 ) / dy
				ENDDO
			ENDDO
			! Buoyant production
			! RHO flux (m2/s3) is the buoyancy flux for full compressible system.
			! Perturbation of RHO has to be computed here.
			DO k = 1, nzm
				DO i = 1, nx
					km1 = km1s02s(i,k,j)
					k00 = k00s02s(i,k,j)
					kp1 = kp1s02s(i,k,j)
					kp2 = kp2s02s(i,k,j)
					f0 = mono_cubic_interp( zs_eta(i,k,j), km1, k00, kp1, kp2, &
							zsh(km1),  zsh(k00),  zsh(kp1),  zsh(kp2), &
							rho0(km1), rho0(k00), rho0(kp1), rho0(kp2) )
					rhop_eta(i,k) = rhom_eta(i,k,j) - f0
				ENDDO
			ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! buoy_eta = - g mu' / (mu' + mub), !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! where mu = p_dry@surface - p_dry@top, mub = reference state of mu (hydrostatic balance), !!!!
!!!!! has to be converted to -g rho' / [rho], where rho is moist air. At this point, it is !!!!!!!!
!!!!! unclear how to do this conversion. Thus, the following lines are uncommented and !!!!!!!!!!!!
!!!!! buoyant production is diagnosed with rho and w below. Also, the inserted lines in !!!!!!!!!!! !!!!! module_big_step_utilities_em.F to fill buoy_eta are uncommented. !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!			DO k = 1, nzm
!!				DO i = 1, nx
!!					iwrf = i0 + i
!!					BUOYA = w_eta(i,k,j) * buoy_eta(iwrf,k,jwrf) <=== W level
!!				ENDDO
!!			ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			
			
			! Interpolation
			DO k = 1, nzm
				DO i = 1, nx
					iwrf = i0 + i
					
					! @ U
					km1 = km1u(i,k,j)
					k00 = k00u(i,k,j)
					kp1 = kp1u(i,k,j)
					kp2 = kp2u(i,k,j)
					! ADVTR
					tkeadv(k) = tkeadv(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zu_eta(i,km1,j), zu_eta(i,k00,j), zu_eta(i,kp1,j), zu_eta(i,kp2,j), &
							adv_ueta(i,km1), adv_ueta(i,k00), adv_ueta(i,kp1), adv_ueta(i,kp2) ) / rho(k)
					! DIFTR
					tkediff(k) = tkediff(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zu_eta(i,km1,j), zu_eta(i,k00,j), zu_eta(i,kp1,j), zu_eta(i,kp2,j), &
							dif_ueta(i,km1), dif_ueta(i,k00), dif_ueta(i,kp1), dif_ueta(i,kp2) ) / rho(k)
					! PRESSTR
					tkepress(k) = tkepress(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zu_eta(i,km1,j), zu_eta(i,k00,j), zu_eta(i,kp1,j), zu_eta(i,kp2,j), &
							prs_ueta(i,km1), prs_ueta(i,k00), prs_ueta(i,kp1), prs_ueta(i,kp2) )
							
					! @ V
					km1 = km1v(i,k,j)
					k00 = k00v(i,k,j)
					kp1 = kp1v(i,k,j)
					kp2 = kp2v(i,k,j)
					! ADVTR
					tkeadv(k) = tkeadv(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zv_eta(i,km1,j), zv_eta(i,k00,j), zv_eta(i,kp1,j), zv_eta(i,kp2,j), &
							adv_veta(i,km1), adv_veta(i,k00), adv_veta(i,kp1), adv_veta(i,kp2) ) / rho(k)
					! DIFTR
					tkediff(k) = tkediff(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zv_eta(i,km1,j), zv_eta(i,k00,j), zv_eta(i,kp1,j), zv_eta(i,kp2,j), &
							dif_veta(i,km1), dif_veta(i,k00), dif_veta(i,kp1), dif_veta(i,kp2) ) / rho(k)
					! PRESSTR
					tkepress(k) = tkepress(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zv_eta(i,km1,j), zv_eta(i,k00,j), zv_eta(i,kp1,j), zv_eta(i,kp2,j), &
							prs_veta(i,km1), prs_veta(i,k00), prs_veta(i,kp1), prs_veta(i,kp2) )
					
					! @ W
					km1 = km1w(i,k,j)
					k00 = k00w(i,k,j)
					kp1 = kp1w(i,k,j)
					kp2 = kp2w(i,k,j)
					! DIFTR
					tkediff(k) = tkediff(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zw_eta(i,km1,j), zw_eta(i,k00,j), zw_eta(i,kp1,j), zw_eta(i,kp2,j), &
							dif_weta(i,km1), dif_weta(i,k00), dif_weta(i,kp1), dif_weta(i,kp2) ) / rho(k)
					! PRESSTR
					tkepress(k) = tkepress(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zw_eta(i,km1,j), zw_eta(i,k00,j), zw_eta(i,kp1,j), zw_eta(i,kp2,j), &
							prs_weta(i,km1), prs_weta(i,k00), prs_weta(i,kp1), prs_weta(i,kp2) )
					
					! @ UW
					km1 = km1uw(i,k,j)
					k00 = k00uw(i,k,j)
					kp1 = kp1uw(i,k,j)
					kp2 = kp2uw(i,k,j)
					! SHEAR
					tkeshear(k) = tkeshear(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zuw_eta(i,km1,j),  zuw_eta(i,k00,j),  zuw_eta(i,kp1,j),  zuw_eta(i,kp2,j), &
							shear_uweta(i,km1),shear_uweta(i,k00),shear_uweta(i,kp1),shear_uweta(i,kp2) ) &
							* 0.5 * ( mweight(i-1,k,j) + mweight(i,k,j) )
					! ADVTR
					tkeadv(k) = tkeadv(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zuw_eta(i,km1,j), zuw_eta(i,k00,j), zuw_eta(i,kp1,j), zuw_eta(i,kp2,j), &
							adv_uweta(i,km1), adv_uweta(i,k00), adv_uweta(i,kp1), adv_uweta(i,kp2) ) &
							* 0.5 * ( mweight(i-1,k,j) + mweight(i,k,j) )
					
					! @ VW
					km1 = km1vw(i,k,j)
					k00 = k00vw(i,k,j)
					kp1 = kp1vw(i,k,j)
					kp2 = kp2vw(i,k,j)
					! SHEAR
					tkeshear(k) = tkeshear(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zvw_eta(i,km1,j), zvw_eta(i,k00,j), zvw_eta(i,kp1,j), zvw_eta(i,kp2,j), &
							shear_vweta(i,km1),shear_vweta(i,k00),shear_vweta(i,kp1),shear_vweta(i,kp2) ) &
							* 0.5 * ( mweight(i,k,j-1) + mweight(i,k,j) )
					! ADVTR
					tkeadv(k) = tkeadv(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zvw_eta(i,km1,j), zvw_eta(i,k00,j), zvw_eta(i,kp1,j), zvw_eta(i,kp2,j), &
							adv_vweta(i,km1), adv_vweta(i,k00), adv_vweta(i,kp1), adv_vweta(i,kp2) ) &
							* 0.5 * ( mweight(i,k,j-1) + mweight(i,k,j) )
					
					! @ S
					km1 = km1s(i,k,j)
					k00 = k00s(i,k,j)
					kp1 = kp1s(i,k,j)
					kp2 = kp2s(i,k,j)
					! ADVTR
					tkeadv(k) = tkeadv(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							adv_seta(i,km1), adv_seta(i,k00), adv_seta(i,kp1), adv_seta(i,kp2) ) / rho(k)
					! BUOYA
					! RHO flux (m2/s3) is the buoyancy flux for full compressible system.
					fprime = mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							rhop_eta(i,km1), rhop_eta(i,k00), rhop_eta(i,kp1), rhop_eta(i,kp2) )
					tkebuoy(k) = tkebuoy(k) - ggr * w3d(i,k,j) * fprime / rho(k)
					! DIFTR, DISSIP
					tmp = mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							e_shear_eta(iwrf,km1,jwrf), e_shear_eta(iwrf,k00,jwrf), &
							e_shear_eta(iwrf,kp1,jwrf), e_shear_eta(iwrf,kp2,jwrf) ) * mweight(i,k,j)
					tkediff(k) = tkediff(k) + tmp
					tkediss(k) = tkediss(k) - tmp
				ENDDO ! i
			ENDDO ! k
		ENDDO ! j
		
		CALL hbuf_put( 'SHEAR',   tkeshear, factor_xy )   ! m2/s3
		CALL hbuf_put( 'BUOYA',   tkebuoy,  factor_xy )   ! m2/s3
		CALL hbuf_put( 'ADVTR',   tkeadv,   factor_xy )   ! m2/s3
		CALL hbuf_put( 'PRESSTR', tkepress, factor_xy*100. )   ! m2/s3, 100 for hPa => Pa (kg/m/s2)
		CALL hbuf_put( 'DISSIP',  tkediss,  factor_xy )   ! m2/s3
		CALL hbuf_put( 'DIFTR',   tkediff,  factor_xy )   ! m2/s3
		
		
		!========== SGS TKE budget ==========
		! Advective tendency due to resolve scale + SGS wind = eadv_eta + edif_eta
		! eadv_eta, edif_eta, e_shear_eta, e_buoy_eta, e_diss_eta are defined as a WRF 3D array
		tkesb(:) = 0.
		tkesbadv(:) = 0.
		tkesbshear(:) = 0.
		tkesbbuoy(:) = 0.
		tkesbdiss(:) = 0.
		tkm(:) = 0.
		tkh(:) = 0.
		DO j = 1, ny
			jwrf = j0 + j
			
			! Diagnose transport
			DO k = 1, nzm-1
				DO i = 1, nx
					iwrf = i0 + i
					adv_seta(i,k) &
					   = - ( rhow_eta(i,k+1,j) * ( ewle_eta(iwrf,k+1,jwrf) + ewsb_eta(iwrf,k+1,jwrf) ) &
					       - rhow_eta(i,k,j) * ( ewle_eta(iwrf,k,jwrf) + ewsb_eta(iwrf,k,jwrf) ) ) &
					   / ( zw_eta(i,k+1,j) - zw_eta(i,k,j) )
				ENDDO
			ENDDO
			
			! Interpolation
			DO k = 1, nzm
				DO i = 1, nx
					iwrf = i0 + i
					
					km1 = km1s(i,k,j)
					k00 = k00s(i,k,j)
					kp1 = kp1s(i,k,j)
					kp2 = kp2s(i,k,j)
					
					tkesb(k) = tkesb(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j),  zs_eta(i,k00,j),  zs_eta(i,kp1,j),  zs_eta(i,kp2,j), &
							tke_eta(i,km1,j), tke_eta(i,k00,j), tke_eta(i,kp1,j), tke_eta(i,kp2,j) ) &
							* mweight(i,k,j)
					
					tkesbadv(k) = tkesbadv(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							adv_seta(i,km1), adv_seta(i,k00), adv_seta(i,kp1), adv_seta(i,kp2) ) / rho(k)
					
					tkesbshear(k) = tkesbshear(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							e_shear_eta(iwrf,km1,jwrf), e_shear_eta(iwrf,k00,jwrf), &
							e_shear_eta(iwrf,kp1,jwrf), e_shear_eta(iwrf,kp2,jwrf) ) * mweight(i,k,j)
					
					tkesbbuoy(k) = tkesbbuoy(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							e_buoy_eta(iwrf,km1,jwrf), e_buoy_eta(iwrf,k00,jwrf), &
							e_buoy_eta(iwrf,kp1,jwrf), e_buoy_eta(iwrf,kp2,jwrf) ) * mweight(i,k,j)
					
					tkesbdiss(k) = tkesbdiss(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							e_diss_eta(iwrf,km1,jwrf), e_diss_eta(iwrf,k00,jwrf), &
							e_diss_eta(iwrf,kp1,jwrf), e_diss_eta(iwrf,kp2,jwrf) ) * mweight(i,k,j)
					
					tkm(k) = tkm(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j),  zs_eta(i,k00,j),  zs_eta(i,kp1,j),  zs_eta(i,kp2,j), &
							tkm_eta(i,km1,j), tkm_eta(i,k00,j), tkm_eta(i,kp1,j), tkm_eta(i,kp2,j) ) &
							* mweight(i,k,j)
					
					tkh(k) = tkh(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j),  zs_eta(i,k00,j),  zs_eta(i,kp1,j),  zs_eta(i,kp2,j), &
							tkh_eta(i,km1,j), tkh_eta(i,k00,j), tkh_eta(i,kp1,j), tkh_eta(i,kp2,j) ) &
							* mweight(i,k,j)
				ENDDO
			ENDDO
		ENDDO
		
		CALL hbuf_put( 'TKES', tkesb, factor_xy )   ! m2/s2
		CALL hbuf_put( 'SHEARS',  tkesbshear, factor_xy )   ! m2/s3
		CALL hbuf_put( 'BUOYAS',  tkesbbuoy,  factor_xy )   ! m2/s3
		CALL hbuf_put( 'ADVTRS',  tkesbadv,   factor_xy )   ! m2/s3
		CALL hbuf_put( 'DISSIPS', tkesbdiss,  factor_xy )   ! m2/s3
		CALL hbuf_put( 'TKM', tkm, factor_xy )   ! m2/s
		CALL hbuf_put( 'TKH', tkh, factor_xy )   ! m2/s
		
		
		!========== Radiation ==========
		! Radiative flux (radlwup_eta, radlwdn_eta, radswup_eta, radswdn_eta ) have to be collected
		! in each radiation code. These arrays are defined as a WRF 3D array
		! Implemented: module_ra_dyclw.F
		radlwup(:) = 0.
		radlwdn(:) = 0.
		radswup(:) = 0.
		radswdn(:) = 0.
		radqrlw(:) = 0.
		radqrsw(:) = 0.
		radqr(:)  = 0.
		radqrs(:) = 0.
		radqrc(:) = 0.
		DO j = 1, ny
			jwrf = j0 + j
			DO k = 1, nzm
				DO i = 1, nx
					iwrf = i0 + i
					
					! Flux (W/m2) <== mass is included
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Should I divid flux by rhow_eta, then mweight after interpolation? !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					km1 = km1w(i,k,j)
					k00 = k00w(i,k,j)
					kp1 = kp1w(i,k,j)
					kp2 = kp2w(i,k,j)
					
					radlwup(k) = radlwup(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zw_eta(i,km1,j), zw_eta(i,k00,j), zw_eta(i,kp1,j), zw_eta(i,kp2,j), &
							radlwup_eta(iwrf,km1,jwrf), radlwup_eta(iwrf,k00,jwrf), &
							radlwup_eta(iwrf,kp1,jwrf), radlwup_eta(iwrf,kp2,jwrf) )
					
					radlwdn(k) = radlwdn(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zw_eta(i,km1,j), zw_eta(i,k00,j), zw_eta(i,kp1,j), zw_eta(i,kp2,j), &
							radlwdn_eta(iwrf,km1,jwrf), radlwdn_eta(iwrf,k00,jwrf), &
							radlwdn_eta(iwrf,kp1,jwrf), radlwdn_eta(iwrf,kp2,jwrf) )
					
					radswup(k) = radswup(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zw_eta(i,km1,j), zw_eta(i,k00,j), zw_eta(i,kp1,j), zw_eta(i,kp2,j), &
							radswup_eta(iwrf,km1,jwrf), radswup_eta(iwrf,k00,jwrf), &
							radswup_eta(iwrf,kp1,jwrf), radswup_eta(iwrf,kp2,jwrf) )
					
					radswdn(k) = radswdn(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zw_eta(i,km1,j), zw_eta(i,k00,j), zw_eta(i,kp1,j), zw_eta(i,kp2,j), &
							radswdn_eta(iwrf,km1,jwrf), radswdn_eta(iwrf,k00,jwrf), &
							radswdn_eta(iwrf,kp1,jwrf), radswdn_eta(iwrf,kp2,jwrf) )
					
					! Heating (K/s)
					km1 = km1s(i,k,j)
					k00 = k00s(i,k,j)
					kp1 = kp1s(i,k,j)
					kp2 = kp2s(i,k,j)
					
					radqrlw(k) = radqrlw(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							radqrlw_eta(i,km1,j), radqrlw_eta(i,k00,j), &
							radqrlw_eta(i,kp1,j), radqrlw_eta(i,kp2,j) ) * mweight(i,k,j)
					
					radqrsw(k) = radqrsw(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							radqrsw_eta(i,km1,j), radqrsw_eta(i,k00,j), &
							radqrsw_eta(i,kp1,j), radqrsw_eta(i,kp2,j) ) * mweight(i,k,j)
					
					radqr3d(i,k,j) = mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							radqrlw_eta(i,km1,j)+radqrsw_eta(i,km1,j), &
							radqrlw_eta(i,k00,j)+radqrsw_eta(i,k00,j), &
							radqrlw_eta(i,kp1,j)+radqrsw_eta(i,kp1,j), &
							radqrlw_eta(i,kp2,j)+radqrsw_eta(i,kp2,j) )
					radqr(k) = radqr(k) + radqr3d(i,k,j) * mweight(i,k,j)
				ENDDO
			ENDDO
		ENDDO
		num = 0.
		DO j = 1, ny
			DO i = 1, nx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Clear and cloudy sky reating rates are not mass-weighted. Leave now for future version. !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				IF ( cwp(i,j) > cwpmax ) THEN
					num = num + 1.
					radqrc(:) = radqrc(:) + radqr3d(i,:,j)
				ELSE
					radqrs(:) = radqrs(:) + radqr3d(i,:,j)
				ENDIF
			ENDDO
		ENDDO
		
		CALL hbuf_put( 'RADLWUP', radlwup, factor_xy )   ! W/m2
		CALL hbuf_put( 'RADLWDN', radlwdn, factor_xy )   ! W/m2
		CALL hbuf_put( 'RADSWUP', radswup, factor_xy )   ! W/m2
		CALL hbuf_put( 'RADSWDN', radswdn, factor_xy )   ! W/m2
		CALL hbuf_put( 'RADQRLW', radqrlw, 3600.*factor_xy )   ! K/s ==> K/hour
		CALL hbuf_put( 'RADQRSW', radqrsw, 3600.*factor_xy )   ! K/s ==> K/hour
		CALL hbuf_put( 'RADQR',   radqr,   3600.*factor_xy )   ! K/s ==> K/hour
		CALL hbuf_put( 'RADQRS',  radqrs, 3600./(REAL(nx*ny)-num+1.E-5) )   ! clear sky column mean
		CALL hbuf_put( 'RADQRC',  radqrc, 3600./(num+1.E-5) )   ! cloudy sky column mean
		
		
		!========== Output from module_wrfles.F ==========
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Output from module_wrfles.F is not mass-weighted. Leave now for future version. !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Observed U, V, W, TP, QV
		CALL hbuf_put( 'UOBS', uobs, 1. )   ! m/s
		CALL hbuf_put( 'VOBS', vobs, 1. )   ! m/s
		CALL hbuf_put( 'WOBS', wobs, 1. )   ! m/s
		CALL hbuf_put( 'TPOBS', tpobs, 1. )   ! K for TP
		CALL hbuf_put( 'QVOBS', qvobs, 1.E3 )   ! kg/kg ==> g/kg
		
		! Large scale horizontal advective tendency
		CALL hbuf_put( 'TPLSADVH', tplshadv, 3600. )   ! K/s ==> K/hour
		CALL hbuf_put( 'QVLSADVH', qvlshadv, 3600.*1.E3 )   ! kg/kg/s ==> g/kg/hour
		
		! Large scale vertical advective tendency = subsidence
		ulsvadv(:) = 0.
		vlsvadv(:) = 0.
		tplsvadv(:) = 0.
		qvlsvadv(:) = 0.
		DO j = 1, ny
			jwrf = j0 + j
			DO k = 1, nzm
				DO i = 1, nx
					iwrf = i0 + i
					! U
					km1 = km1u(i,k,j)
					k00 = k00u(i,k,j)
					kp1 = kp1u(i,k,j)
					kp2 = kp2u(i,k,j)
					ulsvadv(k) = ulsvadv(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zu_eta(i,km1,j), zu_eta(i,k00,j), zu_eta(i,kp1,j), zu_eta(i,kp2,j), &
							ulsvadv_eta(iwrf,km1,jwrf), ulsvadv_eta(iwrf,k00,jwrf), &
							ulsvadv_eta(iwrf,kp1,jwrf), ulsvadv_eta(iwrf,kp2,jwrf) ) &
							* 0.5 * ( mweight(i-1,k,j) + mweight(i,k,j) )
					! V
					km1 = km1v(i,k,j)
					k00 = k00v(i,k,j)
					kp1 = kp1v(i,k,j)
					kp2 = kp2v(i,k,j)
					vlsvadv(k) = vlsvadv(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zv_eta(i,km1,j), zv_eta(i,k00,j), zv_eta(i,kp1,j), zv_eta(i,kp2,j), &
							vlsvadv_eta(iwrf,km1,jwrf), vlsvadv_eta(iwrf,k00,jwrf), &
							vlsvadv_eta(iwrf,kp1,jwrf), vlsvadv_eta(iwrf,kp2,jwrf) ) &
							* 0.5 * ( mweight(i,k,j-1) + mweight(i,k,j) )
					!S
					km1 = km1s(i,k,j)
					k00 = k00s(i,k,j)
					kp1 = kp1s(i,k,j)
					kp2 = kp2s(i,k,j)
					tplsvadv(k) = tplsvadv(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							tplsvadv_eta(iwrf,km1,jwrf), tplsvadv_eta(iwrf,k00,jwrf), &
							tplsvadv_eta(iwrf,kp1,jwrf), tplsvadv_eta(iwrf,kp2,jwrf) ) * mweight(i,k,j)
					qvlsvadv(k) = qvlsvadv(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
							zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
							qvlsvadv_eta(iwrf,km1,jwrf), qvlsvadv_eta(iwrf,k00,jwrf), &
							qvlsvadv_eta(iwrf,kp1,jwrf), qvlsvadv_eta(iwrf,kp2,jwrf) ) * mweight(i,k,j)
				ENDDO
			ENDDO
		ENDDO
		CALL hbuf_put( 'ULSADVV', ulsvadv, 3600.*factor_xy )   ! m/s2 ==> m/s/hour
		CALL hbuf_put( 'VLSADVV', vlsvadv, 3600.*factor_xy )   ! m/s2 ==> m/s/hour
		CALL hbuf_put( 'TPLSADVV', tplsvadv, 3600.*factor_xy )   ! K/s ==> K/hour
		CALL hbuf_put( 'QVLSADVV', qvlsvadv, 3600.*1.E3*factor_xy )   ! kg/kg/s ==> g/kg/hour
		
		! Nudging
		CALL hbuf_put( 'UNUDGE', unudge, 3600. )   ! m/s2 ==> m/s/hour
		CALL hbuf_put( 'VNUDGE', vnudge, 3600. )   ! m/s2 ==> m/s/hour
		CALL hbuf_put( 'TPNUDGE', tpnudge, 3600. )   ! K/s ==> K/hour
		CALL hbuf_put( 'QVNUDGE', qvnudge, 3600.*1.E3 )   ! kg/kg/s ==> g/kg/hour
		CALL hbuf_put( 'QNTNUDGE', qntnudge, 3600. )   ! #/mg/s ==> #/mg/hour
		
		
		!========== Tendencies for LWSE and QT ========== Not implemented yet
		!========== Momentum flux balance ========== Not implemented yet
		!========== LWSE2 and QT2 budget ========== Not implemented yet
		!========== TLW and QTW budget ========== Not implemented yet
		!========== Apparent heat/moisture sources/sinks ========== Not implemented yet
		!========== ISCCP simulater ========== Not implemented yet !CALL isccp_get()
		!========== Mass flux, hydrometeor fraction statistics ========== Not implemented yet
		!========== Updraft Core statistics ==========
		!========== Cloud Downdraft Core statistics ==========
		!========== Convective velocity scale, wstar3 ========== Not implemented yet
		!========== UW ADDITIONS ========== Not implemented yet
		!===== Average surface wind stress used for conditional statistics ===== Not implemented yet
		
		
		!========== Conditional statistics ==========
		! Initialize mask array, which will be one only where the conditional is satisfied, and
		! zero everywhere else.
		!
		! Tak: Conditional avarage will be summed only if condavg_mask = 1 for k00 and kp1
		condavg_mask(:,:,:,:) = 0.
		DO j = 1, ny
			DO k = 1, nzm
				kp1 = MIN( k+1, nzm ) ! use nzm if k+1 > nzm
				DO i = 1, nx
					IF (LES) THEN
						coef = 0.
					ELSE
						coef = MIN( 1.E-5, 0.01*qsatw(tabs_eta(i,k,j),ps_eta(i,k,j)) )
					ENDIF
					IF ( icondavg_cld > 0 .AND. qn_eta(i,k,j) > coef ) THEN
						condavg_mask(i,k,j,icondavg_cld) = 1. ! cloud
					ENDIF
					IF ( icondavg_cldl > 0 .AND. qc_eta(i,k,j) > coef ) THEN
						condavg_mask(i,k,j,icondavg_cldl) = 1. ! liquid cloud
					ENDIF
					IF ( icondavg_cldi > 0 .AND. qi_eta(i,k,j) > coef ) THEN
						condavg_mask(i,k,j,icondavg_cldi) = 1. ! ice cloud
					ENDIF
					IF ( icondavg_cndl > 0 .AND. qcl_eta(i,k,j) > coef ) THEN
						condavg_mask(i,k,j,icondavg_cndl) = 1. ! liquid water
					ENDIF
					IF ( icondavg_cndi > 0 .AND. qci_eta(i,k,j) > coef ) THEN
						condavg_mask(i,k,j,icondavg_cndi) = 1. ! ice water
					ENDIF
					IF ( icondavg_cor > 0 ) THEN
						! updraft (w>1) core (tv'>0) statistics in LES, buoyant cloudy statistics
						condition_cl = qn_eta(i,k,j) > coef
						condition = tpv_eta(i,k,j) > tpv0_eta(i,k,j)
						IF (CSRM) condition = condition .AND. w_eta(i,k,j)+w_eta(i,kp1,j) > 2.
						IF (LES) condition = condition_cl .AND. condition
						IF (condition) condavg_mask(i,k,j,icondavg_cor) = 1. ! core
					ENDIF
					IF ( icondavg_cordn > 0 ) THEN
						! downdraft (w<-1) core (tv'>0) statistics in LES, buoyant, saturated or rainy
						! statistics
						condition_cl = qn_eta(i,k,j) > coef .OR. qp_eta(i,k,j) > 1.E-4
						condition = tpv_eta(i,k,j) < tpv0_eta(i,k,j)
						IF (CSRM) condition = condition .AND. w_eta(i,k,j)+w_eta(i,kp1,j) < -2.
						IF (LES) condition = condition_cl .AND. condition
						IF (condition) condavg_mask(i,k,j,icondavg_cordn) = 1. ! downdraft core
					ENDIF
					condition_cl = qn_eta(i,k,j) > coef
					condition = w_eta(i,k,j)+w_eta(i,kp1,j) >= 0.
					IF ( icondavg_satup > 0 .AND. (condition_cl.AND.condition) ) THEN
						condavg_mask(i,k,j,icondavg_satup) = 1. ! saturated updraft
					ENDIF
					IF ( icondavg_satdn > 0 .AND. (condition_cl.AND.(.NOT.condition)) ) THEN
						condavg_mask(i,k,j,icondavg_satdn) = 1. ! saturated downdraft
					ENDIF
					IF( icondavg_env > 0 .AND. (.NOT.condition_cl) ) THEN
						condavg_mask(i,k,j,icondavg_env) = 1. ! cloud-free environment
					ENDIF
					IF( icondavg_wup > 0 .AND. condition ) THEN
						condavg_mask(i,k,j,icondavg_wup) = 1. ! updraft
					ENDIF
					IF( icondavg_wdn > 0 .AND. (.NOT.condition) ) THEN
						condavg_mask(i,k,j,icondavg_wdn) = 1. ! downdraft
					ENDIF
				ENDDO
			ENDDO
		ENDDO
		
		! Compute RHO for mass-weighted conditional average
		
		
		!
		DO ncond = 1, ncondavg
			cld(:) = 0.
			rhocl(:) = 0.
			ucl(:) = 0.
			vcl(:) = 0.
			wcl(:) = 0.
			dsecl(:)  = 0.
			liwsecl(:) = 0.
			msecl(:)  = 0.
			tpvcl(:)  = 0.
			tpvcla(:) = 0.
			tabscl(:) = 0.
			qtcl(:) = 0.
			qccl(:) = 0.
			qicl(:) = 0.
			qncl(:) = 0.
			qpcl(:) = 0.
			etacl(:) = 0.
			ssrcl(:) = 0.
			qnccl(:) = 0.
			qnicl(:) = 0.
			qnrcl(:) = 0.
			qnscl(:) = 0.
			qngcl(:) = 0.
			qnncl(:) = 0.
			wacl(:) = 0.
			tvwcl(:) = 0.
			twcl(:)  = 0.
			qwcl(:)  = 0.
			qcwcl(:) = 0.
			qiwcl(:) = 0.
			DO j = 1, ny
				DO k = 1, nzm
					DO i = 1, nx
						km1 = km1s(i,k,j)
						k00 = k00s(i,k,j)
						kp1 = kp1s(i,k,j)
						kp2 = kp2s(i,k,j)
						IF (condavg_mask(i,k00,j,ncond) > 0 .AND. condavg_mask(i,kp1,j,ncond) > 0) THEN
							cld(k) = cld(k) + 1.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Conditional average is not mass-weighted. This requires conditional mean mass before !!!!!!!!
!!!!! summation. Leave now for future version. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
							! U, V, W
							ucl(k) = ucl(k) + u3d(i,k,j)
							vcl(k) = vcl(k) + v3d(i,k,j)
							wcl(k) = wcl(k) + w3d(i,k,j)
							
							! Scalar
							rhocl(k)  = rhocl(k) + rho3d(i,k,j)
							dsecl(k)  = dsecl(k) + dse3d(i,k,j)
							liwsecl(k) = liwsecl(k) + liwse3d(i,k,j)
							msecl(k)  = msecl(k) + mse3d(i,k,j)
							tpvcl(k)  = tpvcl(k) + tpv3d(i,k,j)
							tabscl(k) = tabscl(k) + tabs3d(i,k,j)
							qtcl(k) = qtcl(k) + qt3d(i,k,j)
							qccl(k) = qccl(k) + qc3d(i,k,j)
							qicl(k) = qicl(k) + qi3d(i,k,j)
							qncl(k) = qncl(k) + qn3d(i,k,j)
							qpcl(k) = qpcl(k) + qp3d(i,k,j)
							etacl(k) = etacl(k) + eta3d(i,k,j)
							ssrcl(k) = ssrcl(k) + ssr3d(i,k,j)
							qnccl(k) = qnccl(k) + qnc3d(i,k,j)
							qnicl(k) = qnicl(k) + qni3d(i,k,j)
							qnrcl(k) = qnrcl(k) + qnr3d(i,k,j)
							qnscl(k) = qnscl(k) + qns3d(i,k,j)
							qngcl(k) = qngcl(k) + qng3d(i,k,j)
							qnncl(k) = qnncl(k) + qnn3d(i,k,j)
							
							! Anomaly
							tpvcla(k) = tpvcla(k) + mono_cubic_interp( zs(k), km1, k00, kp1, kp2, &
									zs_eta(i,km1,j), zs_eta(i,k00,j), zs_eta(i,kp1,j), zs_eta(i,kp2,j), &
									tpv_eta(i,km1,j)-tpv0_eta(i,km1,j), tpv_eta(i,k00,j)-tpv0_eta(i,k00,j), &
									tpv_eta(i,kp1,j)-tpv0_eta(i,kp1,j), tpv_eta(i,kp2,j)-tpv0_eta(i,kp2,j) )
							
							! Turbulence
							twcl(k)  = twcl(k)  + w3d(i,k,j) * liwse3d(i,k,j)
							tvwcl(k) = tvwcl(k) + w3d(i,k,j) * tpv3d(i,k,j)
							qwcl(k)  = qwcl(k)  + w3d(i,k,j) * qt3d(i,k,j)
							qcwcl(k) = qcwcl(k) + w3d(i,k,j) * qc3d(i,k,j)
							qiwcl(k) = qiwcl(k) + w3d(i,k,j) * qi3d(i,k,j)
						ENDIF
					ENDDO
				ENDDO
			ENDDO
			DO k = 1, nzm
				condavg_factor(k,ncond) = condavg_factor(k,ncond) + cld(k)
				wacl(k) = wcl(k)
			ENDDO
			
			CALL hbuf_put( TRIM(condavgname(ncond)), cld, factor_xy )
			CALL hbuf_put( 'RHO'//TRIM(condavgname(ncond)), rhocl, 1. )
			CALL hbuf_put( 'U'//TRIM(condavgname(ncond)), ucl, 1. )
			CALL hbuf_put( 'V'//TRIM(condavgname(ncond)), vcl, 1. )
			CALL hbuf_put( 'W'//TRIM(condavgname(ncond)), wcl, 1. )
			CALL hbuf_put( 'DSE'//TRIM(condavgname(ncond)), dsecl, 1. )
			CALL hbuf_put( 'LIWSE'//TRIM(condavgname(ncond)), liwsecl, 1. )
			CALL hbuf_put( 'MSE'//TRIM(condavgname(ncond)), msecl, 1. )
			CALL hbuf_put( 'TABS'//TRIM(condavgname(ncond)),tabscl, 1. )
			CALL hbuf_put( 'TPV'//TRIM(condavgname(ncond)), tpvcl, 1. )
			CALL hbuf_put( 'TPV'//TRIM(condavgname(ncond))//'A', tpvcla, 1. )
			CALL hbuf_put( 'QT'//TRIM(condavgname(ncond)), qtcl, 1.E3 )
			CALL hbuf_put( 'QC'//TRIM(condavgname(ncond)), qccl, 1.E3 )
			CALL hbuf_put( 'QI'//TRIM(condavgname(ncond)), qicl, 1.E3 )
			CALL hbuf_put( 'QN'//TRIM(condavgname(ncond)), qncl, 1.E3 )
			CALL hbuf_put( 'QP'//TRIM(condavgname(ncond)), qpcl, 1.E3 )
			CALL hbuf_put( 'ETA'//TRIM(condavgname(ncond)), etacl, 1.E3 )
			CALL hbuf_put( 'SSR'//TRIM(condavgname(ncond)), ssrcl, 1. )
			CALL hbuf_put( 'QNC'//TRIM(condavgname(ncond)), qnccl, 1. )
			CALL hbuf_put( 'QNI'//TRIM(condavgname(ncond)), qnicl, 1. )
			CALL hbuf_put( 'QNR'//TRIM(condavgname(ncond)), qnrcl, 1. )
			CALL hbuf_put( 'QNS'//TRIM(condavgname(ncond)), qnscl, 1. )
			CALL hbuf_put( 'QNG'//TRIM(condavgname(ncond)), qngcl, 1. )
			CALL hbuf_put( 'QNN'//TRIM(condavgname(ncond)), qnncl, 1. )
			CALL hbuf_put( 'W'//TRIM(condavgname(ncond))//'A', wacl, factor_xy )
			CALL hbuf_put( 'LIWSEW'//TRIM(condavgname(ncond)), twcl, factor_xy )
			CALL hbuf_put( 'TPVW'//TRIM(condavgname(ncond)), tvwcl, factor_xy )
			CALL hbuf_put( 'QTW'//TRIM(condavgname(ncond)), qwcl, factor_xy*1.E3 )
			CALL hbuf_put( 'QCW'//TRIM(condavgname(ncond)), qcwcl, factor_xy*1.E3 )
			CALL hbuf_put( 'QIW'//TRIM(condavgname(ncond)), qiwcl, factor_xy*1.E3 )
		ENDDO ! ncond = 1,ncondstats
		
	END SUBROUTINE statistics
	
	!======================= Interpolation subroutine =======================
	SUBROUTINE find_indices( z_interp_height, z_local_height, nz_in, km1, k00, kp1, kp2, klo )
	
		IMPLICIT NONE
		! IN
		INTEGER, INTENT(IN) :: nz_in
		REAL, INTENT(IN) :: z_interp_height
		REAL, DIMENSION(nz_in), INTENT(IN) :: z_local_height
		! OUT
		INTEGER, INTENT(OUT) :: km1, k00, kp1, kp2
		! INOUT
		INTEGER, INTENT(INOUT) :: klo
		! Local
		INTEGER :: khi, km, inc
		
		! Use the subroutine hunt in Numerical recipe to find k00.
		! Hunt has been modified so that it only works when
		! - z_local_height(nz_in) is ascend
		! - z_interp_height is highre or equal to the preiously z_interp_heigh
		! Both are true for this module.
		! Input klo has to be zero as an initial guess for the first time.
		! This subroutine is the most efficient if it is used in vertical level loop.
		
		!====================
		! Hunt
		! Given an array z_local_height(1:nz_in), and given a value z_interp_height, returns a value
		! klo such that z_interp_height is between z_local_height(klo) and z_local_height(klo+1).
		! z_local_height must be monotonically increasing. klo = 0 or klo = nz_in is returned to
		! indicate that z_interp_height is out of range. klo on input is taken as the initial guess
		! for klo on output.
		IF ( klo <= 0 .OR. klo > nz_in ) THEN
			! Input guess not useful. Go immediately to bisection.
			klo = 0
			khi = nz_in + 1
		ELSE
			! Set the hunting increment.
			inc = 1
			! Hunt up since new jlo is either current jlo or larger
			DO
				khi = klo + inc
				IF ( khi > nz_in ) THEN   ! Done hunting, since off end of table.
					khi = nz_in + 1
					EXIT
				ELSE
					IF ( z_interp_height < z_local_height(khi) ) EXIT
					! Not done hunting, so double the increment
					klo = khi
					inc = inc + inc
				ENDIF
			ENDDO
		ENDIF
		! Hunt is done, so begin the final bisection phase:
		DO
			IF ( khi - klo <= 1 ) THEN
				IF ( z_interp_height == z_local_height(nz_in) ) klo = nz_in - 1
				IF ( z_interp_height == z_local_height(1) )     klo = 1
				EXIT
			ELSE
				km = ( khi + klo ) / 2
				IF ( z_interp_height >= z_local_height(km) ) THEN
					klo = km
				ELSE
					khi = km
				ENDIF
			ENDIF
		ENDDO
		!====================
		
		! Assign k00 = klo
		k00 = klo
		
		! Assign km1, k00, kp1, kp2 to comunicate with mono_cubic_interp
		! At this point, 
		! - k00 = 0 [z_interp_height < z_local_height(1)]
		! - k00 = nz_in [z_interp_height > z_local_height(nz_in)]
		! - 1 <= k00 <= nz_in-1
		IF ( k00 == 0 ) THEN
		! interpolation height is below the lowest local height, indicate extrapolation for
		! mono_cubic_interp by setting km1 = nz_in
			km1 = nz_in
			k00 = 1
			kp1 = 2
			kp2 = 3
		ELSE IF ( k00 == nz_in ) THEN
		! interpolation height is above the highest local height, indicate extrapolation for
		! mono_cubic_interp by setting km1 = nz_in
			km1 = nz_in
			k00 = nz_in - 1
			kp1 = nz_in
			kp2 = nz_in
		ELSE IF ( k00 == 1 ) THEN
			km1 = 1
			kp1 = 2
			kp2 = 3
		ELSE IF ( k00 == nz_in-1 ) THEN
			km1 = nz_in - 2
			kp1 = nz_in
			kp2 = nz_in
		ELSE
			km1 = k00 - 1
			kp1 = k00 + 1
			kp2 = k00 + 2
		ENDIF
		
	END SUBROUTINE find_indices

	!-------------------------------------------------------
	
	FUNCTION mono_cubic_interp( z_in, km1, k00, kp1, kp2, zm1, z00, zp1, zp2, fm1, f00, fp1, fp2 ) &
		RESULT( f_out )
	
	! Steffen's monotonic cubic interpolation method
	! Returns monotonic cubic interpolated value between x00 and xp1
	
		IMPLICIT NONE
		! IN
		REAL, INTENT(IN) :: z_in
		INTEGER, INTENT(IN) :: km1, k00, kp1, kp2
		REAL, INTENT(IN) :: zm1, z00, zp1, zp2
		REAL, INTENT(IN) :: fm1, f00, fp1, fp2
		! OUT
		REAL :: f_out
		! Local
		REAL :: hm1, h00, hp1
		REAL :: sm1, s00, sp1
		REAL :: p00, pp1
		REAL :: dfdx00, dfdxp1
		REAL :: c1, c2, c3, c4
		REAL :: w00, wp1
		
		IF ( km1 <= k00 ) THEN
			hm1 = z00 - zm1
			h00 = zp1 - z00
			hp1 = zp2 - zp1
			IF ( km1 == k00 ) THEN
				s00 = ( fp1 - f00 ) / ( zp1 - z00 )
				sp1 = ( fp2 - fp1 ) / ( zp2 - zp1 )
				dfdx00 = s00
				pp1 = ( s00 * hp1 + sp1 * h00 ) / ( h00 + hp1 )
				dfdxp1 = ( SIGN(1.0,s00) + SIGN(1.0,sp1) ) * MIN(ABS(s00),ABS(sp1),0.5*ABS(pp1))
			ELSE IF ( kp1 == kp2 ) THEN
				sm1 = ( f00 - fm1 ) / ( z00 - zm1 )
				s00 = ( fp1 - f00 ) / ( zp1 - z00 )
				p00 = ( sm1 * h00 + s00 * hm1 ) / ( hm1 + h00 )
				dfdx00 = ( SIGN(1.0,sm1) + SIGN(1.0,s00) ) * MIN(ABS(sm1),ABS(s00),0.5*ABS(p00))
				dfdxp1 = s00
			ELSE
				sm1 = ( f00 - fm1 ) / ( z00 - zm1 )
				s00 = ( fp1 - f00 ) / ( zp1 - z00 )
				sp1 = ( fp2 - fp1 ) / ( zp2 - zp1 )
				p00 = ( sm1 * h00 + s00 * hm1 ) / ( hm1 + h00 )
				pp1 = ( s00 * hp1 + sp1 * h00 ) / ( h00 + hp1 )
				dfdx00 = ( SIGN(1.0,sm1) + SIGN(1.0,s00) ) * MIN(ABS(sm1),ABS(s00),0.5*ABS(p00))
				dfdxp1 = ( SIGN(1.0,s00) + SIGN(1.0,sp1) ) * MIN(ABS(s00),ABS(sp1),0.5*ABS(pp1))
			ENDIF
			c1 = ( dfdx00 + dfdxp1 - 2. * s00 ) / ( h00 ** 2 )
			c2 = ( 3. * s00 - 2. * dfdx00 - dfdxp1 ) / h00
			c3 = dfdx00
			c4 = f00
			f_out = c1 * ( (z_in - z00) ** 3 ) + c2 * ( (z_in - z00) ** 2 ) + c3 * (z_in - z00) + c4
		ELSE
			! Linear extrapolation
			wp1 = ( z_in - z00 ) / ( zp1 - z00 )
			w00 = 1. - wp1
			f_out = wp1 * fp1 + w00 * f00
		ENDIF
	
	END FUNCTION mono_cubic_interp
	
	!=============================== Functions in sat.f90 ================================
	REAL FUNCTION esatw(t)
		IMPLICIT NONE
		REAL, INTENT(IN) :: t ! temperature (K)
		REAL, PARAMETER :: a0 = 6.11239921
		REAL, PARAMETER :: a1 = 0.443987641
		REAL, PARAMETER :: a2 = 0.142986287E-1
		REAL, PARAMETER :: a3 = 0.264847430E-3
		REAL, PARAMETER :: a4 = 0.302950461E-5
		REAL, PARAMETER :: a5 = 0.206739458E-7
		REAL, PARAMETER :: a6 = 0.640689451E-10
		REAL, PARAMETER :: a7 = -0.952447341E-13
		REAL, PARAMETER :: a8 = -0.976195544E-15
		REAL :: dtmp
		dtmp = MAX(-80.,t-273.16)
		esatw = a0 + dtmp*(a1+dtmp*(a2+dtmp*(a3+dtmp*(a4+dtmp*(a5+dtmp*(a6+dtmp*(a7+a8*dtmp)))))))
	END FUNCTION esatw
	
	!-------------------------------------------------------
	REAL FUNCTION qsatw(t,p)
		IMPLICIT NONE
		REAL, INTENT(IN) :: t   ! temperature (K)
		REAL, INTENT(IN) :: p   ! pressure    (mb)
		REAL :: esat
		esat = esatw(t)
		qsatw = 0.622 * esat/MAX(esat,p-esat)
	END FUNCTION qsatw
	
	!============================ Subroutines in hbuffer.f90 =============================
	SUBROUTINE hbuf_average
		! Average the profiles in buffer
		IMPLICIT NONE
		INTEGER :: l
		REAL(8) :: coef
		coef = 1.0D0 / DBLE(nstatsteps)
		DO l = 1, hbuf_length*nzm
			hbuf(l) = hbuf(l) * coef
		ENDDO
	END SUBROUTINE hbuf_average
	
	!-------------------------------------------------------
	SUBROUTINE hbuf_flush
		! Flush the buffer
		IMPLICIT NONE
		INTEGER :: l, k, n
		DO l = 1, hbuf_length*nzm
			hbuf(l) = 0.
		ENDDO
		DO n = 1, ncondavg
			DO k = 1, nzm
				condavg_factor(k,n) = 0.
			ENDDO
		ENDDO
		s_sst   = 0.
		s_psfc  = 0.
		s_acld  = 0.
		s_acldl = 0.
		s_acldm = 0.
		s_acldh = 0.
		s_acldcold = 0.
		s_ar   = 0.
		s_prec = 0.
		s_precl = 0.
		s_preci = 0.
		s_pw   = 0.
		s_cwp  = 0.
		s_iwp  = 0.
		s_rwp  = 0.
		s_swp  = 0.
		s_gwp  = 0.
		s_lwp2 = 0.
		s_iwp2 = 0.
		w_max  = 0.
		w_min  = 0.
		u_max  = 0.
		u_min  = 0.
		v_max  = 0.
		v_min  = 0.
		s_lhf  = 0.
		s_shf  = 0.
		s_flns = 0.
		s_flnt = 0.
		s_flnsc = 0.
		s_flntc = 0.
		s_flds  = 0.
		s_fsns  = 0.
		s_fsnt  = 0.
		s_fsnsc = 0.
		s_fsntc = 0.
		s_fsds  = 0.
		s_solin = 0.
	END SUBROUTINE hbuf_flush
	
	!-------------------------------------------------------
	SUBROUTINE hbuf_init
		IMPLICIT NONE
		CHARACTER (LEN=15) :: nm
		CHARACTER (LEN=80) :: def
		CHARACTER (LEN=10) :: un
		INTEGER :: stat, count, type, ios, trcount, i, m, n
		LOGICAL :: duplicate_entries, file_exist
		
		! Read list
		OPEN(66, FILE='./stat/lst', STATUS='old', FORM='formatted')
		! Determine how many entries will be:
		hbuf_length = 0
		DO i = 1, hbuf_max_length
			READ(66, IOSTAT=ios, FMT=*) stat, type, nm
			IF (ios /= 0) EXIT
			IF (stat > 0) hbuf_length = hbuf_length + 1
		ENDDO
		IF (hbuf_length > hbuf_max_length) THEN
			PRINT *,'Error: hbuf_length > hbuf_max_length.'
			CALL task_abort()
		ENDIF
		
		! Fill the buffers:
		REWIND(66)
		count = 0
		DO i = 1, hbuf_length
			READ(66, IOSTAT=ios, FMT=*) stat, type, nm, def, un
			IF (ios /= 0 ) EXIT
			IF (stat > 0) THEN
				count = count + 1
				namelist(count) = nm
				deflist(count) = def
				unitlist(count) = un
				status(count) = stat
				average_type(count) = type
			ENDIF
		ENDDO
		
		! Add conditional sampled variables
		trcount = 0
		CALL hbuf_conditionals_init(count,trcount)
		hbuf_length = hbuf_length + trcount
		
		! Check if there are duplicate entries in the stat list:
		duplicate_entries = .FALSE.
		DO n = 1, count-1
			DO m = n+1, count
				IF ( TRIM(namelist(n)) == TRIM(namelist(m)) ) THEN
					duplicate_entries = .TRUE.
					PRINT *,'Error: Multiple definition of '//namelist(n)// ' variable in stat list'
				ENDIF
			ENDDO
		ENDDO
		
		! Halt simulation if duplicate entries appear in namelist.
		IF (duplicate_entries) CALL task_abort()
		
		! Print list of stat variables
		PRINT *,'Number of statistics profiles:', hbuf_length
		PRINT *,'Statistics profiles to save:'
		WRITE(*,'(5(A,3X))')(namelist(i), i=1,hbuf_length)
		IF (masterproc) THEN
			! Make sure that the stat file doesn't exist if a new run to prevent accidental overwrite
			INQUIRE( FILE='./stat/'//TRIM(stat_file_name), EXIST=file_exist )
			IF (file_exist.AND.nstep == 0) THEN
				PRINT*,'A same data file exists. Abort.'
				CALL task_abort
			ELSE IF (file_exist.AND.nstep /= 0) THEN
				PRINT*,'Restart run. A data file will be appended.'
			ELSE IF (.NOT.file_exist.AND.nstep == 0) THEN
				PRINT*,'Initial run. A data file will be created.'
				OPEN (55, FILE='./stat/'//TRIM(stat_file_name), STATUS='new', FORM='unformatted')
				CLOSE(55)
			ELSE IF (.NOT.file_exist.AND.nstep /= 0) THEN
				PRINT*,'A data file to be appended does not exist. Abort.'
				CALL task_abort
			ENDIF
		ENDIF
		
		! Close lst
		CLOSE(66)
		
		! Flush the buffer
		CALL hbuf_flush()
	
	END SUBROUTINE hbuf_init
	
	!-------------------------------------------------------
	SUBROUTINE hbuf_put(name, f, factor)
		! Write to buffer (not to file yet)
		IMPLICIT NONE
		REAL, DIMENSION(nzm), INTENT(IN) :: f
		REAL, INTENT(IN) :: factor
		CHARACTER *(*) name
		INTEGER :: l, begin, k
		LOGICAL :: flag
		flag = .FALSE.
		DO l = 1, hbuf_length
			IF (.NOT.(lgt(name,namelist(l)).OR.llt(name,namelist(l)))) THEN
				flag = .TRUE.
				IF (status(l) > 0) THEN
					status(l) = 999
					begin = (l-1)*nzm
					DO k = 1, nzm
						hbuf(begin+k) = hbuf(begin+k) + f(k)*factor
					ENDDO
				ENDIF
			ENDIF
		ENDDO
		IF (.NOT.flag.AND.masterproc) print*,'name ', name,' is missing in "lst" file.'
	END SUBROUTINE hbuf_put
	
	!-------------------------------------------------------
	SUBROUTINE hbuf_write
	
		! Write buffer to file
		
		IMPLICIT NONE
		
		INTEGER :: l, k, ntape, length, nbuf
		REAL(8) :: coef, aver, factor
		REAL, DIMENSION(nzm) :: tmp, tmp1
		REAL :: day
		REAL(8), DIMENSION(HBUF_MAX_LENGTH*nzm+40) :: hbuf1
		REAL(4), DIMENSION(HBUF_MAX_LENGTH*nzm+40) :: hbuf_real4
		REAL, DIMENSION(nzm,ncondavg) :: cloud_f   !UNUSED=>!, tmp_f !bloss
		REAL, DIMENSION(nzm*ncondavg) :: condavg_factor1D, cloud_f1D
		INTEGER :: ncond !bloss
!UNUSED!		INTEGER :: nsteplast
!UNUSED!		REAL(4) :: dummy_real4 !bloss
		REAL(4) :: time_real4(2) !bloss
		REAL(4) :: deltas_real4(3) !bloss
		REAL(4) :: adz_real4(nzm) !bloss
		REAL(4) :: z_real4(nz) !bloss
		REAL(4) :: pres_real4(nzm) !bloss
		REAL(4) :: misc_real4(40) !bloss
		
		DATA ntape/55/
		
		aver = 1.0D0 / DBLE(nstatsteps)
		factor = 1.0D0 / DBLE(nx*ny)
		
		IF (dompi) THEN
			! average condavg_factor across domains.  This will sum the weighting of all of the
			! conditional statistics across the processors and allow us to normalize the statistics
			! on each processor accordingly.  In the end, these normalized statistics on each
			! processor will sum to give the value of the conditional statistic.
			
			! sum across processors
			! Put data into a 1D array
			DO ncond = 1, ncondavg
				condavg_factor1D( 1+(ncond-1)*nzm:ncond*nzm ) = condavg_factor(1:nzm,ncond)
			ENDDO
			! MPI
			CALL task_sum_real(condavg_factor1D, cloud_f1D, nzm*ncondavg)
			! Put data into a 2D array
			DO ncond = 1, ncondavg
				cloud_f(1:nzm,ncond) = cloud_f1D( 1+(ncond-1)*nzm:ncond*nzm )
			ENDDO
			! The following code crashes sometimes, especially with larger number of processors.
			! Passing 2D arrays may be problematic. The other MPI calls with 1D arrays has never
			! crashed... All subroutines, task_#, are modified so that they only accept 1D array.
			! This is probably good way since assumed shaped array usd before is obsolete.
		!!	CALL task_sum_real(condavg_factor, cloud_f, nzm*ncondavg)
		ELSE
			cloud_f(1:nzm, 1:ncondavg) = condavg_factor(1:nzm, 1:ncondavg)
		ENDIF
		
		! Create normalization/weighting factor for conditional statistics
		! Here, cloud_f holds the sum of condavg_factor across all processors.
		condavg_factor(:,:) = 0.
		DO ncond = 1, ncondavg
			DO k = 1, nzm
				IF (ABS(cloud_f(k,ncond)) > EPSILON(1.)) THEN
					condavg_factor(k,ncond) = FLOAT(nsubdomains)*FLOAT(nstatsteps)/cloud_f(k,ncond)
				ENDIF
			ENDDO
		ENDDO
		
		! Compute each processor's component of the total conditional average.
		length = 0
		DO l = 1, hbuf_length
			IF (status(l) == 999) THEN
				length = length + 1
				DO ncond = 1, ncondavg
					IF (average_type(l) == ncond) THEN
						DO k = 1, nzm
							hbuf((l-1)*nzm+k) = hbuf((l-1)*nzm+k)*condavg_factor(k,ncond)
						ENDDO
						DO k = 1, nzm
							! Insert a missing value if there are no samples of this conditional
							! statistic at this level.
							IF (ABS(cloud_f(k,ncond)) < EPSILON(1.)) hbuf((l-1)*nzm+k) = -9999.
						ENDDO
					ENDIF
				ENDDO
			ENDIF
		ENDDO
		
		! Get statistics buffer from different processes, add them together and average
		IF (dompi) THEN
			coef = 1.0D0 / DBLE(nsubdomains)
			
			tmp1(1) = w_max
			tmp1(2) = u_max
			tmp1(3) = v_max
			CALL task_max_real(tmp1,tmp,3)
			w_max = tmp(1)
			u_max = tmp(2)
			v_max = tmp(3)
			
			tmp1(1) = w_min
			tmp1(2) = u_min
			tmp1(3) = v_min
			CALL task_min_real(tmp1,tmp,3)
			w_min = tmp(1)
			u_min = tmp(2)
			v_min = tmp(3)
			
			k = hbuf_length*nzm
			hbuf(k+1)  = s_sst
			hbuf(k+2)  = s_psfc
			hbuf(k+3)  = s_acld
			hbuf(k+4)  = s_acldl
			hbuf(k+5)  = s_acldm
			hbuf(k+6)  = s_acldh
			hbuf(k+7)  = s_acldcold
			hbuf(k+8)  = s_ar
			hbuf(k+9)  = s_prec
			hbuf(k+10) = s_precl
			hbuf(k+11) = s_preci
			hbuf(k+12) = s_pw
			hbuf(k+13) = s_cwp
			hbuf(k+14) = s_iwp
			hbuf(k+15) = s_rwp
			hbuf(k+16) = s_swp
			hbuf(k+17) = s_gwp
			hbuf(k+18) = s_lwp2
			hbuf(k+19) = s_iwp2
			hbuf(k+20) = s_lhf
			hbuf(k+21) = s_shf
			hbuf(k+22) = s_flns
			hbuf(k+23) = s_flnt
			hbuf(k+24) = s_flnsc
			hbuf(k+25) = s_flntc
			hbuf(k+26) = s_flds
			hbuf(k+27) = s_fsns
			hbuf(k+28) = s_fsnt
			hbuf(k+29) = s_fsnsc
			hbuf(k+30) = s_fsntc
			hbuf(k+31) = s_fsds
			hbuf(k+32) = s_solin
			nbuf = k + 32   ! change value whenever the number of hbuf is changed
			CALL task_sum_real8(hbuf,hbuf1,nbuf)
			hbuf(1:nbuf) = hbuf1(1:nbuf)*coef
			s_sst   = hbuf(k+1)
			s_psfc  = hbuf(k+2)
			s_acld  = hbuf(k+3)
			s_acldl = hbuf(k+4)
			s_acldm = hbuf(k+5)
			s_acldh = hbuf(k+6)
			s_acldcold = hbuf(k+7)
			s_ar    = hbuf(k+8)
			s_prec  = hbuf(k+9)
			s_precl = hbuf(k+10)
			s_preci = hbuf(k+11)
			s_pw    = hbuf(k+12)
			s_cwp   = hbuf(k+13)
			s_iwp   = hbuf(k+14)
			s_rwp   = hbuf(k+15)
			s_swp   = hbuf(k+16)
			s_gwp   = hbuf(k+17)
			s_lwp2  = hbuf(k+18)
			s_iwp2  = hbuf(k+19)
			s_lhf   = hbuf(k+20)
			s_shf   = hbuf(k+21)
			s_flns  = hbuf(k+22)
			s_flnt  = hbuf(k+23)
			s_flnsc = hbuf(k+24)
			s_flntc = hbuf(k+25)
			s_flds  = hbuf(k+26)
			s_fsns  = hbuf(k+27)
			s_fsnt  = hbuf(k+28)
			s_fsnsc = hbuf(k+29)
			s_fsntc = hbuf(k+30)
			s_fsds  = hbuf(k+31)
			s_solin = hbuf(k+32)
		ENDIF
		
		IF (masterproc) THEN
			
			IF (nstep == 0) THEN
				! WRF saves initial profile, so save statistics, too
				OPEN(ntape, FILE='./stat/'//TRIM(stat_file_name), STATUS='old', FORM='unformatted')
				day = 0.
			ELSE
				OPEN(ntape, FILE='./stat/'//TRIM(stat_file_name), STATUS='old', FORM='unformatted', &
				     POSITION='append')
				day = REAL(nstep) * dt/86400.
				IF (nstatfrq > 1) THEN
					! time should be middle for the time average period.
					day = day - REAL(nstat)/2.*dt/86400.
				ENDIF
			ENDIF
			PRINT *,'Writing history file ', TRIM(stat_file_name)
			WRITE(ntape) stat_file_name
			
			!bloss: output in single precision.
			time_real4(1) = day   ! day-nstat*dt/2./86400., & ! day is computed above
			time_real4(2) = dt
			deltas_real4(1) = dx
			deltas_real4(2) = dy
			deltas_real4(3) = 1.   ! dz, & ! write 1 instead of dz
			adz_real4(1:nzm) = zs(1:nzm)   ! adz, & there is not adz(1:nzm)
			z_real4(1:nz) = zs(1:nz)   ! Vertical height for output
			pres_real4(1:nzm) = press(1:nzm)   ! initial hydrostatic pressure level
			
			! 1D time series
			! Index number corresponds to npar in stat2nc.f
			! Whenever the list below is modified, stat2nc.F in ./stat/stat2nc/SRC has to be modified.
			! 1. SST
			misc_real4(1) = s_sst * aver * factor
			
			! 2. Surface pressure
			misc_real4(2) = s_psfc * aver * factor
			
			! 3. Shaded cloud fraction
			misc_real4(3) = s_acld * aver * factor
			
			! 4. Low cloud fraction
			misc_real4(4) = s_acldl * aver * factor
			
			! 5. Middle cloud fraction
			misc_real4(5) = s_acldm * aver * factor
			
			! 6. High cloud fraction
			misc_real4(6) = s_acldh * aver * factor
			
			! 7. Cloud fraction above 245 K level
			misc_real4(7) = s_acldcold * aver * factor
			
			! 8. Surface precipitation fraction
			misc_real4(8) = s_ar * aver * factor
			
			! 9-11. Surface precipitation
			misc_real4(9) = s_prec * aver * factor
			misc_real4(10) = s_precl * aver * factor
			misc_real4(11) = s_preci * aver * factor
			
			! 12. Convective precipitation rate
			! 13. Stratiform precipitation rate
!!!!! Not yet implemented !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			misc_real4(12) = -1.
			misc_real4(13) = -1.
			
			! 14. Precipitable water
			misc_real4(14) = s_pw * aver * factor
			
			! 15. Cloud water path
			misc_real4(15) = s_cwp * aver * factor
			
			! 16. Ice water path
			misc_real4(16) = s_iwp * aver * factor
			
			! 17. Rain water path
			misc_real4(17) = s_rwp * aver * factor
			
			! 18. Snow water path
			misc_real4(18) = s_swp * aver * factor
			
			! 19. Grauple water path
			misc_real4(19) = s_gwp * aver * factor
			
			! 20. Variance of liquid water path
			misc_real4(20) = s_lwp2 * aver * factor
			
			! 21. Variance of ice water path
			misc_real4(21) = s_iwp2 * aver * factor
			
			! 22-27. Maximum updraft, downdraft, x wind, negative x wind, y wind, negative y wind
			misc_real4(22) = w_max
			misc_real4(23) = w_min
			misc_real4(24) = u_max
			misc_real4(25) = u_min
			misc_real4(26) = v_max
			misc_real4(27) = v_min
				
			! 28. Surface latent heat flux
			misc_real4(28) = s_lhf * aver * factor
				
			! 29. Surface sensible heat flux
			misc_real4(29) = s_shf * aver * factor
				
			! 30. Net LW flux at sfc
!!!!! Not yet implemented !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! See module_first_rk_step_part1.F ==> module_radiation_driver.F !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			misc_real4(30) = -1.   ! s_flns*aver*factor
				
			! 31. Net LW flux at TOA
!!!!! Not yet implemented !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			misc_real4(31) = -1.   ! s_flnt*aver*factor
				
			! 32. Net LW flux at sfc (clear sky)
!!!!! Not yet implemented !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			misc_real4(32) = -1.   ! s_flnsc*aver*factor
				
			! 33. Net LW flux at TOA (clear sky)
!!!!! Not yet implemented !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			misc_real4(33) = -1.   ! s_flntc*aver*factor
				
			! 34. Downward LW flux at sfc
!!!!! Not yet implemented !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			misc_real4(34) = -1.   ! s_flds*aver*factor
				
			! 35. Net SW flux at sfc
!!!!! Not yet implemented !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			misc_real4(35) = -1.   ! s_fsns*aver*factor
				
			! 36. Net SW flux at TOA
!!!!! Not yet implemented !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			misc_real4(36) = -1.   ! s_fsnt*aver*factor
				
			! 37. Net SW flux at sfc (clear sky)
!!!!! Not yet implemented !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			misc_real4(37) = -1.   ! s_fsnsc*aver*factor
			
			! 38. Net SW flux at TOA (clear sky)
!!!!! Not yet implemented !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			misc_real4(38) = -1.   ! s_fsntc*aver*factor
			
			! 39. Downward SW flux at sfc
!!!!! Not yet implemented !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			misc_real4(39) = -1.   ! s_fsds*aver*factor
			
			! 40. Incoming SW flux at TOA
!!!!! Not yet implemented !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			misc_real4(40) = -1.   ! s_solin*aver*factor
			
			WRITE(ntape) (time_real4(k),k=1,2), nstep, nx, ny, nz, nzm, &
			             (deltas_real4(k),k=1,3), adz_real4, z_real4, pres_real4, misc_real4
			
			WRITE(ntape) length
			hbuf_real4(1:hbuf_length*nzm) = hbuf(1:hbuf_length*nzm)
			DO l = 1, hbuf_length
				IF (status(l) == 999) THEN
					WRITE(ntape) namelist(l)
					WRITE(ntape) deflist(l)
					WRITE(ntape) unitlist(l)
					WRITE(ntape) (hbuf_real4((l-1)*nzm+k),k=1,nzm)
				ENDIF
			ENDDO
			CLOSE(ntape)
		
		ENDIF
	
	END SUBROUTINE hbuf_write
	
	!===================== Subroutines in hbuf_conditionals_init.f90 =====================
	SUBROUTINE hbuf_conditionals_init(count, trcount)
		IMPLICIT NONE
		! Initialize the list of UW statistics variables written in statistics.f90
		INTEGER :: count,trcount, n
		!bloss: setup to add an arbitrary number of conditional statistics
		DO n = 1, ncondavg
			!bloss: add all of the conditional statistics here, so that they don't have to be added &
			!       to the lst file
			CALL add_to_namelist(count,trcount,TRIM(condavgname(n)), &
				TRIM(condavglongname(n))//' fraction',' ',0)
			CALL add_to_namelist(count,trcount,'RHO'//TRIM(condavgname(n)), &
				'Mean RHO in '//TRIM(condavglongname(n)),'kg/m3',n)
			CALL add_to_namelist(count,trcount,'U'//TRIM(condavgname(n)), &
				'Mean U in '//TRIM(condavglongname(n)),'m/s',n)
			CALL add_to_namelist(count,trcount,'V'//TRIM(condavgname(n)), &
				'Mean V in '//TRIM(condavglongname(n)),'m/s',n)
			CALL add_to_namelist(count,trcount,'W'//TRIM(condavgname(n)), &
				'Mean W in '//TRIM(condavglongname(n)),'m/s',n)
			CALL add_to_namelist(count,trcount,'DSE'//TRIM(condavgname(n)), &
				'Mean dry static energy in '//TRIM(condavglongname(n)),'K',n)
			CALL add_to_namelist(count,trcount,'LIWSE'//TRIM(condavgname(n)), &
				'Mean liquid-ice static energy in '//TRIM(condavglongname(n)),'K',n)
			CALL add_to_namelist(count,trcount,'MSE'//TRIM(condavgname(n)), &
				'Mean moist static energy in '//TRIM(condavglongname(n)),'K',n)
			CALL add_to_namelist(count,trcount,'TABS'//TRIM(condavgname(n)), &
				'Mean TABS in '//TRIM(condavglongname(n)),'K',n)
			CALL add_to_namelist(count,trcount,'TPV'//TRIM(condavgname(n)), &
				'Mean TPV in '//TRIM(condavglongname(n)),'K',n)
			CALL add_to_namelist(count,trcount,'TPV'//TRIM(condavgname(n))//'A', &
				'Mean THETAV anomaly in '//TRIM(condavglongname(n)),'K',n)
			CALL add_to_namelist(count,trcount,'QT'//TRIM(condavgname(n)), &
				'Mean QT in '//TRIM(condavglongname(n)),'g/kg',n)
			CALL add_to_namelist(count,trcount,'QC'//TRIM(condavgname(n)), &
				'Mean QC in '//TRIM(condavglongname(n)),'g/kg',n)
			CALL add_to_namelist(count,trcount,'QI'//TRIM(condavgname(n)), &
				'Mean QI in '//TRIM(condavglongname(n)),'g/kg',n)
			CALL add_to_namelist(count,trcount,'QN'//TRIM(condavgname(n)), &
				'Mean QN in '//TRIM(condavglongname(n)),'g/kg',n)
			CALL add_to_namelist(count,trcount,'QP'//TRIM(condavgname(n)), &
				'Mean QP in '//TRIM(condavglongname(n)),'g/kg',n)
			CALL add_to_namelist(count,trcount,'ETA'//TRIM(condavgname(n)), &
				'Mean ETA in '//TRIM(condavglongname(n)),'g/kg',n)
			CALL add_to_namelist(count,trcount,'SSR'//TRIM(condavgname(n)), &
				'Mean SSR in '//TRIM(condavglongname(n)),'per cent',n)
			CALL add_to_namelist(count,trcount,'QNC'//TRIM(condavgname(n)), &
				'Mean QNC in '//TRIM(condavglongname(n)),'#/mg',n)
			CALL add_to_namelist(count,trcount,'QNI'//TRIM(condavgname(n)), &
				'Mean QNI in '//TRIM(condavglongname(n)),'#/mg',n)
			CALL add_to_namelist(count,trcount,'QNR'//TRIM(condavgname(n)), &
				'Mean QNR in '//TRIM(condavglongname(n)),'#/mg',n)
			CALL add_to_namelist(count,trcount,'QNS'//TRIM(condavgname(n)), &
				'Mean QNS in '//TRIM(condavglongname(n)),'#/mg',n)
			CALL add_to_namelist(count,trcount,'QNG'//TRIM(condavgname(n)), &
				'Mean QNG in '//TRIM(condavglongname(n)),'#/mg',n)
			CALL add_to_namelist(count,trcount,'QNN'//TRIM(condavgname(n)), &
				'Mean QNN in '//TRIM(condavglongname(n)),'#/mg',n)
			CALL add_to_namelist(count,trcount,'W'//TRIM(condavgname(n))//'A', &
				'W in '//TRIM(condavglongname(n))//' averaged over the whole domain','m/s',0)
			CALL add_to_namelist(count,trcount,'LIWSEW'//TRIM(condavgname(n)), &
				'LIWSExW in '//TRIM(condavglongname(n))//' averaged over the whole domain','Km/s',0)
			CALL add_to_namelist(count,trcount,'TPVW'//TRIM(condavgname(n)), &
				'TPVxW in '//TRIM(condavglongname(n))//' averaged over the whole domain','Km/s',0)
			CALL add_to_namelist(count,trcount,'QTW'//TRIM(condavgname(n)), &
				'QTxW in '//TRIM(condavglongname(n))//' averaged over the whole domain','g/kg m/s',0)
			CALL add_to_namelist(count,trcount,'QCW'//TRIM(condavgname(n)), &
				'QCxW in '//TRIM(condavglongname(n))//' averaged over the whole domain','g/kg m/s',0)
			CALL add_to_namelist(count,trcount,'QIW'//TRIM(condavgname(n)), &
				'QIxW in '//TRIM(condavglongname(n))//' averaged over the whole domain','g/kg m/s',0)
	  ENDDO ! n = 1,ncondavg
	END SUBROUTINE hbuf_conditionals_init
	
	!-------------------------------------------------------
	SUBROUTINE add_to_namelist(count,trcount,varname,varlongname,varunits,varavg)
		! Add variable to namelist
		IMPLICIT NONE
		INTEGER count, trcount, varavg
		CHARACTER(*) varname
		CHARACTER(*) varlongname
		CHARACTER(*) varunits
		count = count + 1
		trcount = trcount + 1
		namelist(count) = trim(varname)
		deflist(count) = trim(varlongname)
		unitlist(count) = trim(varunits)
		status(count) = 1
		average_type(count) = varavg
	END SUBROUTINE add_to_namelist
	
	!============== MPI utility copied form SAM6.8.2/SRC/task_util_MPI.f90 ===============
	SUBROUTINE task_init
		USE module_dm, ONLY: local_communicator
		IMPLICIT NONE
		INCLUDE 'mpif.h'
		INTEGER :: ierr
#ifndef DM_PARALLEL
		dompi = .FALSE.
		rank = 0
		nsubdomains = 1
#else
		mpi_comm = local_communicator
		dompi = .TRUE.
		CALL MPI_COMM_RANK(mpi_comm, rank, ierr)
		CALL MPI_COMM_SIZE(mpi_comm, nsubdomains, ierr)
#endif
		masterproc = rank == 0
	END SUBROUTINE task_init
	!-------------------------------------------------------
	SUBROUTINE task_abort
		IMPLICIT NONE
		INCLUDE 'mpif.h'
		INTEGER :: ierr, rc
		IF (dompi) THEN
			CALL MPI_ABORT(mpi_comm, rc, ierr)
			CALL MPI_FINALIZE(ierr)
		ENDIF
		CALL EXIT(999)   ! to avolid resubmission when finished
	END SUBROUTINE task_abort
	!-------------------------------------------------------
	SUBROUTINE task_barrier
		IMPLICIT NONE
		INCLUDE 'mpif.h'
		INTEGER :: ierr
		IF (dompi) CALL MPI_BARRIER(mpi_comm, ierr)
	END SUBROUTINE task_barrier
	!-------------------------------------------------------
	SUBROUTINE task_sum_real(buffer_in, buffer_out, length)
		IMPLICIT NONE
		INCLUDE 'mpif.h'
		INTEGER, INTENT(IN) :: length   ! buffers' length
		REAL, DIMENSION(length), INTENT(IN)  :: buffer_in
		REAL, DIMENSION(length), INTENT(OUT) :: buffer_out   ! buffer of data
		INTEGER :: ierr
		CALL MPI_ALLREDUCE(buffer_in, buffer_out, length, MPI_REAL, MPI_SUM, mpi_comm, ierr)
	END SUBROUTINE task_sum_real
	!-------------------------------------------------------
	SUBROUTINE task_sum_real8(buffer_in, buffer_out, length)
		IMPLICIT NONE
		INCLUDE 'mpif.h'
		INTEGER, INTENT(IN) :: length
		DOUBLE PRECISION, DIMENSION(length), INTENT(IN)  :: buffer_in
		DOUBLE PRECISION, DIMENSION(length), INTENT(OUT) :: buffer_out
		INTEGER :: ierr
		CALL MPI_ALLREDUCE(buffer_in,buffer_out,length,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm,ierr)
	END SUBROUTINE task_sum_real8
	!-------------------------------------------------------
	SUBROUTINE task_max_real(buffer_in, buffer_out, length)
		IMPLICIT NONE
		INCLUDE 'mpif.h'
		INTEGER, INTENT(IN) :: length
		REAL, DIMENSION(length), INTENT(IN)  :: buffer_in
		REAL, DIMENSION(length), INTENT(OUT) :: buffer_out
		INTEGER :: ierr
		CALL MPI_ALLREDUCE(buffer_in, buffer_out, length, MPI_REAL, MPI_MAX, mpi_comm, ierr)
	END SUBROUTINE task_max_real
	!-------------------------------------------------------
	SUBROUTINE task_min_real(buffer_in, buffer_out, length)
		IMPLICIT NONE
		INCLUDE 'mpif.h'
		INTEGER, INTENT(IN) :: length
		REAL, DIMENSION(length), INTENT(IN)  :: buffer_in
		REAL, DIMENSION(length), INTENT(OUT) :: buffer_out
		INTEGER :: ierr
		CALL MPI_ALLREDUCE(buffer_in, buffer_out, length, MPI_REAL, MPI_MIN, mpi_comm, ierr)
	END SUBROUTINE task_min_real
	
	!=====================================================================================
	
END MODULE module_statistics
#endif /*WRFSTAT*/
