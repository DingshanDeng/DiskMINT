   MODULE global_variables

   PUBLIC
   
   INTEGER, PUBLIC, PARAMETER :: dp = SELECTED_REAL_KIND(p=15)
   
   INTEGER :: nspe
   INTEGER :: nspegas
   INTEGER :: nspegr
   
   INTEGER :: nreac
   INTEGER :: nreacgas
   INTEGER :: nreacgr
   
   !INTEGER, PUBLIC, PARAMETER :: nelem = 13
   INTEGER, PUBLIC, PARAMETER :: nelem = 15 ! When isotopes chemistry is taken into account (13C, 18O)
   INTEGER, PUBLIC, PARAMETER :: ncomp = 8
   
   LOGICAL :: first_step_done = .false.
   
   ! For DLSODES
   INTEGER                                  :: lrw                             ! declared length of RWORK (in user dimension).
   INTEGER                                  :: liw                             ! declared length of IWORK (in user dimension).
   INTEGER, DIMENSION(:), allocatable       :: iwork                           ! integer work array of length at least 30.
   REAL(kind=dp), DIMENSION(:), allocatable :: rwork                           ! real work array
   INTEGER                                  :: nb_nonzeros_values              ! number of non-zeros values in the jacobian. This is usefull for ODEPACK, to increase speed
   REAL(kind=dp)                            :: relative_tolerance = 1.0e-4_dp
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: satol ! Array that contain the absolute tolerance,
                                                     ! one value per species, either a minimum value
                                                     ! or a value related to its abundance.
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: ian, jan
   INTEGER :: itol = 2                ! 1 or 2 according as ATOL (below) is a scalar or array.
   INTEGER :: itask = 1               ! 1 for normal computation of output values of Y at t = TOUT.
   INTEGER :: iopt = 1                ! 0 to indicate no optional inputs used.
   INTEGER :: mf = 121                ! Method flag
   REAL(kind=dp) :: atol = 1.0e-40_dp ! integrator tolerance
   INTEGER :: istate
   
   INTEGER, DIMENSION(:), ALLOCATABLE :: idx_const_rates
   INTEGER, DIMENSION(:), ALLOCATABLE :: idx_var_rates
   
   !----------------------------------------------------------
   ! FLAGS OR SWITCHES
   !----------------------------------------------------------
   INTEGER            :: is_alloutput              = 0! If -1 then we save evrything
   INTEGER            :: is_grain_reactions        = 1
   INTEGER, PARAMETER :: is_absorption             = 1
   INTEGER, PARAMETER :: is_h2_adhoc_form          = 1
   INTEGER, PARAMETER :: grain_tunneling_diffusion = 0
   INTEGER, PARAMETER :: conservation_type         = 0
   INTEGER, PARAMETER :: is_er_cir                 = 0
   INTEGER, PARAMETER :: is_photodesorb            = 1
   INTEGER, PARAMETER :: is_reac_diff              = 1
   INTEGER            :: is_growth
   !----------------------------------------------------------
   
   !----------------------------------------------------------
   ! PHYSICAL CONSTANTS
   !----------------------------------------------------------
   REAL(KIND=dp), PARAMETER :: k_b                      = 1.3806488e-16_dp     !  Boltzmann constant in CGS (cm^2 g s^⁻2 K-1)
   REAL(KIND=dp), PARAMETER :: planck_constant          = 6.62565e-27_dp       !  Planck constant in CGS (erg.s or g cm2 s-1)
   REAL(KIND=dp), PARAMETER :: clight                   = 2.99792458e10_dp     !  speed of light in CGS (cm/s)
   REAL(KIND=dp), PARAMETER :: pi                       = 3.1415926535898_dp   !  The number Pi
   REAL(KIND=dp), PARAMETER :: sqpi                     = 1.77245385090551602_dp  ! sqrt(pi)
   REAL(KIND=dp), PARAMETER :: h_barre                  = 1.054571628E-27_dp   !  Reduced Planck constant h/2*pi in CGS (g cm2 s-1)
   REAL(KIND=dp), PARAMETER :: hplanck                  = 6.6260755E-27_dp
   REAL(KIND=dp), PARAMETER :: amu                      = 1.66053892e-24_dp    !  Atomic mass unit in g
   REAL(KIND=dp), PARAMETER :: electron_mass            = 0.000548579909_dp    !  Electron mass in AMU (close to 1/1836)
   REAL(KIND=dp), PARAMETER :: avogadro                 = 6.02214129e+23_dp    !  avogadro number : number of atom in 1 mol
   REAL(KIND=dp), PARAMETER :: year                     = 3.15576e+7_dp        !  one year in seconds
   REAL(KIND=dp), PARAMETER :: au                       = 1.49597871e+13_dp    !  Astronomical unit in cm (mean earth-sun distance)
   REAL(KIND=DP), PARAMETER :: EVerg                    = 1.60218e-12_dp! 1 eV       =  1.60218D-12 erg
   REAL(KIND=DP), PARAMETER :: me                       = 9.10938291e-28_dp       ! mass of the electron (g)
   REAL(KIND=DP), PARAMETER :: qe                       = 4.80320450571347e-10_dp ! charge of the electron (esu)
   !----------------------------------------------------------
   
   !----------------------------------------------------------
   ! PHYSICAL PARAMETERS
   !----------------------------------------------------------
   REAL(KIND=dp) :: Rstar = 6.957e10_dp
   REAL(KIND=dp) :: scale_G0 = 1.0_dp
   REAL(KIND=dp) :: density = 2.0e4_dp
   REAL(KIND=dp) :: ndust = 1.0e-99_dp
   REAL(KIND=dp) :: ndusta = 1.0e-99_dp
   REAL(KIND=dp) :: ndusta2 = 1.0e-99_dp
   REAL(KIND=dp) :: Tgas = 10.0_dp
   REAL(KIND=dp) :: Tdust = 10.0_dp
   REAL(KIND=dp) :: Av = 1.0e1_dp
   REAL(KIND=dp) :: Avup = 1.0e-99_dp
   REAL(KIND=dp) :: Avdn = 1.0e-99_dp
   REAL(KIND=dp) :: CRrate = 1.8e-17_dp
   
   !----------------------------------------------------------
   ! UV FIELD AND TAUS
   !----------------------------------------------------------
   REAL(KIND=dp) :: G0Hab = 1.0e-99_dp      
   REAL(KIND=dp) :: G0upHab = 1.0_dp      
   REAL(KIND=dp) :: G0dnHab = 1.0_dp       
   REAL(KIND=dp) :: Tauh2 = 1.0e-99_dp
   REAL(KIND=dp) :: Tauh2v = 1.0e-99_dp
   REAL(KIND=dp) :: Tauuph2 = 1.0e-99_dp
   REAL(KIND=dp) :: Tauuph2v = 1.0e-99_dp
   REAL(KIND=dp) :: Taudnh2 = 1.0e-99_dp
   REAL(KIND=dp) :: Taudnh2v = 1.0e-99_dp
      
   !----------------------------------------------------------
   
   !----------------------------------------------------------
   ! GRAINS PARAMETERS
   !----------------------------------------------------------
   ! - Alocated parameters
   REAL(KIND=dp), PARAMETER :: grain_density               =  3.0_dp
   REAL(KIND=dp), PARAMETER :: diff_binding_ratio_surf     =  0.4_dp
   REAL(KIND=dp), PARAMETER :: diff_binding_ratio_mant     =  0.8_dp
   REAL(KIND=dp), PARAMETER :: diffusion_barrier_thickness =  1.000e-08_dp ! Barrier thickness [cm]
   REAL(KIND=dp), PARAMETER :: ds                          =  2.6e-08_dp   ! Distance between 2 sites [cm]
   REAL(KIND=dp), PARAMETER :: surface_site_density        =  1.500e+15_dp ! site density [cm-2]
   REAL(KIND=dp), PARAMETER :: chemical_barrier_thickness  =  1.000e-08_dp ! grain reaction activation energy barrier width. [cm]
   REAL(KIND=dp), PARAMETER :: cr_peak_grain_temp          =  7.000E+01_dp ! peak grain temperature [K] (CR heating)
   REAL(KIND=dp), PARAMETER :: cr_peak_duration            =  1.000e-05_dp ! duration [s] of peak grain temperature
   REAL(KIND=dp), PARAMETER :: Fe_ionisation_rate          =  3.000e-14_dp ! (cosmic) Fe-ion--grain encounter [s-1 grain-1]
                                                                           ! (for 0.1 micron grain) For cosmic photo desorptions,
                                                                           ! only Fe-ions are efficient to heat grains.
   REAL(KIND=dp), PARAMETER :: Yphdes                      =  1.0e-4_dp    ! Photodesorption yield
   
   ! - Computed by the model
   REAL(KIND=dp)            :: grain_radius                =  1.0e-5_dp
   REAL(KIND=dp)            :: gtodn                                       !  Gas to dust number ratio. 1/gtodn is equivalent to the grain abundance [no unit]
   REAL(KIND=dp)            :: nb_sites_per_grain
   REAL(KIND=dp)            :: dtg_mass_ratio = 1.0e-2_dp
   REAL(KIND=dp)            :: ndsigma, nsites
   ! - 3-phase model
   REAL(KIND=dp), PARAMETER :: nb_active_lay = 2.0_dp
   real(KIND=dp)            :: rate_tot_acc,rate_tot_des, rate_tot
   real(KIND=dp)            :: Nlaysurfsave,Nlaymantsave
   !----------------------------------------------------------
      
   !----------------------------------------------------------
   ! COLUMN DENSITIES USED TO COMPUTE SELF SHIELDING
   !----------------------------------------------------------
   REAL(KIND=dp)            :: Ntot = 1.0e-99_dp
   REAL(KIND=dp)            :: NCO = 1.0e-99_dp
   REAL(KIND=dp)            :: N13CO = 1.0e-99_dp
   REAL(KIND=dp)            :: NH2 = 1.0e-99_dp
   REAL(KIND=dp)            :: NH2v = 1.0e-99_dp
   REAL(KIND=dp)            :: NN2 = 1.0e-99_dp
   REAL(KIND=dp)            :: NH = 1.0e-99_dp
   REAL(KIND=dp)            :: NCO_up = 1.0e-99_dp
   REAL(KIND=dp)            :: N13CO_up = 1.0e-99_dp
   REAL(KIND=dp)            :: NH2_up = 1.0e-99_dp
   REAL(KIND=dp)            :: NH2v_up = 1.0e-99_dp
   REAL(KIND=dp)            :: NN2_up = 1.0e-99_dp
   REAL(KIND=dp)            :: NH_up = 1.0e-99_dp
   REAL(KIND=dp)            :: NCO_dn = 1.0e-99_dp
   REAL(KIND=dp)            :: N13CO_dn = 1.0e-99_dp
   REAL(KIND=dp)            :: NH2_dn = 1.0e-99_dp
   REAL(KIND=dp)            :: NH2v_dn = 1.0e-99_dp
   REAL(KIND=dp)            :: NN2_dn = 1.0e-99_dp
   REAL(KIND=dp)            :: NH_dn = 1.0e-99_dp
   !----------------------------------------------------------
   
   !----------------------------------------------------------
   ! SELF-SHIELDING TABLES
   !----------------------------------------------------------
   REAL(KIND=dp), DIMENSION(251,226) :: NH2CO_file, NCO_file, thetaCO_file
   REAL(KIND=dp), DIMENSION(226,226) :: NH2N2_file, NN2_file, thetaN2_file
         
   !----------------------------------------------------------
   ! USEFUL INDEXES
   !----------------------------------------------------------
   CHARACTER(len=11),PARAMETER :: speelec     = 'e-         ' !  Gas phase electrons
   CHARACTER(len=11),PARAMETER :: spegr0      = 'GRAIN0     ' !  Grain zero
   CHARACTER(len=11),PARAMETER :: spegrm      = 'GRAIN-     ' !  Grain minus
   CHARACTER(len=11),PARAMETER :: speh        = 'H          ' !  Gas phase Hydrogen
   CHARACTER(len=11),PARAMETER :: spehp       = 'H+         ' !  Gas phase Hydrogen
   CHARACTER(len=11),PARAMETER :: spejh       = 'JH         ' !  Hydrogen on grains
   CHARACTER(len=11),PARAMETER :: spekh       = 'KH         ' !  Hydrogen on grains
   CHARACTER(len=11),PARAMETER :: speh2       = 'H2         ' !  Gas phase Dihydrogen
   CHARACTER(len=11),PARAMETER :: speh2v      = 'H2v        ' !  Gas phase Dihydrogen (excited)
   CHARACTER(len=11),PARAMETER :: spejh2      = 'JH2        ' !  Dihydrogen on grains
   CHARACTER(len=11),PARAMETER :: spekh2      = 'KH2        ' !  Dihydrogen on grains
   CHARACTER(len=11),PARAMETER :: spehe       = 'He         ' !  Helium
   CHARACTER(len=11),PARAMETER :: spehep      = 'He+        ' !  Gas phase Helium+
   CHARACTER(len=11),PARAMETER :: speco       = 'CO         ' !  Gas phase CO
   CHARACTER(len=11),PARAMETER :: spen2       = 'N2         ' !  Gas phase CO
   CHARACTER(len=11),PARAMETER :: speh2o      = 'H2O        ' !  Gas phase CO
   CHARACTER(len=11),PARAMETER :: speo        = 'O          ' !  Gas phase CO
   INTEGER                     :: indel                    !  Index corresponding to e- in nb_species length arrays
   INTEGER                     :: indh
   INTEGER                     :: indhp
   INTEGER                     :: indh2
   INTEGER                     :: indh2v
   INTEGER                     :: indco
   INTEGER                     :: ind13co
   INTEGER                     :: indhe
   INTEGER                     :: indhep
   INTEGER                     :: indn2
   INTEGER                     :: indgr0
   INTEGER                     :: indgrm
   INTEGER                     :: indh2o
   INTEGER                     :: indo
   INTEGER                     :: indcp
   
   INTEGER                     :: indpah
   INTEGER                     :: indpahm
   INTEGER                     :: indpahp
   INTEGER                     :: indpahcp
   INTEGER                     :: indpahhp
   !----------------------------------------------------------
   
   INTEGER, DIMENSION(0:99) :: type_id_start     !  list of id start for each reaction type given their type number
   INTEGER, DIMENSION(0:99) :: type_id_stop      !  list of id stop for each reaction type given their type number
   
   ! About the optimization so that, for each species, we only check the reactions we know the species is involved.
   INTEGER                              :: max_reactions_same_species
   INTEGER, DIMENSION(:,:), ALLOCATABLE :: relevant_reactions
   INTEGER, DIMENSION(:), ALLOCATABLE   :: nb_reactions_using_species
   
   !----------------------------------------------------------
   TYPE element
      CHARACTER(LEN=11)                           :: name
      REAL(KIND=dp)                               :: ab
      INTEGER                                     :: index
      REAL(KIND=dp)                               :: mass
   END TYPE element
   !----------------------------------------------------------
   
   !----------------------------------------------------------
   TYPE chem_species
      CHARACTER(LEN=11)                           :: name
      REAL(KIND=dp)                               :: ab
      INTEGER                                     :: index
      REAL(KIND=dp)                               :: mass
      INTEGER, DIMENSION(nelem)                   :: compo
      REAL(KIND=dp)                               :: form_enthalpy
      REAL(KIND=dp)                               :: bind
      REAL(KIND=dp)                               :: diff
      INTEGER                                     :: charge
      REAL(KIND=dp)                               :: stick
   END TYPE chem_species
   !----------------------------------------------------------
   
   !----------------------------------------------------------
   TYPE reactions
      CHARACTER(LEN=11), DIMENSION(ncomp)         :: comp
      INTEGER, DIMENSION(3)                       :: r
      INTEGER, DIMENSION(5)                       :: p
      INTEGER                                     :: id
      INTEGER                                     :: itype
      REAL(KIND=dp)                               :: tmin
      REAL(KIND=dp)                               :: tmax
      REAL(KIND=dp)                               :: alpha
      REAL(KIND=dp)                               :: beta
      REAL(KIND=dp)                               :: gamma
      REAL(KIND=dp)                               :: tmp
      INTEGER                                     :: formula
      REAL(KIND=dp)                               :: rate
      REAL(KIND=dp)                               :: br
      REAL(KIND=dp)                               :: br1
      REAL(KIND=dp)                               :: br2
   END TYPE reactions
   !----------------------------------------------------------
   
   !----------------------------------------------------------
   TYPE(chem_species), DIMENSION(:), ALLOCATABLE  :: spe
   TYPE(element), DIMENSION(:), ALLOCATABLE       :: elem
   TYPE(reactions), DIMENSION(:), ALLOCATABLE     :: reac
   !----------------------------------------------------------
   
   !----------------------------------------------------------
   ! DISK STRUCTURE
   !----------------------------------------------------------
   INTEGER  :: nr
   INTEGER  :: nz
   REAL(KIND=dp) :: G0disk
   INTEGER, DIMENSION(:), ALLOCATABLE             :: nztmp
   REAL(KIND=dp),DIMENSION(:), ALLOCATABLE        :: rdisk
   REAL(KIND=dp),DIMENSION(:,:), ALLOCATABLE      :: zdisk
   REAL(KIND=dp),DIMENSION(:,:), ALLOCATABLE      :: ndisk
   REAL(KIND=dp),DIMENSION(:,:), ALLOCATABLE      :: Avdisk
   REAL(KIND=dp),DIMENSION(:,:), ALLOCATABLE      :: Avupdisk
   REAL(KIND=dp),DIMENSION(:,:), ALLOCATABLE      :: Avdndisk
   REAL(KIND=dp),DIMENSION(:,:), ALLOCATABLE      :: Tgdisk
   REAL(KIND=dp),DIMENSION(:,:), ALLOCATABLE      :: Tddisk
   REAL(KIND=dp),DIMENSION(:,:,:), ALLOCATABLE    :: ColDensdisk
   REAL(KIND=dp),DIMENSION(:,:,:), ALLOCATABLE    :: ColDensupdisk
   REAL(KIND=dp),DIMENSION(:,:,:), ALLOCATABLE    :: ColDensdndisk
   !REAL(KIND=dp),DIMENSION(:,:,:), ALLOCATABLE    :: Xray_table
   !REAL(KIND=dp),DIMENSION(:,:,:), ALLOCATABLE    :: pah_rate_table
   !REAL(KIND=dp),DIMENSION(:,:), ALLOCATABLE      :: in_ndust
   !REAL(KIND=dp),DIMENSION(:,:), ALLOCATABLE      :: in_ndusta
   !REAL(KIND=dp),DIMENSION(:,:), ALLOCATABLE      :: in_ndusta2
   REAL(KIND=dp) :: in_ndust
   REAL(KIND=dp) :: in_ndusta
   REAL(KIND=dp) :: in_ndusta2
   ! Outputs
   REAL(KIND=dp),DIMENSION(:,:,:), ALLOCATABLE    :: abdisk
   REAL(KIND=dp),DIMENSION(:,:,:), ALLOCATABLE    :: ratesdisk
   !----------------------------------------------------------
   
   REAL(KIND=dp), dimension(47) :: NCO_visser
   REAL(KIND=dp), dimension(42) :: NH2_visser
   REAL(KIND=dp), dimension(47,42) :: shield_12c16o_visser, shield_12c17o_visser, shield_12c18o_visser
   REAL(KIND=dp), dimension(47,42) :: shield_13c16o_visser, shield_13c17o_visser, shield_13c18o_visser   

END MODULE global_variables
