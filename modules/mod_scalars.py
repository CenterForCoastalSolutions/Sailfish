import mod_param

import dateclock



# Time clock structure.
# See ref_clock() in dateclock.py for a description of the possible values of time_ref

time_ref = 0.0                    # YYYYMMDD.dd



Rclock = dateclock.Clock()       # reference/base date







# Control switches.
# -----------------------------------------------------------------------


# If equilibrium tides, switch to apply the 18.6-year lunar nodal cycle correction.
Lnodal = True





# Various model parameters.  Some of these parameters are overwritten
# with the values provided from model standard input script.
# -----------------------------------------------------------------------

spherical  = False       # Switch for spherical grid (lon,lat) configurations.

Lstiffness = True   # Switch to compute the grid stiffness.



VolCons = None   # Lateral open boundary edges volume conservation switches.



#
# !
# !  Switches to read and process climatology fields.
# !
#         logical, allocatable :: CLM_FILE(:)          ! Process NetCDF
#         logical, allocatable :: Lclimatology(:)      ! any field
#         logical, allocatable :: LsshCLM(:)           ! free-surface
#         logical, allocatable :: Lm2CLM(:)            ! 2D momentum
#         logical, allocatable :: Lm3CLM(:)            ! 3D momentum
#         logical, allocatable :: LtracerCLM(:,:)      ! tracers
#
#
# !  Switches to activate point Source/Sinks in an application:
# !    * Horizontal momentum transport (u or v)
# !    * Vertical mass transport (w)
# !    * Tracer transport
# !
#         logical, allocatable :: LuvSrc(:)            ! momentum
#         logical, allocatable :: LwSrc(:)             ! mass
#         logical, allocatable :: LtracerSrc(:,:)      ! tracers
#
#
#         blowup = 0
#         NoError = 0
#         blowup_string = ''
#
#         # thresholds maximum speed (m/s) and density anomaly (kg/m3) to test if the model is blowing-up.
#         max_speed = 20.0        # m/s
#         max_rho   = 200.0       # kg/m3
#
#
#         # Interpolation scheme.
#         linear, cubic = (0, 1)
#         InterpFlag = linear
#
#
#         # Shallowest and Deepest levels to apply bottom momentum stresses as a bodyforce
#         levsfrc = 0.0
#         levbfrc = 0.0
#
#
#
# !  Vertical coordinates transform.  Currently, there are two vertical
# !  transformation equations (see set_scoord.F for details):
# !
# !    Original transform (Vtransform=1):
# !
# !         z_r(x,y,s,t) = Zo_r + zeta(x,y,t) * [1.0 + Zo_r / h(x,y)]
# !
# !                 Zo_r = hc * [s(k) - C(k)] + C(k) * h(x,y)
# !
# !    New transform (Vtransform=2):
# !
# !         z_r(x,y,s,t) = zeta(x,y,t) + [zeta(x,y,t)+ h(x,y)] * Zo_r
# !
# !                 Zo_r = [hc * s(k) + C(k) * h(x,y)] / [hc + h(x,y)]
# !
#         integer, allocatable :: Vtransform(:)
# !
# !  Vertical grid stretching function flag:
# !
# !    Vstretcing = 1   Original function (Song and Haidvogel, 1994)
# !               = 2   A. Shchepetkin (ROMS-UCLA) function
# !               = 3   R. Geyer BBL function
# !
#         integer, allocatable :: Vstretching(:)
# !
# !  Vertical grid stretching parameters.
# !
# !    Tcline        Width (m) of surface or bottom boundary layer in
# !                    which higher vertical resolution is required
# !                    during stretching.
# !    hc            S-coordinate critical depth, hc=MIN(hmin,Tcline).
# !    theta_s       S-coordinate surface control parameter.
# !    theta_b       S-coordinate bottom control parameter.
# !
#         real(dp), allocatable :: Tcline(:)      ! m, positive
#         real(dp), allocatable :: hc(:)          ! m, positive
#         real(dp), allocatable :: theta_s(:)     ! 0 < theta_s < 20
#         real(dp), allocatable :: theta_b(:)     ! 0 < theta_b < 1
# !
# !  Bathymetry range values.
# !
#         real(dp), allocatable :: hmin(:)        ! m, positive
#         real(dp), allocatable :: hmax(:)        ! m, positive
# !
# !  Length (m) of domain box in the XI- and ETA-directions.
# !
#         real(r8), allocatable :: xl(:)          ! m
#         real(r8), allocatable :: el(:)          ! m
# !
# !  Minimum and Maximum longitude and latitude at RHO-points
# !
#         real(r8), allocatable :: LonMin(:)      ! degrees east
#         real(r8), allocatable :: LonMax(:)      ! degrees east
#         real(r8), allocatable :: LatMin(:)      ! degrees north
#         real(r8), allocatable :: LatMax(:)      ! degrees north
#
#
#
# !  Minimun and maximum grid spacing
# !
#         real(dp), allocatable :: DXmin(:)       ! all grid points
#         real(dp), allocatable :: DXmax(:)
#         real(dp), allocatable :: DYmin(:)
#         real(dp), allocatable :: DYmax(:)
#         real(dp), allocatable :: DXminW(:)      ! only water points
#         real(dp), allocatable :: DXmaxW(:)
#         real(dp), allocatable :: DYminW(:)
#         real(dp), allocatable :: DYmaxW(:)
#         real(dp), allocatable :: DZminW(:)      ! only water points
#         real(dp), allocatable :: DZmaxW(:)
#
#
# !  Maximum size of a grid node (m) over the whole curvilinear grid
# !  application. Used for scaling horizontal mixing by the grid size.
# !
#         real(dp), allocatable :: grdmax(:)
#
# !
# !  Courant Numbers due to gravity wave speed limits.
# !
#         real(dp), allocatable :: Cg_min(:)      ! Minimun barotropic
#         real(dp), allocatable :: Cg_max(:)      ! Maximun barotropic
#         real(dp), allocatable :: Cg_Cor(:)      ! Maximun Coriolis
# !
# !  Time dependent Courant Numbers due to velocity components and
# !  indices location of maximum value.
# !
#         integer :: max_Ci = 0                   ! maximum I-location
#         integer :: max_Cj = 0                   ! maximum J-location
#         integer :: max_Ck = 0                   ! maximum K-location
#         real(r8) :: max_C = 0.0_r8              ! maximum total
#         real(r8) :: max_Cu = 0.0_r8             ! maximum I-component
#         real(r8) :: max_Cv = 0.0_r8             ! maximum J-component
#
#
# !  Linear equation of state parameters.
# !
# !    R0            Background constant density anomaly (kg/m3).
# !    Tcoef         Thermal expansion coefficient (1/Celsius).
# !    Scoef         Saline contraction coefficient (1/PSU).
# !
#         real(r8), allocatable :: R0(:)
#         real(r8), allocatable :: Tcoef(:)
#         real(r8), allocatable :: Scoef(:)
# !
# !  Background potential temperature (Celsius) and salinity (PSU) values
# !  used in analytical initializations.
# !
#         real(r8), allocatable :: T0(:)
#         real(r8), allocatable :: S0(:)
# !
# !  Slipperiness variable, either 1.0 (free slip) or -1.0 (no slip).
# !
#         real(r8), allocatable :: gamma2(:)
# !
# !  Weighting coefficient for the newest (implicit) time step derivatives
# !  using either a Crack-Nicolson implicit scheme (lambda=0.5) or a
# !  backward implicit scheme (lambda=1.0).
# !
#         real(r8) :: lambda = 1.0_r8
# !
# !  Jerlov water type to assign everywhere, range values: 1 - 5.
# !
#         integer, allocatable :: lmd_Jwt(:)
# !
# !  Grid r-factor (non-dimensional).
# !
#         real(dp), allocatable :: rx0(:)         ! Beckmann and Haidvogel
#         real(dp), allocatable :: rx1(:)         ! Haney
# !
# !  Linear (m/s) and quadratic (nondimensional) bottom drag coefficients.
# !
#         real(r8), allocatable :: rdrg(:)
#         real(r8), allocatable :: rdrg2(:)
# !
# !  Minimum and maximum threshold for transfer coefficient of momentum.
# !
#         real(dp) :: Cdb_min = 0.000001_dp
#         real(dp) :: Cdb_max = 0.5_dp
# !
# !  Surface and bottom roughness (m)
# !
#         real(r8), allocatable :: Zos(:)
#         real(r8), allocatable :: Zob(:)
# !
# !  Minimum depth for wetting and drying (m).
# !
#         real(r8), allocatable :: Dcrit(:)
# !
# !  Mean density (Kg/m3) used when the Boussinesq approximation is
# !  inferred.
# !
#         real(dp) :: rho0 = 1025.0_dp
# !
# !  Background Brunt-Vaisala frequency (1/s2).
# !
#         real(dp) :: bvf_bak = 0.00001_dp
#
# !  Vector containing USER generic parameters.
# !
#         integer :: Nuser
#         real(r8), dimension(25) :: user(25)
# !
# !  Weights for the time average of 2D fields.
# !
#         real(dp), allocatable :: weight(:,:,:)
# !
# !  Constants.
# !
#         pi = cp.pi
#         deg2rad = pi/180.0
#         rad2deg = 180.0/pi
#         day2sec = 86400.0
#         sec2day = 1.0/86400.0
#         spval   = 1.0E+37      # "special" value
#         Large   = 1.0E+20
#         jul_off = 2440000.0
# !
# !  Set special check value.  Notice that a smaller value is assigned
# !  to account for both NetCDF fill value and roundoff. There are
# !  many Matlab scripts out there that do not inquire correctly
# !  the spval from the _FillValue attribute in single/double
# !  precision.
# !
#         spval_check = 1.0E+35
# !
# !-----------------------------------------------------------------------
# !  Horizontal and vertical constant mixing coefficients.
# !-----------------------------------------------------------------------
# !
# !    Akk_bak       Background vertical mixing coefficient (m2/s) for turbulent energy.
# !    Akp_bak       Background vertical mixing coefficient (m2/s) for generic statistical field "psi".
# !    Akt_bak       Background vertical mixing coefficient (m2/s) for tracers.
# !    Akv_bak       Background vertical mixing coefficient (m2/s) for momentum.
# !    Akt_limit     Upper threshold vertical mixing coefficient (m2/s) for tracers.
# !    Akv_limit     Upper threshold vertical mixing coefficient (m2/s) for momentum.
# !    Kdiff         Isopycnal mixing thickness diffusivity (m2/s) for tracers.
# !    ad_visc2      ADM lateral harmonic constant mixing coefficient (m2/s) for momentum.
# !    nl_visc2      NLM lateral harmonic constant mixing coefficient (m2/s) for momentum.
# !    tl_visc2      TLM lateral harmonic constant mixing coefficient (m2/s) for momentum.
# !    visc2         Current lateral harmonic constant mixing coefficient (m2/s) for momentum.
# !    ad_visc4      ADM lateral biharmonic (squared root) constant mixing coefficient (m2 s^-1/2) for momentum.
# !    nl_visc4      NLM lateral biharmonic (squared root) constant mixing coefficient (m2 s^-1/2) for momentum.
# !    tl_visc4      TLM lateral biharmonic (squared root) constant mixing coefficient (m2 s^-1/2) for momentum.
# !    visc4         Current lateral biharmonic (squared root) constant mixing coefficient (m2 s^-1/2) for momentum.
# !    ad_tnu2       ADM lateral harmonic constant mixing coefficient (m2/s) for tracer type variables.
# !    nl_tnu2       NLM lateral harmonic constant mixing coefficient (m2/s) for tracer type variables.
# !    tl_tnu2       TLM lateral harmonic constant mixing coefficient (m2/s) for tracer type variables.
# !    tnu2          Current lateral harmonic constant mixing coefficient (m2/s) for tracer type variables.
# !    ad_tnu4       ADM lateral biharmonic (squared root) constant mixing coefficient (m2 s^-1/2) for tracers.
# !    nl_tnu4       NLM lateral biharmonic (squared root) constant mixing coefficient (m2 s^-1/2) for tracers.
# !    tl_tnu4       TLM lateral biharmonic (squared root) constant mixing coefficient (m2 s^-1/2) for tracers.
# !    tnu4          Current lateral biharmonic (squared root) constant mixing coefficient (m2 s^-1/2) for tracers.
# !    tkenu2        Lateral harmonic constant mixing coefficient (m2/s) for turbulent energy.
# !    tkenu4        Lateral biharmonic (squared root) constant mixing coefficient (m2 s^-1/2) for turbulent energy.
# !

#
#
#         # Switches to increase/decrease horizontal viscosity and/or diffusion
#         # in specific areas of the application domain (like sponge areas).
#         logical, allocatable :: Lsponge(:)
#         logical, allocatable :: LuvSponge(:)         ! viscosity
#         logical, allocatable :: LtracerSponge(:,:)   ! diffusion
#
#
#         # IO parameters.
#         # -----------------------------------------------------------------------
#
#         # Switches to activate creation and writing of output NetCDF files.
#         logical, allocatable :: LdefAVG(:)       ! Average file
#         logical, allocatable :: LdefHIS(:)       ! History file
#         logical, allocatable :: LdefHSS(:)       ! Hessian file
#         logical, allocatable :: LdefINI(:)       ! Initial file
#         logical, allocatable :: LdefIRP(:)       ! Initial RPM file
#         logical, allocatable :: LdefITL(:)       ! Initial TLM file
#         logical, allocatable :: LdefLCZ(:)       ! Lanczos Vectors file
#         logical, allocatable :: LdefLZE(:)       ! Evolved Lanczos file
#         logical, allocatable :: LdefQCK(:)       ! Quicksave file
#         logical, allocatable :: LdefRST(:)       ! Restart file
#
#         logical, allocatable :: LdefSTA(:)       ! Stations file
#         logical, allocatable :: LdefTIDE(:)      ! tide forcing file
#         logical, allocatable :: LdefTLF(:)       ! TLM/RPM impulse file
#
#         logical, allocatable :: LreadADM(:)      ! Read ADM multi-files
#         logical, allocatable :: LreadBLK(:)      ! Read NLM bulk fluxes
#         logical, allocatable :: LreadFRC(:)      ! Read FRC files
#         logical, allocatable :: LreadFWD(:)      ! Read FWD trajectory
#         logical, allocatable :: LreadQCK(:)      ! Read QCK trajectory
#         logical, allocatable :: LreadTLM(:)      ! Read TLM multi-files
#
#         logical, allocatable :: LwrtAVG(:)       ! Write average file
#         logical, allocatable :: LwrtHIS(:)       ! Write history file
#         logical, allocatable :: LwrtPER(:)       ! Write during ensemble
#         logical, allocatable :: LwrtQCK(:)       ! write quicksave file
#         logical, allocatable :: LwrtRST(:)       ! Write restart file
#         logical, allocatable :: LwrtTLF(:)       ! Write impulse file
#
#         logical, allocatable :: LdefNRM(:,:)     ! Norm file
#         logical, allocatable :: LwrtNRM(:,:)     ! Write norm file
#
# !
# !  Switch to append information to an existing ROMS standard output
# !  log file.
# !
#         logical :: Lappend = .FALSE.
# !
# !  Switch to read input open boundary conditions data.
# !
#         logical, allocatable :: LprocessOBC(:)
# !
# !  Switch to read input tidal forcing data.
# !
#         logical, allocatable :: LprocessTides(:)
# !
# !  Switch to write application set-up information to standard output.
# !
#         logical, allocatable :: LwrtInfo(:)
# !
# !  Switch used to create new output NetCDF files. If TRUE, new output
# !  files are created. If FALSE, data is appended to an existing output
# !  files.  Used only for history, average and station files.
# !
#         logical, allocatable :: ldefout(:)       ! New output files
# !
# !  Number of timesteps between creation of new output files.
# !
#         integer, allocatable :: ndefHIS(:)       ! History file
#         integer, allocatable :: ndefQCK(:)       ! Quicksave file
#
# !
# !  Number of timesteps between writing of output data.
# !
#         integer, allocatable :: nHIS(:)          ! History file
#         integer, allocatable :: nQCK(:)          ! Quicksave file
#         integer, allocatable :: nRST(:)          ! Restart file
#         integer, allocatable :: nSTA(:)          ! Stations file
#
#
# !
# !  Number of timesteps between print of single line information to
# !  standard output.
# !
#         integer, allocatable :: ninfo(:)
# !
# !  Here, it is assumed that nOBC is a multiple of NTIMES or greater
# !  than NTIMES. If nOBC > NTIMES, only one record is stored in the
# !  output history NetCDF files and the adjustment is for constant
# !  open boundaries with constant correction.
# !
#         integer, allocatable :: nOBC(:)          ! number of timesteps
#         integer, allocatable :: Nbrec(:)         ! number of records
#         integer, allocatable :: OBCcount(:)      ! record counter
#
#
# !
# !  Restart time record to read from disk and use as the initial
# !  conditions. Use nrrec=0 for new solutions. If nrrec is negative
# !  (say, nrrec=-1), the model will restart from the most recent
# !  time record. That is, the initialization record is assigned
# !  internally.
# !
#         integer, allocatable :: nrrec(:)
# !
# !  Switch to activate processing of input data.  This switch becomes
# !  very useful when reading input data serially in parallel
# !  applications.
# !
#         logical, allocatable :: synchro_flag(:)
# !
# !  Switch to inialize model with latest time record from initial
# !  (restart/history) NetCDF file.
# !
#         logical, allocatable :: LastRec(:)
#
# !
# !  Switches used to recycle time records in some output file. If TRUE,
# !  only the latest two time records are maintained.  If FALSE, all
# !  field records are saved.
# !

#
#
#         real(r8) :: lmd_nuwm = 1.0E-5_r8      ! m2/s
#         real(r8) :: lmd_nuws = 1.0E-6_r8      ! m2/s
#
#         real(r8) :: lmd_sdd1 = 0.15_r8        ! non-dimensional
#         real(r8) :: lmd_sdd2 = 1.85_r8        ! non-dimensional
#         real(r8) :: lmd_sdd3 = 0.85_r8        ! non-dimensional
#         real(r8) :: lmd_tdd1 = 0.909_r8       ! non-dimensional
#         real(r8) :: lmd_tdd2 = 4.6_r8         ! non-dimensional
#         real(r8) :: lmd_tdd3 = 0.54_r8        ! non-dimensional
#
# !
# !-----------------------------------------------------------------------
# !  Generic Length Scale parameters.

#
# !
# ! Constants used in the various formulation of surface flux boundary
# ! conditions for the GLS vertical turbulence closure in terms of
# ! Charnok surface roughness (CHARNOK_ALPHA), roughness from wave
# ! amplitude (zos_hsig_alpha), wave dissipation (SZ_ALPHA), and
# ! Craig and Banner wave breaking (CRGBAN_CW).
# ! Wec_alpha partitions energy to roller or breaking.
#         real(r8), allocatable :: charnok_alpha(:)
#         real(r8), allocatable :: zos_hsig_alpha(:)
#         real(r8), allocatable :: sz_alpha(:)
#         real(r8), allocatable :: crgban_cw(:)
#         real(r8), allocatable :: wec_alpha(:)
#
#
#
#       SUBROUTINE initialize_scalars
#
# !=======================================================================
# !                                                                      !
# !  This routine initializes several variables in module for all nested !
# !  grids.                                                              !
# !                                                                      !
# !=======================================================================
#
#         IniVal = 0.0
#
#
#
#         # Tracer identification indices.
#         # ---------------------------------------------------------------------
#         itemp, isalt = (1,2)
#
#
#
#
#
#         # Activate all computation control switches.
#         # -----------------------------------------------------------------------
#         CompositeGrid[0:3] = False
#         LastRec      = False
#         RefinedGrid  = False
#         GetDonorData = False
#         Lbiology     = True
#         Lfloats      = True
#         Lsediment    = True
#         Lstations    = True
#
#
#
#         # Initialize several scalar variables.
#         # -----------------------------------------------------------------------
#         Co = 1.0/(2.0 + cp.sqrt(2.0))
#         gorho0 = g/rho0

#         SetGridConfig = True
#         RefineScale = 0
#         gamma2 = -1.0
#         Vtransform = 1
#         Vstretching = 1
#
#         first_time = 0
#
#         maxspeed = -Large
#         maxrho   =-Large
#         TotVolume = 0.0
#         MinVolume =  Large
#         MaxVolume = -Large
#         DXmin =  Large
#         DXmax = -Large
#         DYmin =  Large
#         DYmax = -Large
#
#         DXminW =  Large
#         DXmaxW = -Large
#         DYminW =  Large
#         DYmaxW = -Large
#
#         DZminW =  Large
#         DZmaxW = -Large
#
#         grdmax = -Large
#
#         Cg_min =  Large
#         Cg_max = -Large
#         Cg_Cor = -Large
#
#         rx0 = -Large
#         rx1 = -Large
#         CLM_FILE = .FALSE.
#         Lnudging = .FALSE.
#         LnudgeM2CLM = .FALSE.
#         LnudgeM3CLM = .FALSE.
#         Lclimatology = .FALSE.
#         Lm2CLM = .FALSE.
#         Lm3CLM = .FALSE.
#         LsshCLM = .FALSE.
#         Lsponge = .FALSE.
#         LuvSponge = .FALSE.
#         LuvSrc = .FALSE.
#         LwSrc = .FALSE.
#
#
#         # Tracers
#         DO itrc=1,MT
#           LnudgeTCLM(itrc,ng)=.FALSE.
#           LtracerCLM(itrc,ng)=.FALSE.
#           LtracerSrc(itrc,ng)=.FALSE.
#           LtracerSponge(itrc,ng)=.FALSE.
#           ad_Akt_fac(itrc,ng)=1.0_r8
#           tl_Akt_fac(itrc,ng)=1.0_r8
#           ad_tnu2(itrc,ng)=IniVal
#           nl_tnu2(itrc,ng)=IniVal
#           tl_tnu2(itrc,ng)=IniVal
#           tnu2(itrc,ng)=IniVal
#           ad_tnu4(itrc,ng)=IniVal
#           nl_tnu4(itrc,ng)=IniVal
#           tl_tnu4(itrc,ng)=IniVal
#           tnu4(itrc,ng)=IniVal
#         END DO
#
#         DO itrc=1,NAT
#           Akt_limit(itrc,ng)=1.0E-3_r8
#         END DO
#
#         Akv_limit(ng)=1.0E-3_r8
#         ad_Akv_fac(ng)=1.0_r8
#         tl_Akv_fac(ng)=1.0_r8
#         ad_visc2(ng)=IniVal
#         nl_visc2(ng)=IniVal
#         tl_visc2(ng)=IniVal
#         visc2(ng)=IniVal
#         ad_visc4(ng)=IniVal
#         nl_visc4(ng)=IniVal
#         tl_visc4(ng)=IniVal
#         visc4(ng)=IniVal
#
#
#         DO i=1,4
#           VolCons(i,ng)=.FALSE.
#           FSobc_in (ng,i)=0.0_dp
#           FSobc_out(ng,i)=0.0_dp
#           M2obc_in (ng,i)=0.0_dp
#           M2obc_out(ng,i)=0.0_dp
#
#         END DO
#       END DO
#
#

#
#         # Initialize thread private variables.
#         synchro_flag=.FALSE.
#
#         ntfirst=1
#         ntstart=1
#         ntend=0
#         step_counter=0
#
#
# !
# !  Initialize several IO flags.
# !
#       LmultiGST=.FALSE.
#       LrstGST=.FALSE.
#
#       DO ng=1,Ngrids
#
#         PerfectRST(ng)=.FALSE.
#
#         Ladjusted(ng)=.FALSE.
#
#         LprocessOBC(ng)=.FALSE.
#         LprocessTides(ng)=.FALSE.
#
#         LdefADJ(ng)=.FALSE.
#         LdefAVG(ng)=.TRUE.
#         LdefDAI(ng)=.FALSE.
#         LdefDIA(ng)=.TRUE.
#         LdefERR(ng)=.FALSE.
#         LdefHIS(ng)=.TRUE.
#         LdefINI(ng)=.FALSE.
#
#         LdefIRP(ng)=.FALSE.
#         LdefITL(ng)=.FALSE.
#         LdefMOD(ng)=.FALSE.
#         LdefQCK(ng)=.FALSE.
#         LdefRST(ng)=.TRUE.
#         LdefSTA(ng)=.TRUE.
#         LdefTLM(ng)=.FALSE.
#
#         LdefTIDE(ng)=.FALSE.
#
#         LreadADM(ng)=.FALSE.
#         LreadBLK(ng)=.FALSE.
#         LreadFRC(ng)=.FALSE.
#
#         LreadFWD(ng)=.FALSE.
#         LreadQCK(ng)=.FALSE.
#         LreadTLM(ng)=.FALSE.
#         LwrtADJ(ng)=.FALSE.
#         LwrtAVG(ng)=.FALSE.
#         LwrtDIA(ng)=.FALSE.
#         LwrtHIS(ng)=.FALSE.
#         LwrtPER(ng)=.FALSE.
#         LwrtQCK(ng)=.FALSE.
#         LwrtRST(ng)=.FALSE.
#         LwrtTLM(ng)=.FALSE.
#         LwrtInfo(ng)=.TRUE.
#         LwrtState2d(ng)=.FALSE.
#         LwrtTime(ng)=.TRUE.
#         LwrtCost(ng)=.FALSE.
#         ldefout(ng)=.FALSE.
#
#       END DO
#
#
#         # Initialize the NLM initial conditions time to a negative number to
#         # check if its value was assigned elsewhere.  It can be used during
#         # the initialization of the adjoint model when DSTART is not the
#         # same as the start of the simulation.
#
#         INItime  = -1.0
#         INItimeS = -1.0
#
