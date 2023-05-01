

# Structure, T_IO, used to store information about the input and output files.
# ----------------------------------------------------------------------------
#  This structure is used to declare the variables associated with
#  input and output files, like TYPE(IO) :: HIS(Ngrids). It is a
#  compact way to store a lot of information.  The "Fcount" variable
#  is used during data processing in "check_multifile" and "inquiry"
#  whereas "load" is used storing information in the structure.


class T_IO:
    def __init__(self):
        pass
        # integer :: Nfiles                        ! number of multi-files
        # integer :: Fcount                        ! multi-file counter
        # integer :: load                          ! filename load counter
        # integer :: Rindex                        ! NetCDF record index
        # integer :: ncid                          ! NetCDF file ID
        # integer,  pointer :: Nrec(:)             ! NetCDF record size
        # integer,  pointer :: Vid(:)              ! NetCDF variables IDs
        # integer,  pointer :: Tid(:)              ! NetCDF tracers IDs
        # real(dp), pointer :: time_min(:)         ! starting time
        # real(dp), pointer :: time_max(:)         ! ending time
        # character (len=50 ) :: label             ! structure label
        # character (len=256) :: head              ! head filename
        # character (len=256) :: base              ! base filename
        # character (len=256) :: name              ! current name
        # character (len=256), pointer :: files(:) ! multi-file names

# !-----------------------------------------------------------------------
# !  Lower and upper bounds indices per domain partition for all grids.
# !-----------------------------------------------------------------------
# !
# !  Notice that these indices have different values in periodic and
# !  nesting applications, and on tiles next to the boundaries. Special
# !  indices are required to process overlap regions (suffices P and T)
# !  lateral boundary conditions (suffices B and M) in nested grid
# !  applications. The halo indices are used in private computations
# !  which include ghost-points and are limited by MAX/MIN functions
# !  on the tiles next to the  model boundaries. For more information
# !  about all these indices, see routine "var_bounds" in file
# !  "Utility/get_bounds.F".
# !
# !  All the 1D array indices are of size -1:NtileI(ng)*NtileJ(ng)-1. The
# !  -1 index include the values for the full (no partitions) grid.
# !
# !  Notice that the starting (Imin, Jmin) and ending (Imax, Jmax) indices
# !  for I/O processing are 3D arrays. The first dimension (1:4) is for
# !  1=PSI, 2=RHO, 3=u, 4=v points; the second dimension (0:1) is number
# !  of ghost points (0: no ghost points, 1: Nghost points), and the
# !  the third dimension is for 0:NtileI(ng)*NtileJ(ng)-1.
# !
#       TYPE T_BOUNDS
#         integer, pointer :: tile(:)  ! tile partition
#
#         integer, pointer :: LBi(:)   ! lower bound I-dimension
#         integer, pointer :: UBi(:)   ! upper bound I-dimension
#         integer, pointer :: LBj(:)   ! lower bound J-dimension
#         integer, pointer :: UBj(:)   ! upper bound J-dimension
#
#         integer :: LBij              ! lower bound MIN(I,J)-dimension
#         integer :: UBij              ! upper bound MAX(I,J)-dimension
#
#         integer :: edge(4,4)         ! boundary edges I- or J-indices
#
#
#
#
#       END TYPE T_BOUNDS

      # TYPE (T_BOUNDS), allocatable :: BOUNDS(:)

#  Lower and upper bounds in NetCDF files.


# !
# !-----------------------------------------------------------------------
# !  Domain boundary edges switches and tiles minimum and maximum
# !  fractional grid coordinates.
# !-----------------------------------------------------------------------
# !
#       TYPE T_DOMAIN
#         real(r8) :: Xmin_psi(:)
#         real(r8) :: Xmax_psi(:)
#         real(r8) :: Ymin_psi(:)
#         real(r8) :: Ymax_psi(:)
#
#         real(r8) :: Xmin_rho(:)
#         real(r8) :: Xmax_rho(:)
#         real(r8) :: Ymin_rho(:)
#         real(r8) :: Ymax_rho(:)
#
#         real(r8) :: Xmin_u(:)
#         real(r8) :: Xmax_u(:)
#         real(r8) :: Ymin_u(:)
#         real(r8) :: Ymax_u(:)
#
#         real(r8) :: Xmin_v(:)
#         real(r8) :: Xmax_v(:)
#         real(r8) :: Ymin_v(:)
#         real(r8) :: Ymax_v(:)

DOMAIN = T_DOMAIN


# Lateral Boundary Conditions (LBC) switches structure.
# -----------------------------------------------------------------------

# The lateral boundary conditions are specified by boolean switches.
# The LBC structure is allocated as:
#
#       LBC[:, nLBCvar]
#
#  where : are the number boundary edges and nLBCvar are the number of state variables.
#  For example, for free-surface gradient boundary conditions we have:
#
#       LBC[idx, isFsur].gradient


class T_LBC:
    def __init__(self):
        self.acquire          = False        # process lateral boundary data
        self.Chapman_explicit = False
        self.Chapman_implicit = False
        self.clamped          = False
        self.closed           = False
        self.Flather          = False
        self.gradient         = False
        self.mixed            = False
        self.nested           = False
        self.nudging          = False
        self.periodic         = False
        self.radiation        = False
        self.reduced          = False
        self.Shchepetkin      = False





def initialize_param():
    """ This routine initializes several parameters in module 'mod_param' """

    # Model grid(s) parameters.
    # -----------------------------------------------------------------------

    #   Number of interior RHO-points in the XI- and ETA-directions. The
    #   size of models state variables (C-grid) at input and output are:
    #     RH0-type variables:  [0:L, 0:M]        ----v(i,j+1)----
    #     PSI-type variables:  [1:L, 1:M]        |              |
    #       U-type variables:  [1:L, 0:M]     u(i,j)  r(i,j)  u(i+1,j)
    #       V-type variables:  [0:L, 1:M]        |              |
    #                                            -----v(i,j)-----
    #  Lm = L - 1
    #  Mm = M - 1

    # Derived dimension parameters.
    # -------------------------------


    # Global horizontal size of model arrays including padding.  All the
    # model state arrays are of same size to facilitate parallelization.

    # Lm, Mm Number of interior grid points in the (ξ,η)-directions
    # Im, Jm Number of global   grid points in the (ξ,η)-directions
    Im = Lm | 1   # Im is the closest odd number larger or equal than Lm
    Jm = Mm | 1   # Jm is the closest odd number larger or equal than Mm (I don't know why an odd number is required)



    # Total number of tracer variables.
    # NAT: Number of active tracers. Usually, NAT=2 for potential temperature and salinity.
    # NBT: Number of biological tracers.
    # NST: Number of sediment tracer type variables (NCS+NNS).
    # NPT: Number of inert passive tracers to advect and diffuse only (like dyes, etc). This parameter is independent of the number of biological and/or sediment tracers
    NT = NAT + NBT + NST + NPT

    # Number of model state variables.
    NSV = NT + 5

    # Set the maximum number of tracers.
    MT = cp.max(NAT, NT)   # This seems unnecessary


    # Allocate Lateral Boundary Conditions switches structure.
    # Each 2D node has an independent T_LBC structure.
    # -----------------------------------------------------------------------
    nLBCvar=3


    for ivar in range(nLBCvar):
        LBC[:, ivar] = T_LBC()




