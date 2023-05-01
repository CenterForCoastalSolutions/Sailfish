import cupy as cp
import mod_param
import mod_ncparam
import mod_scalars


#   Open boundary conditions arrays:
#
#   zetaBC: Free-surface  (m)   boundary conditions. Substitutes old zeta_east, zeta_west...
#   ubarBC: 2D u-momentum (m/s) boundary conditions. Substitutes old ubar_east, ubar_west...
#   vbarBC: 2D v-momentum (m/s) boundary conditions. Substitutes old vbar_east, vbar_west...

#
#
# -----------------------------------------------------------------------
#   Lateral boundary condition apply switches.
# -----------------------------------------------------------------------
#
#   The following switches are used to control which grid points are
#   processed by the lateral boundary conditions. These switches are
#   set to TRUE by default.  However in composite grids, the points
#   processed by nesting are set to FALSE to allow mixed boundary
#   conditions along the grid edges.



# Lateral boundary conditions structure.
# --------------------------------------

class Boundary:
    def __init__(self):

        # self.LBi = BOUNDS.LBi
        # self.UBi = BOUNDS.UBi
        # self.LBj = BOUNDS.LBj
        # self.UBj = BOUNDS.UBj
        #
        # self.Xsize = UBi - LBi
        # self.Ysize = UBj - LBj

        self.countBCNodes = 1000  # Total number of 2D boundary nodes.

        # Lateral boundary conditions apply switches.
        # These switches are part of the application grid and will be set to
        # FALSE elsewhere, if the boundary point is assigned by a nested grid.
        LBC_apply = cp.ones(self.countBCNodes, dtype = cp.Bool)



# WARNING:
#     We have to compute: zetaFSBCIdx, ubarBCIdx, vbarBCIdx, tBCIdx
#####################
    def initialize_boundary():

        IniVal = 0.0

        setBoundaryIdx = (iwest, ieast, inorth, isouth)

        # 2D boundary info.
        # 2D here means that this information has no dependency on Z (or k index).

        BoundaryTypeIdx = cp.zeros(Nbounds, cp.Int32)
        Boundary2DIndex = cp.zeros(Nbounds, cp.Int64)

        # !  Initialize module variables.

        # Acquire means "process lateral boundary data"
        BOUNDARY.zetaFSBC[zetaFSBCIdx] = IniVal     # zetaFSBC: Free-surface (m)  boundary conditions.
        BOUNDARY.ubarBC  [ubarBCIdx  ] = IniVal     # 2D u-momentum (m/s) boundary conditions.
        BOUNDARY.vbarBC  [vbarBCIdx  ] = IniVal     # 2D v-momentum (m/s) boundary conditions
        BOUNDARY.tBC     [tBCIdx     ] = IniVal     # TRACES boundary conditions. The T (and t) in isTvar, tBC and t_west means 'tracers'.




        # Defines aliases, using the original code's naming.
        BOUNDARY.zeta_west, BOUNDARY.zeta_east, BOUNDARY.zeta_north, BOUNDARY.zeta_south = BOUNDARY.zetaFSBC
        BOUNDARY.ubar_west, BOUNDARY.ubar_east, BOUNDARY.ubar_north, BOUNDARY.ubar_south = BOUNDARY.ubarBC
        BOUNDARY.vbar_west, BOUNDARY.vbar_east, BOUNDARY.vbar_north, BOUNDARY.vbar_south = BOUNDARY.vbarBC
        BOUNDARY.t_west,    BOUNDARY.t_east,    BOUNDARY.t_north,    BOUNDARY.t_south    = BOUNDARY.tBC



