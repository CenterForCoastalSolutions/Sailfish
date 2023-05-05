import cupy as cp
import mod_param
import mod_ncparam
import mod_scalars


#   Open boundary conditions arrays:
#
#   zetaBC: Free-surface  (m)   boundary conditions. Substitutes old zeta_east, zeta_west...
#   ubarBC: 2D u-momentum (m/s) boundary conditions. Substitutes old ubar_east, ubar_west...
#   vbarBC: 2D v-momentum (m/s) boundary conditions. Substitutes old vbar_east, vbar_west...


# Lateral boundary conditions structure.
# --------------------------------------

class Boundary:
    def __init__(self):

        self.zetaFSBCIdx = decodeLBC(getInputVal('LBC(isFsur)'))
        self.ubarBCIdx   = decodeLBC(getInputVal('LBC(isUbar)'))
        self.vbarBCIdx   = decodeLBC(getInputVal('LBC(isVbar)'))

        self.VolCons     = decodeVolCon(west  = getInputVal('VolCons(west)'), east  = getInputVal('VolCons(east)'),
                                        north = getInputVal('VolCons(west)'), south = getInputVal('VolCons(west)'))


        self.nBCfiles    = getInputVal('NBCFILES', dtype = int, minVal = 0)
        self.max_Ffiles  = self.nBCfiles


        self.BRYids      = cp.zeros(self.max_Ffiles) - 1
        self.NBCcount    = cp.zeros(self.max_Ffiles)


        XXXXX # REWRITE IN PYTHON
        line = getInputVal('NBCFILES')
        label='BRY - lateral open boundary conditions'
        BRY =load_s2d(Nval, Cval, Cdim, line, label, ibcfile, igrid, nBCfiles, NBCcount, max_Ffiles)





        self.countBCNodes = 1000  # Total number of 2D boundary nodes.

        # Lateral boundary conditions apply switches.
        # These switches are part of the application grid and will be set to
        # FALSE elsewhere, if the boundary point is assigned by a nested grid.
        LBC_apply = cp.ones(self.countBCNodes, dtype = cp.Bool)

# WARNING:
#     We have to compute: zetaFSBCIdx, ubarBCIdx, vbarBCIdx, tBCIdx
#####################
#####################
#####################


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







