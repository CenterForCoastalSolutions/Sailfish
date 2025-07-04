import cupy as cp
import os
from misc import *
# import mod_param
# import mod_ncparam
# import mod_scalars




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

setBC = None
copyBC = None

# REMEMBER: This code has to be syncronized with the same definitions in C++ code.
bcNone               = 0
bcAcquire            = 1 << 0       # Read lateral boundary data from files
bcChapman_explicit   = 1 << 1
bcChapman_implicit   = 1 << 2
bcClamped            = 1 << 3
bcClosed             = 1 << 4
bcFlather            = 1 << 5       # (2D momentum)
bcGradient           = 1 << 6
bcMixed              = 1 << 7
bcNested             = 1 << 8
bcNudging            = 1 << 9
bcPeriodic           = 1 << 10
bcRadiation          = 1 << 11
bcReduced            = 1 << 12      # (2D momentum)
bcShchepetkin        = 1 << 13      # (2D momentum)


dictLBC = {
    'CHA'   : bcChapman_explicit,
    'CHE'   : bcChapman_implicit,
    'CLA'   : bcClamped + bcAcquire,
    'CLO'   : bcClosed,
    'FLA'   : bcFlather + bcAcquire,
    'GRA'   : bcGradient,
    'NES'   : bcNested,
    'RAD'   : bcRadiation,
    'RADNUD': bcRadiation + bcNudging + bcAcquire,
    'RED'   : bcReduced,
    'SHC'   : bcShchepetkin,
}

def decodeLBC(val):
    # Converts a LBC string (like 'CHA') into an integer that describes it.
    val = val.upper()

    try:
        return dictLBC[val]
    except:
        msgError('LBC type "%s" is not recognized' % val, 2)


def checkLBC(bc, bcComp):
    # Checks if an LBC codified as an integer contains one "component" (like bcAcquire or bcClosed).
    return (bc & bcComp) != 0




#   Open boundary conditions arrays:
#
#   zetaBC: Free-surface  (m)   boundary conditions. Substitutes old zeta_east, zeta_west...
#   ubarBC: 2D u-momentum (m/s) boundary conditions. Substitutes old ubar_east, ubar_west...
#   vbarBC: 2D v-momentum (m/s) boundary conditions. Substitutes old vbar_east, vbar_west...


# Lateral boundary conditions structure.
# --------------------------------------

iW, iS, iE, iN = (0, 1, 2, 3)


class BCIdx:
    def __init__(self, bcIdx, bcIdx1, bcIdx2, bcIdxFieldIdx2, bcIdxFieldType, LBC):
        self.bcIdx          = bcIdx              # See description of idxFlat below.
        self.bcIdx1         = bcIdx1             # idxFlat1
        self.bcIdx2         = bcIdx2             # idxFlat2
        self.bcIdxFieldIdx2 = bcIdxFieldIdx2     # Same as bcIdx2, but as field (indexed by the node number).
        self.bcIdxFieldType = bcIdxFieldType     # Same as LBC, but as field (indexed by the node number).
        self.LBC            = LBC

    def unpack(self):
        return self.bcIdx, self.bcIdx1, self.bcIdx2, self.LBC


class Boundary:
    def buildLBC(self, val, i0, j0):
        # LBC's are stored in this order: W  S  E  N
        #
        # idxFlat : Index of the *cell* where the BC is being set (corners cells can have two BC nodes)
        # idxFlat1: Index of the node where the BC value is modified (it may be a ghost cell)
        # idxFlat2: index of the node where the BC value is copied (it mau be a fully interior cell).
        #
        # idxFlat is a *cell* index, meaning that U, V and R BC's in the same cell would have the same idxFlat index.
        # It follows the same indexing as R nodes.
        #
        # Note the difference between a cell, indexed by their R node, and a node (with their own indexing). For example,
        # the right U node in cell i has index i+1.

        L = self.GRID.L
        M = self.GRID.M

        bcIdx  = []
        bcIdx1 = []
        bcIdx2 = []
        LBC = []
        shp = (M + 1, L + 1)

        bcW = decodeLBC(val[iW])
        bcE = decodeLBC(val[iE])
        bcN = decodeLBC(val[iN])
        bcS = decodeLBC(val[iS])

        iRange = cp.arange(i0, L)
        jRange = cp.arange(j0, M)

        idxFlat  = cp.ravel_multi_index((jRange, cp.array([i0    ])), shp)
        idxFlat1 = cp.ravel_multi_index((jRange, cp.array([i0    ])), shp)
        idxFlat2 = cp.ravel_multi_index((jRange, cp.array([i0 + 1])), shp)
        bcIdx  += idxFlat .tolist()
        bcIdx1 += idxFlat1.tolist()
        bcIdx2 += idxFlat2.tolist()
        LBC += [bcW]*(M-j0)

        idxFlat  = cp.ravel_multi_index((jRange, cp.array([L - 1])), shp)
        idxFlat1 = cp.ravel_multi_index((jRange, cp.array([L    ])), shp)
        idxFlat2 = cp.ravel_multi_index((jRange, cp.array([L - 1])), shp)
        bcIdx  += idxFlat .tolist()
        bcIdx1 += idxFlat1.tolist()
        bcIdx2 += idxFlat2.tolist()
        LBC += [bcE]*(M-j0)

        idxFlat  = cp.ravel_multi_index((cp.array([i0    ]), iRange), shp)
        idxFlat1 = cp.ravel_multi_index((cp.array([i0    ]), iRange), shp)
        idxFlat2 = cp.ravel_multi_index((cp.array([i0 + 1]), iRange), shp)
        bcIdx  += idxFlat .tolist()
        bcIdx1 += idxFlat1.tolist()
        bcIdx2 += idxFlat2.tolist()
        LBC += [bcN]*(L-i0)

        idxFlat  = cp.ravel_multi_index((cp.array([M - 1]), iRange), shp)
        idxFlat1 = cp.ravel_multi_index((cp.array([M    ]), iRange), shp)
        idxFlat2 = cp.ravel_multi_index((cp.array([M - 1]), iRange), shp)
        bcIdx  += idxFlat .tolist()
        bcIdx1 += idxFlat1.tolist()
        bcIdx2 += idxFlat2.tolist()
        LBC += [bcS]*(L-i0)






        # TODO: Remove???
        # # Corners.
        # idxFlat  = cp.ravel_multi_index((j0,     i0    ), shp)
        # idxFlat1 = cp.ravel_multi_index((j0,     i0    ), shp)
        # idxFlat2 = cp.ravel_multi_index((j0 + 1, i0 + 1), shp)
        # bcIdx  += [idxFlat ]
        # bcIdx1 += [idxFlat1]
        # bcIdx2 += [idxFlat2]
        # LBC += [decodeLBC([val[iN], val[iW]])]


        bcIdx  = cp.array(bcIdx , dtype = cp.int32)
        bcIdx1 = cp.array(bcIdx1, dtype = cp.int32)
        bcIdx2 = cp.array(bcIdx2, dtype = cp.int32)
        LBC    = cp.array(LBC,    dtype = cp.int32)

        # Sorts the indices to improve cache access. TODO: Recover?
        # sortIdx = cp.argsort(bcIdx.flatten())

        # bcIdx  = bcIdx [sortIdx]
        # bcIdx1 = bcIdx1[sortIdx]
        # bcIdx2 = bcIdx2[sortIdx]
        # LBC    = LBC   [sortIdx]

        bcIdxFieldIdx2 = -cp.ones(shp, cp.int32)
        bcIdxFieldType = -cp.ones(shp, cp.int32)
        bcIdxFieldIdx2.ravel()[bcIdx1] = bcIdx2
        bcIdxFieldType.ravel()[bcIdx1] = LBC


        return BCIdx(bcIdx, bcIdx1, bcIdx2, bcIdxFieldIdx2, bcIdxFieldType, LBC)

    def getLocalBCIdx(self, idxBC, flatIdx):
        # Given an array of flat indices and a flat index (that must be in the array), returns the location
        # of that index in the array.
        # Uses binary search.
        flatIdxArray = idxBC.bcIdx
        L = 0
        R = flatIdxArray.size - 1
        M = (R + L)//2

        if flatIdx < flatIdxArray[L]:
            msgError('Index %i not in the list (function Boundary.getBCIdx)' % flatIdx)

        if flatIdx > flatIdxArray[R]:
            msgError('Index %i not in the list (function Boundary.getBCIdx)' % flatIdx)

        if flatIdx == flatIdxArray[L]:
            return L

        if flatIdx == flatIdxArray[R]:
            return R

        while R - L >= 1:

            if flatIdx >= flatIdxArray[M]:
                L = M
            else:
                R = M

            M = (R + L) // 2
            if flatIdx == flatIdxArray[M]:
                return M

        msgError('Index %i not in the list (function Boundary.getBCIdx)' % flatIdx)


    def __init__(self, input, GRID):

        global setBC, copyBC

        filename = os.path.join(exePath, r'modules/mod_BC.c')
        with open(filename, 'r') as file:
            code = file.read()
        moduleBC = cp.RawModule(code=code, options=compilationOptions)
        setBC  = moduleBC.get_function('setBC')
        copyBC = moduleBC.get_function('copyBC')

        self.GRID = GRID

        self.zetaBC = self.buildLBC(input.getVal('LBC(isFsur)').split(), 0, 0)
        self.ubarBC = self.buildLBC(input.getVal('LBC(isUbar)').split(), 1, 0)
        self.vbarBC = self.buildLBC(input.getVal('LBC(isVbar)').split(), 0, 1)
        self.uvelBC = self.buildLBC(input.getVal('LBC(isUvel)').split(), 1, 0)
        self.vvelBC = self.buildLBC(input.getVal('LBC(isVvel)').split(), 0, 1)

        # Performs some adjustments
        # When Flather or Shchepetkin conditions are used for the velocities, we will need to "acquire"
        # the values for the free surfaces too.
        for idx, idx1, idx2, bc in zip(*(self.ubarBC.unpack())):
            if checkLBC(bc, bcFlather) or checkLBC(bc, bcShchepetkin):
                idxLocalZetaBC = self.getLocalBCIdx(self.zetaBC, idx)
                self.zetaBC.LBC[idxLocalZetaBC] |= bcAcquire

        for idx, idx1, idx2, bc in zip(*(self.vbarBC.unpack())):
            if checkLBC(bc, bcFlather) or checkLBC(bc, bcShchepetkin):
                idxLocalZetaBC = self.getLocalBCIdx(self.zetaBC, idx)
                self.zetaBC.LBC[idxLocalZetaBC] |= bcAcquire


        # Save arrays of BCs indices to simplify computations later.
        idx = (((self.zetaBC.LBC & bcClosed) | (self.zetaBC.LBC & bcGradient)) != 0)
        self.zetaClosedOrGradientBCIdx1 = self.zetaBC.bcIdx1[idx]
        self.zetaClosedOrGradientBCIdx2 = self.zetaBC.bcIdx2[idx]

        idx = ((self.zetaBC.LBC & bcClamped) != 0)
        self.zetaClampedBCIdx1 = self.zetaBC.bcIdx1[idx]

        idx = ((self.ubarBC.LBC & bcClosed) != 0)
        self.ubarClosedBCIdx1 = self.ubarBC.bcIdx1[idx]

        idx = ((self.vbarBC.LBC & bcClosed) != 0)
        self.vbarClosedBCIdx1 = self.vbarBC.bcIdx1[idx]

        idx = ((self.ubarBC.LBC & bcClamped) != 0)
        self.ubarClampedBCIdx1 = self.ubarBC.bcIdx1[idx]

        idx = ((self.vbarBC.LBC & bcClamped) != 0)
        self.vbarClampedBCIdx1 = self.vbarBC.bcIdx1[idx]

        idx = ((self.ubarBC.LBC & bcGradient) != 0)
        self.ubarGradientBCIdx1 = self.ubarBC.bcIdx1[idx]
        self.ubarGradientBCIdx2 = self.ubarBC.bcIdx2[idx]

        idx = ((self.vbarBC.LBC & bcGradient) != 0)
        self.vbarGradientBCIdx1 = self.vbarBC.bcIdx1[idx]
        self.vbarGradientBCIdx2 = self.vbarBC.bcIdx2[idx]





        #  XXXXX DO THIS XXXXXXXXXX
        # self.VolCons     = decodeVolCon(west  = input.getVal('VolCons(west)'), east  = input.getVal('VolCons(east)'),
        #                                 north = input.getVal('VolCons(west)'), south = input.getVal('VolCons(west)'))


        self.nBCfiles    = input.getVal('NBCFILES', dtype = int, minVal = 0)
        self.max_Ffiles  = self.nBCfiles


        self.BRYids      = cp.zeros(self.max_Ffiles) - 1
        self.NBCcount    = cp.zeros(self.max_Ffiles)


        # XXXXX # REWRITE IN PYTHON
        # line = input.getVal('NBCFILES')
        # label='BRY - lateral open boundary conditions'
        # BRY =load_s2d(Nval, Cval, Cdim, line, label, ibcfile, igrid, nBCfiles, NBCcount, max_Ffiles)





        self.countBCNodes = 1000  # Total number of 2D boundary nodes.

        # Lateral boundary conditions apply switches.
        # These switches are part of the application grid and will be set to
        # FALSE elsewhere, if the boundary point is assigned by a nested grid.
        LBC_apply = cp.ones(self.countBCNodes, dtype = bool)

# WARNING:
#     We have to compute: zetaFSBCIdx, ubarBCIdx, vbarBCIdx, tBCIdx
#####################
#####################
#####################


        # IniVal = 0.0
        #  XXXXX DO THIS XXXXXXXXXX
        #
        # setBoundaryIdx = (iwest, ieast, inorth, isouth)
        #
        # # 2D boundary info.
        # # 2D here means that this information has no dependency on Z (or k index).
        #
        # BoundaryTypeIdx = cp.zeros(Nbounds, cp.Int32)
        # Boundary2DIndex = cp.zeros(Nbounds, cp.Int64)
        #
        # # !  Initialize module variables.
        #
        # # Acquire means "process lateral boundary data"
        # BOUNDARY.zetaFSBC[zetaFSBCIdx] = IniVal     # zetaFSBC: Free-surface (m)  boundary conditions.
        # BOUNDARY.ubarBC  [ubarBCIdx  ] = IniVal     # 2D u-momentum (m/s) boundary conditions.
        # BOUNDARY.vbarBC  [vbarBCIdx  ] = IniVal     # 2D v-momentum (m/s) boundary conditions
        # BOUNDARY.tBC     [tBCIdx     ] = IniVal     # TRACES boundary conditions. The T (and t) in isTvar, tBC and t_west means 'tracers'.







