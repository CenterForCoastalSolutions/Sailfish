# !================================================== Hernan G. Arango ===
# !  Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
# !    Licensed under a MIT/X style license                              !
# !    See License_ROMS.txt                                              !
# !=======================================================================
# !                                                                      !
# !  These routines apply close, gradient or periodic boundary           !
# !  conditions to generic 2D fields.                                    !
# !                                                                      !
# !  On Input:                                                           !
# !                                                                      !
# !     ng                Nested grid number.                            !
# !     tile              Domain partition.                              !
# !     LBi               I-dimension Lower bound.                       !
# !     UBi               I-dimension Upper bound.                       !
# !     LBj               J-dimension Lower bound.                       !
# !     UBj               J-dimension Upper bound.                       !
# !     A                 2D field.                                      !
# !                                                                      !
# !  On Output:                                                          !
# !                                                                      !
# !     A                 Processed 2D field.                            !
# !                                                                      !
# !  Routines:                                                           !
# !                                                                      !
# !     bc_r:-       Boundary conditions for field at RHO-points    !
# !     bc_u2d_tile       Boundary conditions for field at U-points      !
# !     bc_v2d_tile       Boundary conditions for field at V-points      !
# !                                                                      !
# !=======================================================================

import cupy as cp

def precomputeBC(ng):
    """ Precomputes indices (idx prefix) and weight coefficients (w prefix) for setting up the BC.
        The main idea is that the complicated stuff (computing indices and coefficients
        depending on the different BC choices and masks) is done at the beginning, leaving
        a very simple operation for the actual BC call.
    """

    shpR = R.shape
    # I use a method that is not very efficient but it is more clear. It consists in
    # creating helper arrays (h prefix) that we initially fill with zeros and then with the indices of the
    # nodes they depend on.
    # These helper arrays are only used during this preprocessing stage and then disposed.
    indices = cp.indices(shp)
    flatIdx = cp.ravel_multi_index(indices. shp)  # This is an array to convert between flat and 3-D arrays.

    hRZGradBC = cp.zeros(shp, np.int32) - 1


    # Mask of BC nodes at each mesh side (E, W, N, S). If there are periodic boundary conditions in a pair of sides
    # These take precedence, meaning that these mask will be all false in these sides.
    maskLBCe = LBC_apply[ng].east [Jstr:Jend] & ~EWperiodic[ng]
    maskLBCw = LBC_apply[ng].west [Jstr:Jend] & ~EWperiodic[ng]
    maskLBCn = LBC_apply[ng].north[Istr:Iend] & ~NSperiodic[ng]
    maskLBCs = LBC_apply[ng].south[Istr:Iend] & ~NSperiodic[ng]

    hRZGradBC[Iend + 1, maskLBCe] = flatIdx[Iend, maskLBCe]
    hRZGradBC[Istr - 1, maskLBCw] = flatIdx[Istr, maskLBCw]
    hRZGradBC[maskLBCn, Jend + 1] = flatIdx[maskLBCn, Jend]
    hRZGradBC[maskLBCs, Jstr - 1] = flatIdx[maskLBCs, Jstr]


    # BC at the corners. In the original code they did an average of two cell that should
    # have identical values, so we are just taking the one value.
    if hRZGradBC[Istr - 1, Jstr - 1] >=0:
        hRZGradBC[Istr - 1, Jstr - 1] = flatIdx[Istr, Jstr]

    if hRZGradBC[Iend + 1, Jstr - 1] >= 0:
        hRZGradBC[Iend - 1, Jstr - 1] = flatIdx[Iend, Jstr]

    if hRZGradBC[Istr - 1, Jend + 1] >=0:
        hRZGradBC[Istr - 1, Jend + 1] = flatIdx[Istr, Jend]

    if hRZGradBC[Iend - 1, Jend + 1] >=0:
        hRZGradBC[Iend - 1, Jend + 1] = flatIdx[Iend, Jend]


    # U-Points


    # This masks for Zero flow conditions are all True or all False.
    maskZe    = np.full(Jend-Jstr,  LBC[ieast,  isBu2d, ng].closed and not EWperiodic[ng]
    maskZw    = np.full(Jend-Jstr,  LBC[iwest,  isBu2d, ng].closed and not EWperiodic[ng]
    maskZn    = np.full(Iend-Istr,  LBC[inorth, isBu2d, ng].closed and not NWperiodic[ng]
    maskZs    = np.full(Iend-Istr,  LBC[isouth, isBu2d, ng].closed and not NWperiodic[ng]

    maskGrade = np.full(Jend-Jstr, ~LBC[ieast,  isBu2d, ng].closed and not EWperiodic[ng]
    maskGradw = np.full(Jend-Jstr, ~LBC[iwest,  isBu2d, ng].closed and not EWperiodic[ng]
    maskGradn = np.full(Iend-Istr, ~LBC[inorth, isBu2d, ng].closed and not NWperiodic[ng]
    maskGrads = np.full(Iend-Istr, ~LBC[isouth, isBu2d, ng].closed and not NWperiodic[ng]

    hUClosedBC[Iend+1, maskZe & maskLBCe ] = flatIdx[Iend,   maskZe & maskLBCe]
    hUClosedBC[Istr,   maskZw & maskLBCw ] = flatIdx[Istr-1, maskZw & maskLBCw]
    hUClosedBC[maskZn & maskLBCn, jend+1 ] = flatIdx[maskZn & maskLBCn, Jend  ]
    hUClosedBC[maskZs & maskLBCn, jstr   ] = flatIdx[maskZs & maskLBCs, Jstr-1]

    hUGradBC[Iend+1, maskZe & maskGrade  ] = flatIdx[Iend,   maskGrade & maskLBCe]
    hUGradBC[Istr,   maskZw & maskGradw  ] = flatIdx[Istr-1, maskGradw & maskLBCw]
    hUGradBC[maskGradn & maskLBCn, jend+1] = flatIdx[maskGradn & maskLBCn, Jend  ]
    hUGradBC[maskGrads & maskLBCn, jstr  ] = flatIdx[maskGrads & maskLBCs, Jstr-1]
    if (EWperiodic(ng)):
        Imin = IstrU
        Imax = Iend
        ELSE
        Imin = Istr
        Imax = IendR

    SORT Indices




def bc_r2d_tile(ng, tile, LBi, UBi, LBj, UBj, vars):
    """ BC for rho-type cells: Boundary conditions are "imposed" by setting values to ghost  nodes"""

    for var in vars:
        # Sets zero-gradient boundary conditions.
        var.flat[idxDstZGradBCR] = A.flat[idSrcZGradBCR]

        # Sets periodic boundary conditions.
        var.flat[idxDstPeriodicBCR] = A.flat[idxSrcPeriodicBCR]



def bc_u2d_tile(ng, tile, LBi, UBi, LBj, UBj, vars):
    """ BC for U-type cells:  Boundary conditions are "imposed" by setting values to ghost  nodes"""

    for var in vars:
        # Sets zero-gradient boundary conditions.
        var.flat[idxDstZGradBCU] = A.flat[idSrcZGradBCR]

        # Sets free-slip boundary conditions.
        var.flat[idxDstZGradBCU] = gamma2*A.flat[idSrcZGradBCU]

        # Sets no-slip boundary conditions.
        var.flat[idxDstClosedBCU] = 0.0

        # Sets periodic boundary conditions.
        var.flat[idxDstPeriodicBCU] = A.flat[idxSrcPeriodicBCU]


def bc_u2d_tile(ng, tile, LBi, UBi, LBj, UBj, vars):
    """ BC for U-type cells:  Boundary conditions are "imposed" by setting values to ghost  nodes"""

    for var in vars:
        # Sets zero-gradient boundary conditions.
        var.flat[idxDstZGradBCV] = A.flat[idSrcZGradBCV]

        # Sets free-slip boundary conditions.
        var.flat[idxDstZGradBCV] = gamma2*A.flat[idSrcZGradBCV]

        # Sets no-slip boundary conditions.
        var.flat[idxDstClosedBCV] = 0.0

        # Sets periodic boundary conditions.
        var.flat[idxDstPeriodicBCV] = A.flat[idxSrcPeriodicBCV]

