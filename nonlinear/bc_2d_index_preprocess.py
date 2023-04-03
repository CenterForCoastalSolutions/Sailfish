import cupy as cp


# Rho Points
idxDstZGradBCR    = None
idxSrcZGradBCR    = None

idxDstPeriodicBCR = None
idxSrcPeriodicBCR = None


# U Points
idxDstZGradBCU    = None
idxSrcZGradBCR    = None

idxDstZGradBCU    = None
idxSrcZGradBCU    = None

idxDstClosedBCU   = None

idxDstPeriodicBCU = None
idxSrcPeriodicBCU = None


# V Points
idxDstZGradBCV    = None
idxSrcZGradBCV    = None

idxDstZGradBCV    = None
idxSrcZGradBCV    = None

idxDstClosedBCV   = None

idxDstPeriodicBCV = None
idxSrcPeriodicBCV = None




def precomputeBCIndices(ng):
    """ Precomputes indices (idx prefix) and weight coefficients (w prefix) for setting up the BC.
        The main idea is that the complicated stuff (computing indices and coefficients
        depending on the different BC choices and masks) is done at the beginning, leaving
        a very simple operation for the actual BC call.
    """

    global idxDstZGradBCR, idxSrcZGradBCR, idxDstPeriodicBCR, idxSrcPeriodicBCR, idxDstZGradBCU,
           idxSrcZGradBCR, idxDstZGradBCU, idxSrcZGradBCU, idxDstClosedBCU, idxDstPeriodicBCU,
           idxSrcPeriodicBCU, idxDstZGradBCV, idxSrcZGradBCV, idxDstZGradBCV, idxSrcZGradBCV,
           idxDstClosedBCV, idxDstPeriodicBCV, idxSrcPeriodicBCV

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

