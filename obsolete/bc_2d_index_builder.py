
def buildBCIndices(Istr, Iend, Jstr, Jend, mask):
    # Zero gradient BC is implemented by setting one node value equal to another.
    # In other words: dA/dx = 0  ->  (A[i+1] - A[i])/Dx = 0  ->  A[i+1] = A[i]
    # At the corners, the derivative is taken in the direction of the grid cell diagonal.

    domain = DOMAIN
    lbc_apply = LBC_apply

    zeroGradBCIdx    = []
    zeroGradBCIdxSrc = []

    # East-West gradient boundary conditions.
    if not EWperiodic[ng]:
        zeroGradBCIdx    += flattenedIdx(lbc_apply.east, Iend + 1, mask)
        zeroGradBCIdxSrc += flattenedIdx(lbc_apply.east, Iend    , mask)

        zeroGradBCIdx    += flattenedIdx(lbc_apply.west, Istr - 1, mask)
        zeroGradBCIdxSrc += flattenedIdx(lbc_apply.west, Istr    , mask)


    # North-South gradient boundary conditions.
    if not EWperiodic[ng]:
        zeroGradBCIdx    += flattenedIdx(Jend + 1, lbc_apply.north, mask)
        zeroGradBCIdxSrc += flattenedIdx(Jend    , lbc_apply.north, mask)

        zeroGradBCIdx    += flattenedIdx(Jstr - 1, lbc_apply.south, mask)
        zeroGradBCIdxSrc += flattenedIdx(Jend    , lbc_apply.south, mask)


    # Boundary corners.
    # In the cases where gradient bc is applied to both sides of a corner,
    # the value in the corner is taken as the closest node in the
    # diagonal direction.
    if not (EWperiodic[ng] or NSperiodic[ng]):
        if lbc_apply.south[Istr - 1] and lbc_apply.west[Jstr - 1]:
            zeroGradBCIdxSrc += flattenedIdx(Jstr,     Istr,     mask)
            zeroGradBCIdx    += flattenedIdx(Jstr - 1, Istr - 1, mask)

        if lbc_apply.south[Istr - 1] and lbc_apply.east[Jend + 1]:
            zeroGradBCIdxSrc += flattenedIdx(Jend,     Istr,     mask)
            zeroGradBCIdx    += flattenedIdx(Jend + 1, Istr - 1, mask)

        if lbc_apply.north[Istr + 1] and lbc_apply.west[Jstr - 1]:
            zeroGradBCIdxSrc += flattenedIdx(Jstr,     Iend,     mask)
            zeroGradBCIdx    += flattenedIdx(Jstr - 1, Iend + 1, mask)

        if lbc_apply.north[Istr + 1] and lbc_apply.east[Jend + 1]:
            zeroGradBCIdxSrc += flattenedIdx(Jend,     Iend,     mask)
            zeroGradBCIdx    += flattenedIdx(Jend + 1, Iend + 1, mask)



    # Build indices for periodic boundary conditions.
    # -----------------------------------------------

    L = Lm + 1
    M = Mm + 1

    periodicBCIdx = []
    periodicBCIdxSrc = []

    iBCRim =    [*range(-2, 1) ] + \
                [*range(M + 1, M + NghostPoints + 1)]
    iBCRimSrc = [*range(M - 2, M + 1)] + \
                [*range(1, NghostPoints + 1)]

    jBCRim =    [*range(-2, 1) ] + \
                [*range(L + 1, L + NghostPoints + 1)]
    jBCRimSrc = [*range(L - 2, L + 1)] + \
                [*range(1, NghostPoints + 1)]


    # East-West periodic boundary conditions.
    if EWperiodic[ng]:
        if NSperiodic[ng]:
            jRange = [*range(Jstr, jEnd)]
        else:
            jRange = [*range(Jstr, jEndR)]


        for i, iSrc in zip(iBCRim, iBCRimSrc):
            periodicBCIdx    += flattenedIdx(jRange, i)
            periodicBCIdxSrc += flattenedIdx(jRange, iSrc)


    # North-South periodic boundary conditions.
    if NSperiodic[ng]:
        if EWperiodic[ng]:
            if NSperiodic[ng]:
                iRange = [*range(Istr, Iend)]
            else:
                iRange = [*range(Istr, IendR)]

        for j, jSrc in zip(jBCRim, jBCRimSrc):
            periodicBCIdx    += flattenedIdx(j,    iRange)
            periodicBCIdxSrc += flattenedIdx(jSrc, iRange)



    # Boundary corners.
    if EWperiodic[ng] and NSperiodic[ng]:
        if EW_exchange[ng] and NS_exchange[ng]:

            for j, jSrc in zip(jBCRim, jBCRimSrc):
                for i, iSrc in zip(iBCRim, iBCRimSrc):
                    periodicBCIdx    += flattenedIdx(j,    i)
                    periodicBCIdxSrc += flattenedIdx(jSrc, iSrc)



