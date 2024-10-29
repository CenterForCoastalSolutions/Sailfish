from misc import *
from math import pi, cos
import cupy as cp


# Power-law shape filter parameters for time-averaging of barotropic
# Fields.  The power-law shape filters are given by:
#
#    F(xi) = xi^Falpha*(1 - xi^Fbeta) - Fgamma*xi
#
# Possible settings of parameters to yield the second-order accuracy:
#
#    Falpha  Fbeta      Fgamma
#    ------------------------------
#     2.0     1.0    0.1181  0.169     The problem here is setting
#     2.0     2.0    0.1576  0.234     Fgamma. Its value here is
#     2.0     3.0    0.1772  0.266     understood as the MAXIMUM
#     2.0     4.0    0.1892  0.284     allowed. It is computed using
#     2.0     5.0    0.1976  0.296     a Newton iteration scheme.
#     2.0     6.0    0.2039  0.304
#     2.0     8.0    0.2129  0.314
#
# NOTE: Theoretical values of Fgamma presented in the table above are
# derived assuming "exact" barotropic mode stepping. Consequently, it
# does not account for effects caused by Forward-Euler (FE) startup
# of the barotropic mode at every 3D time step.  As the result, the
# code may become unstable if the theoretical value of Fgamma is used
# when mode splitting ratio "ndtfast" is small, thus yielding non-
# negligible start up effects.  To compensate this, the accepted
# value of Fgamma is reduced relatively to theoretical one, depending
# on splitting ratio "ndtfast".  This measure is empirical. It is
# shown to work with setting of "ndtfast" as low as 15, which is
# more robust that the Hamming Window the squared cosine weights
# options in "set_weights".
#

Falpha = 2.0
Fbeta  = 4.0
Fgamma = 0.284

POWER_LAW = True
COSINE2 = False


def set_weights(compTimes):
# This routine sets the weight functions for the time averaging of
# 2D fields over all short time-steps.


    # Compute time-averaging filter for barotropic fields.
    # =======================================================================

    msgInfo('Computing time averaging weights')

    # Initialize both sets of weights to zero.
    compTimes.nfast = 0
    compTimes.weight1[:] = 0.0
    compTimes.weight2[:] = 0.0

    ndtfast = compTimes.ndtfast


    if POWER_LAW:

        # Power-law shape filters.
        # -----------------------------------------------------------------------
        #
        #   The power-law shape filters are given by:
        #
        #      F(xi) = xi^Falpha*(1 - xi^Fbeta) - Fgamma*xi
        #
        #   where xi=scale*i/ndtfast; and scale, Falpha, Fbeta, Fgamma, and
        #   normalization are chosen to yield the correct zeroth-order
        #   (normalization), first-order (consistency), and second-order moments,
        #   resulting in overall second-order temporal accuracy for time-averaged
        #   barotropic motions resolved by baroclinic time step. There parameters
        #   are set in "mod_scalars".
        scale = (Falpha + 1.0)*(Falpha + Fbeta + 1.0)/((Falpha + 2.0)*(Falpha + Fbeta + 2.0)*ndtfast)

        # Find center of gravity of the primary weighting shape function and
        # iteratively adjust "scale" to place the  centroid exactly at
        # "ndtfast".

        gamma = Fgamma*max(0.0, 1.0 - 10.0/ndtfast)
        for iter in range(16):
            compTimes.nfast = 0
            for i in range(2*ndtfast):
                cff = scale*(i+1)
                compTimes.weight1[i] = (cff**Falpha) - (cff**(Falpha+Fbeta)) - gamma*cff
                if compTimes.weight1[i] > 0.0:
                    compTimes.nfast = i + 1

                if (compTimes.nfast > 0) and (compTimes.weight1[i] < 0.0):
                    compTimes.weight1[i] = 0.0

            wsum  = 0.0
            shift = 0.0
            for i in range(compTimes.nfast):
                wsum  +=       compTimes.weight1[i]
                shift += (i+1)*compTimes.weight1[i]

            scale = scale*shift/(wsum*ndtfast)


    elif COSINE2:

        # Cosine-squared shaped filter.
        # -----------------------------------------------------------------------

        cff = pi/ndtfast
        for i in range(2*ndtfast):
            arg = cff*(i + 1 - ndtfast)
            if (2.0*abs(arg) < pi):
                compTimes.weight1[i] = 0.0882 + cos(arg)**2    # hamming window
                compTimes.nfast = i + 1


    # Post-processing of primary weights.
    # -----------------------------------------------------------------------

    # Although it is assumed that the initial settings of the primary
    # weights has its center of gravity "reasonably close" to NDTFAST,
    # it may be not so according to the discrete rules of integration.
    # The following procedure is designed to put the center of gravity
    # exactly to NDTFAST by computing mismatch (NDTFAST-shift) and
    # applying basically an upstream advection of weights to eliminate
    # the mismatch iteratively. Once this procedure is complete primary
    # weights are normalized.

    # Find center of gravity of the primary weights and subsequently
    # calculate the mismatch to be compensated.

    for iter in range(ndtfast):
        wsum  = 0.0
        shift = 0.0
        for i in range(compTimes.nfast):
            wsum  +=       compTimes.weight1[i]
            shift += (i+1)*compTimes.weight1[i]

        shift = shift/wsum
        cff = ndtfast - shift

        # Apply advection step using either whole, or fractional shifts.
        # Notice that none of the four loops here is reversible.

        if cff > 1.0:
            compTimes.nfast += 1
            for i in range(compTimes.nfast-1,1,-1):
                compTimes.weight1[i] = compTimes.weight1[i-1]

            compTimes.weight1[0] = 0.0

        elif cff > 0.0:
            wsum = 1.0 - cff

            for i in range(compTimes.nfast-1,1,-1):
                compTimes.weight1[i] = wsum*compTimes.weight1[i] + cff*compTimes.weight1[i-1]

            compTimes.weight1[0] = wsum*compTimes.weight1[0]

        elif cff < -1.0:
            compTimes.nfast -= 1
            for i in range(compTimes):
                compTimes.weight1[i] = compTimes.weight1[i+1]

            compTimes.weight1[compTimes.nfast] = 0.0

        elif cff < 0.0:
            wsum = 1.0 + cff
            for i in range(compTimes.nfast-1):
                compTimes.weight1[i] = wsum*compTimes.weight1[i] - cff*compTimes.weight1[i+1]

            compTimes.weight1[compTimes.nfast] = wsum*compTimes.weight1[compTimes.nfast]


    # Set SECONDARY weights assuming that backward Euler time step is used
    # for free surface.  Notice that array weight(2,i,ng) is assumed to
    # have all-zero status at entry in this segment of code.

    for j in range(compTimes.nfast):
        cff = compTimes.weight1[j]
        for i in range(j):
          compTimes.weight2[i] += cff

    # Normalize both set of weights.
    wsum = 0.0
    cff = 0.0
    for i in range(compTimes.nfast):
        wsum += compTimes.weight1[i]
        cff  += compTimes.weight2[i]

    wsum = 1.0/wsum
    cff  = 1.0/cff
    for i in range(compTimes.nfast):
        compTimes.weight1[i] = wsum*compTimes.weight1[i]
        compTimes.weight2[i] =  cff*compTimes.weight2[i]


    # Report weights.
    if LwrtInfo:
        msgInfo('Time Splitting Weights. ndtfast = %i, nfast = %i' % (ndtfast, compTimes.nfast))
        msgInfo('==================================')
        msgInfo('Primary      Secondary     Accumulated to Current Step')

        cff  = 0.0
        cff1 = 0.0
        cff2 = 0.0
        wsum = 0.0
        shift = 0.0
        for i in range(compTimes.nfast):
            cff  +=             compTimes.weight1[i]
            cff1 +=       (i+1)*compTimes.weight1[i]
            cff2 += (i+1)*(i+1)*compTimes.weight1[i]

            wsum  +=         compTimes.weight2[i]
            shift += (i+0.5)*compTimes.weight2[i]
            msgInfo('%i  %f %f  %f  %f'  % (i, compTimes.weight1[i], compTimes.weight2[i], cff, wsum))

        cff1 = cff1/ndtfast
        cff2 = cff2/(ndtfast*ndtfast)
        shift = shift/ndtfast

        msgInfo('ndtfast, nfast = %i, %i    nfast/ndtfast = %i' % (ndtfast, compTimes.nfast, compTimes.nfast/ndtfast))

        if POWER_LAW:
            msgInfo ('Centers of gravity and integrals (values must be 1, 1, approx 1/2, 1, 1): %f, %f, %f, %f, %f' %
                     (cff1, cff2, shift, cff, wsum))
            msgInfo ('Power filter parameters, Fgamma, gamma = %f, %f' % (Fgamma, gamma))

        if (cff2 < 1.0001):
            msgInfo('WARNING: unstable weights, reduce parameter Fgamma in mod_scalars.F')
