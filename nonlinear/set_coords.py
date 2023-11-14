from misc import *
import cupy as cp

cosh = cp.cosh
sinh = cp.sinh
tanh = cp.tanh
log  = cp.log
exp  = cp.exp

def calcCs1(theta_b, theta_s, sc):
    # Original vertical strectching function, Song and Haidvogel (1994).
    # -----------------------------------------------------------------------
    # This vertical stretching function is defined as:
    #
    #     C(s) = (1 - b) * [sinh(s * a) / sinh(a)] + b * [-0.5 + 0.5 * tanh(a * (s + 0.5)) / tanh(0.5 * a)]
    #
    # where the stretching parameters (a, b) are specify at input:
    #
    #        a = theta_s               0 <  theta_s <= 8
    #        b = theta_b               0 <= theta_b <= 1
    #
    # If theta_b=0, the refinement is surface intensified as theta_s is increased.
    # If theta_b=1, the refinement in both bottom and surface is intensified as theta_s is increased.

    if theta_s != 0.0:
        cff1=1.0/cp.sinh(theta_s)
        cff2=0.5/cp.tanh(theta_s*0.5)

        w2 = 1.0 - theta_b
        w1 = theta_b
        return w2*cff1*sinh(theta_s*sc) + w1*(cff2*tanh(theta_s*(sc + 0.5)) - 0.5)

    else:
        return sc

def calcCs2(theta_b, theta_s, sc):
    # A. Shchepetkin vertical stretching function. This function was
    # improved further to allow bottom refinement (see Vstretching=4).
    #----------------------------------------------------------------------

    # This vertical stretching function is defined, in the simplest form,
    # as:
    #
    #     C(s) = [1.0 - cosh(theta_s * s)] / [cosh(theta_s) - 1.0]
    #
    # it is similar in meaning to the original vertical stretcing function
    # (Song and Haidvogel, 1994), but note that hyperbolic functions are
    # cosh, and not sinh.
    #
    # Note that the above definition results in
    #
    #        -1 <= C(s) <= 0
    #
    # as long as
    #
    #        -1 <= s <= 0
    #
    # and, unlike in any previous definition
    #
    #        d[C(s)]/ds  -->  0      if  s -->  0
    #
    # For the purpose of bottom boundary layer C(s) is further modified
    # to allow near-bottom refinement.  This is done by blending it with
    # another function.

    Aweight = 1.0
    Bweight = 1.0

    if sc == 0.0:
        return 0.0

    if theta_s > 0.0:
        Csur = (1.0 - cosh(theta_s*sc))/(cosh(theta_s) - 1.0)

        if theta_b > 0.0:
            Cbot = sinh(theta_b*(sc + 1.0))/sinh(theta_b) - 1.0
            Cweight = ((sc + 1.0)**Aweight)*(1.0 + (Aweight/Bweight)*(1.0 - (sc + 1.0)**Bweight))
            return Cweight*Csur + (1.0 - Cweight)*Cbot
        else:
            return Csur
    else:
        return sc

def calcCs3(theta_b, theta_s, sc):
    #   R. Geyer stretching function for high bottom boundary layer
    #   resolution.
    # -----------------------------------------------------------------------
    #
    #   This stretching function is intended for very shallow coastal
    #   applications, like gravity sediment flows.
    #
    #   At the surface, C(s=0)=0
    #
    #       C(s) = - LOG(cosh(Hscale * ABS(s) ** alpha)) /
    #                LOG(cosh(Hscale))
    #
    #   At the bottom, C(s=-1)=-1
    #
    #       C(s) = LOG(cosh(Hscale * (s + 1) ** beta)) /
    #              LOG(cosh(Hscale)) - 1
    #
    #   where
    #
    #        Hscale : scale value for all hyperbolic functions
    #                   Hscale = 3.0    set internally here
    #         alpha : surface stretching exponent
    #                   alpha = 0.65   minimal increase of surface resolution
    #                           1.0    significant amplification
    #          beta : bottoom stretching exponent
    #                   beta  = 0.58   no amplification
    #                           1.0    significant amplification
    #                           3.0    super-high bottom resolution
    #             s : stretched vertical coordinate, -1 <= s <= 0
    #                   s[k] = (k-N)/N       k=0:N,    W-points  (s_w)
    #                   s[k] = (k-N-0.5)/N   k=1:N,  RHO-points  (s_rho)
    #
    #   The stretching exponents (alpha, beta) are specify at input:
    #
    #          alpha = theta_s
    #          beta  = theta_b
    #

    if sc == 0.0:
        return 0.0

    exp_sur = theta_s
    exp_bot = theta_b
    Hscale = 3.0

    Cbot =  log(cosh(Hscale*(sc + 1.0)**exp_bot))/log(cosh(Hscale)) - 1.0
    Csur = -log(cosh(Hscale*abs(sc)**exp_sur))/log(cosh(Hscale))
    Cweight = 0.5*(1.0 - tanh(Hscale*(sc + 0.5)))

    return Cweight*Cbot + (1.0 - Cweight)*Csur

def calcCs4(theta_b, theta_s, sc):
    #   A. Shchepetkin improved double vertical stretching functions with
    #   bottom refiment.
    # -----------------------------------------------------------------------
    #
    #   The range of meaningful values for the control parameters are:
    #
    #        0 <  theta_s <= 10.0
    #        0 <= theta_b <=  3.0
    #
    #   Users need to pay attention to extreme r-factor (rx1) values near
    #   the bottom.
    #
    #   This vertical stretching function is defined, in the simplest form,
    #   as:
    #
    #       C(s) = [1.0 - cosh(theta_s * s)] / [cosh(theta_s) - 1.0]
    #
    #   it is similar in meaning to the original vertical stretcing function
    #   (Song and Haidvogel, 1994), but note that hyperbolic functions are
    #   cosh, and not sinh.
    #
    #   Note that the above definition results in
    #
    #          -1 <= C(s) <= 0
    #
    #   as long as
    #
    #          -1 <= s <= 0
    #
    #   and
    #
    #          d[C(s)]/ds  -->  0      if  s -->  0
    #
    #   For the purpose of bottom boundary layer C(s) is further modified
    #   to allow near-bottom refinement by using a continuous, second
    #   stretching function
    #
    #          C(s) = [EXP(theta_b * C(s)) - 1.0] / [1.0 - EXP(-theta_b)]
    #
    #   This double transformation is continuous with respect to "theta_s"
    #   and "theta_b", as both values approach to zero.

    if sc == 0.0:
        return 0.0

    if theta_s > 0.0:
        Csur = (1.0 - cosh(theta_s*sc))/(cosh(theta_s) - 1.0)
    else:
        Csur = -sc**2

    if theta_b > 0.0:
        return (exp(theta_b*Csur) - 1.0)/(1.0 - exp(-theta_b))

    else:
        return Csur

def calcCs5(theta_b, theta_s, sc):
    #  Stretching 5 case using a quadratic Legendre polynomial function
    #  aproach for the s-coordinate to enhance the surface exchange layer.
    #
    #  J. Souza, B.S. Powell, A.C. Castillo-Trujillo, and P. Flament, 2015:
    #    The Vorticity Balance of the Ocean Surface in Hawaii from a
    #    Regional Reanalysis.'' J. Phys. Oceanogr., 45, 424-440.
    #
    #  Added by Joao Marcos Souza - SOEST - 05/07/2012.
    # -----------------------------------------------------------------------

    if sc == 0.0:
        return 0.0

    if theta_s > 0.0:
        Csur = (1.0 - cosh(theta_s*sc))/(cosh(theta_s) - 1.0)
    else:
        Csur = -sc**2


    if theta_b > 0.0:
        return (exp(theta_b*Csur) - 1.0)/(1.0 - exp(-theta_b))
    else:
        return

def set_coords(GRID):
# This routine sets and initializes relevant variables associated
# with the vertical terrain-following coordinates transformation.
#
# Definitions:
#
#   N : Number of vertical levels for each nested grid.
#
#    zeta : time-varying free-surface, zeta(x,y,t), (m)
#
#       h : bathymetry, h(x,y), (m, positive, maybe time-varying)
#
#      hc : critical (thermocline, pycnocline) depth (m, positive)
#
#       z : vertical depths, z(x,y,s,t), meters, negative
#             z_w(x,y,0:N)      at   W-points  (top/bottom cell)
#             z_r(z,y,1:N)      at RHO-points  (cell center)
#
#             z_w(x,y,0    ) = -h(x,y)
#             z_w(x,y,N) = zeta(x,y,t)
#
#       s : nondimensional stretched vertical coordinate,
#            -1 <= s <= 0
#
#             s = 0   at the free-surface, z(x,y, 0,t) = zeta(x,y,t)
#             s = -1  at the bottom,       z(x,y,-1,t) = - h(x,y,t)
#
#             sc_w[k] = (k-N)/N       k=0:N,    W-points
#             sc_r[k] = (k-N-0.5)/N   k=1:N,  RHO-points
#
#       C : nondimensional vertical stretching function, C(s),
#             -1 <= C(s) <= 0
#
#             C(s) = 0    for s = 0,  at the free-surface
#             C(s) = -1   for s = -1, at the bottom
#
#             Cs_w[k] = F(s,theta_s,theta_b)  k=0:N,    W-points
#             Cs_r[k] = C(s,theta_s,theta_b)  k=1:N,  RHO-points
#
#      Zo : vertical transformation functional, Zo(x,y,s):
#
#             Zo(x,y,s) = H(x,y)C(s)      separable functions
#
#
# Two vertical transformations are supported, z => z(x,y,s,t):
#
# (1) Original transformation (Shchepetkin and McWilliams, 2005): In
#     ROMS since 1999 (version 1.8):
#
#       z(x,y,s,t) = Zo(x,y,s) + zeta(x,y,t) * [1 + Zo(x,y,s)/h(x,y)]
#
#     where
#
#       Zo(x,y,s) = hc * s + [h(x,y) - hc] * C(s)
#
#       Zo(x,y,s) = 0         for s = 0,  C(s) = 0,  at the surface
#       Zo(x,y,s) = -h(x,y)   for s = -1, C(s) = -1, at the bottom
#
# (2) New transformation: In UCLA-ROMS since 2005:
#
#       z(x,y,s,t) = zeta(x,y,t) + [zeta(x,y,t) + h(x,y)] * Zo(x,y,s)
#
#     where
#
#       Zo(x,y,s) = [hc * s[k] + h(x,y) * C[k]] / [hc + h(x,y)]
#
#       Zo(x,y,s) = 0         for s = 0,  C(s) = 0,  at the surface
#       Zo(x,y,s) = -1        for s = -1, C(s) = -1, at the bottom
#
#     At the rest state, corresponding to zero free-surface, this
#     transformation yields the following unperturbed depths, zhat:
#
#       zhat = z(x,y,s,0) = h(x,y) * Zo(x,y,s)
#
#            = h(x,y) * [hc * s[k] + h(x,y) * C[k]] / [hc + h(x,y)]
#
#     and
#
#       d(zhat) = ds * h(x,y) * hc / [hc + h(x,y)]
#
#     As a consequence, the uppermost grid box retains very little
#     dependency from bathymetry in the areas where hc << h(x,y),
#     that is deep areas. For example, if hc=250 m, and  h(x,y)
#     changes from 2000 to 6000 meters, the uppermost grid box
#     changes only by a factor of 1.08 (less than 10.).
#
#     Notice that:
#
#     * Regardless of the design of C(s), transformation (2) behaves
#       like equally-spaced sigma-coordinates in shallow areas, where
#       h(x,y) << hc.  This is advantageous because high vertical
#       resolution and associated CFL limitation is avoided in these
#       areas.
#
#     * Near-surface refinement is close to geopotential coordinates
#       in deep areas (level thickness do not depend or weakly-depend
#       on the bathymetry).  Contrarily,  near-bottom refinement is
#       like sigma-coordinates with thicknesses roughly proportional
#       to depth reducing high r-factors in these areas.
#
#
# This generic transformation design facilitates numerous vertical
# stretching functions, C(s).  These functions are set-up in this
# routine in terms of several stretching parameters specified in
# the standard input file.
#
# C(s) vertical stretching function properties:
#
# * a nonlinear, monotonic function
# * a continuous differentiable function, or
# * a piecewise function with smooth transition and differentiable
# * must be constrained by -1 <= C(s) <= 0, with C(0)=0 at the
#   free-surface and C(-1)=-1 at the bottom (bathymetry).
#
# References:
#
#   Shchepetkin, A.F. and J.C. McWilliams, 2005: The regional oceanic
#        modeling system (ROMS): a split-explicit, free-surface,
#        topography-following-coordinate oceanic model, Ocean
#        Modelling, 9, 347-404.
#
#   Song, Y. and D. Haidvogel, 1994: A semi-implicit ocean
#        circulation model using a generalized topography-
#        following coordinate system,  J.  Comp.  Physics,
#        115, 228-244.
#

    # -----------------------------------------------------------------------
    #   Set thickness controlling vertical coordinate stretching.
    # -----------------------------------------------------------------------
    #
    #   Set hc <= hmin, in the original formulation (Vtransform=1) to avoid
    #   [h(x,y)-hc] to be negative which results in dz/ds to be negative.
    #   Notice that this restriction is REMOVED in the new transformation
    #   (Vtransform=2): hc can be any value. It works for both hc < hmin
    #   and hc > hmin.

    if GRID.Vtransform == 1:
        hc = min(GRID.hmin, GRID.Tcline)

    elif GRID.Vtransform == 2:
        hc = GRID.Tcline

    N = GRID.N


    GRID.sc_w[0] = -1.0
    GRID.Cs_w[0] = -1.0

    ds = 1.0/N

    if GRID.Vstretching == 1:

        for k in range(1,N+1):
            GRID.sc_w[k] = ds*(k - N)
            GRID.sc_r[k] = ds*(k - N - 0.5)

            GRID.Cs_w[k] = calcCs1(GRID.theta_b, GRID.theta_s, GRID.sc_w[k])
            GRID.Cs_r[k] = calcCs1(GRID.theta_b, GRID.theta_s, GRID.sc_r[k])


    elif GRID.Vstretching == 2:

        for k in range(1,N+1):

            GRID.sc_w[k] = ds*(k - N)
            GRID.sc_r[k] = ds*(k - N - 0.5)

            GRID.Cs_w[k] = calcCs2(GRID.theta_b, GRID.theta_s, GRID.sc_w[k])
            GRID.Cs_r[k] = calcCs2(GRID.theta_b, GRID.theta_s, GRID.sc_r[k])


    elif GRID.Vstretching == 3:

        for k in range(1,N+1):

            GRID.sc_w[k] = ds*(k - N)
            GRID.sc_r[k] = ds*(k - N - 0.5)

            GRID.Cs_w[k] = calcCs3(GRID.theta_b, GRID.theta_s, GRID.sc_w[k])
            GRID.Cs_r[k] = calcCs3(GRID.theta_b, GRID.theta_s, GRID.sc_r[k])


    elif GRID.Vstretching == 4:

        for k in range(0, N):
            GRID.sc_w[k] = ds*(k - N + 1)
            GRID.sc_r[k] = ds*(k - N + 0.5)

            GRID.Cs_w[k] = calcCs4(GRID.theta_b, GRID.theta_s, GRID.sc_w[k])
            GRID.Cs_r[k] = calcCs4(GRID.theta_b, GRID.theta_s, GRID.sc_r[k])


    elif GRID.Vstretching == 5:

        for k in range(1, N+1):
            k2 = k + 0.5
            GRID.sc_w[k] = -(k*k   - 2.0*k*N  + k  + N*N - N)/(N*N - N) - 0.01*(k *k  - k *N)/(1.0 - N)
            GRID.sc_r[k] = -(k2*k2 - 2.0*k2*N + k2 + N*N - N)/(N*N - N) - 0.01*(k2*k2 - k2*N)/(1.0 - N)

            GRID.Cs_w[k] = calcCs5(GRID.theta_b, GRID.theta_s, GRID.sc_w[k])
            GRID.Cs_r[k] = calcCs5(GRID.theta_b, GRID.theta_s, GRID.sc_r[k])



    # Report information about vertical transformation.
    # -----------------------------------------------------------------------
    LwrtInfo = True
    if LwrtInfo:

        msgInfo('Vertical S-coordinate System:\nlevel   S-coord     Cs-curve   Z  at hmin       at hc    half way     at hmax')
        cff = 0.5*(GRID.hmax + GRID.hmin)

        for k in range(N-1,-1,-1):
            if GRID.Vtransform == 1:
                zhc = GRID.hc*GRID.sc_w[k]
                z1 = zhc + (GRID.hmin - GRID.hc)*GRID.Cs_w[k]
                z2 = zhc + (cff       - GRID.hc)*GRID.Cs_w[k]
                z3 = zhc + (GRID.hmax - GRID.hc)*GRID.Cs_w[k]

            elif GRID.Vtransform == 2:
                z1 = GRID.hmin*(GRID.hc*GRID.sc_w[k] + GRID.hmin*GRID.Cs_w[k])/(GRID.hc + GRID.hmin)
                z2 = cff      *(GRID.hc*GRID.sc_w[k] + cff      *GRID.Cs_w[k])/(GRID.hc + cff )
                z3 = GRID.hmax*(GRID.hc*GRID.sc_w[k] + GRID.hmax*GRID.Cs_w[k])/(GRID.hc + GRID.hmax)

                if hc > GRID.hmax:
                    zhc = z3      ## same as hmax, other values do not make sense
                else:
                    zhc = 0.5*hc*(GRID.sc_w[k] + GRID.Cs_w[k])



            msgInfo('%.3i     %12.6f     %12.6f     %12.6f     %12.6f     %12.6f     %12.6f' %
                    (k, GRID.sc_w[k], GRID.Cs_w[k], z1, zhc, z2, z3))


