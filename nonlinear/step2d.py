import mod_comptimes
from mod_operators import *
from barotropicVelocityBC import barotropicVelocityBC
from zetabc               import zetabc
from mod_constants import *
from misc          import *


def computeZetaRHS(zeta, h, ubar, vbar, GRID):

    #Apply mass point sources (volume vertical influx), if any
    # TODO: Implement XXXX
    # if LwSrc:
    #     ij = SOURCES.IJsrc   # This has to be precomputed as i + j*width
    #     zeta[knew,:,:].ravel()[ij] += SOURCES.Qbar*pm.ravel()[ij]*pn.ravel()[ij]*dtfast
    # addForcings()

    # compute the water column depth
    D = zeta + h


    DU = ubar*RtoU(D, ubar)   # TODO: Remember to check if we can remove the extra parameter (ubar)
    DV = vbar*RtoV(D, vbar)
    # DU[:] = DU + 1
    # DV[:] = DV + 2
    # print(2, DU.shape, DU)
    # a = divUVtoR(DU, DV, D, GRID)
    # print (3,a)
    return divUVtoR(DU, DV, D, GRID)   # TODO: Remember to check if we can remove the extra parameter (D)


def computeMomentumRHS(h, gzeta, U, V):
    rhs_ubar = -0.5*g*RtoU(h, U)*DξRtoU(gzeta + gzeta*gzeta, U)
    rhs_vbar = -0.5*g*RtoV(h, V)*DηRtoV(gzeta + gzeta*gzeta, V)

    # if UV_ADV:
    #     #!---------------------------------------------------------------------------
    #     #! Contribution of a term corresponding to product of
    #     #! Stokes and Eulerian Velocity Eqn. 26 and 27.
    #     #! This removes terms that were unneccessarily added in flux form.
    #     #!---------------------------------------------------------------------------
    #     cff         = 0.5 * (Drhs(i-1,j) + Drhs(i,j))
    #     DUSon(i,j)  = cff * on_u(i, j) * ubar_stokes(i, j)
    #     DVSon(i, j) = 0.25 * cff * on_u(i, j) * (vbar_stokes(i  ,j  ) + vbar_stokes(i  ,j+1) + vbar_stokes(i-1,j  ) + vbar_stokes(i-1,j+1))
    #
    #     UFx(i,j)    = 0.5 * (ubar(i, j, krhs) + ubar(i+1, j, krhs))
    #     VFx(i,j)    = 0.5 * (vbar(i, j, krhs) + vbar(i, j+1, krhs))
    #
    #     cff         = 0.5 * (Drhs(i,j) + Drhs(i,j-1))
    #     DUSom(i,j)  = cff * 0.25 * om_v(i,j) * (ubar_stokes(i,j) + ubar_stokes(i+1, j) + ubar_stokes(i  , j-1) + ubar_stokes(i+1, j-1))
    #
    #     DVSom(i,j)    = cff * om_v(i, j) * vbar_stokes(i, j)
    #     cff           = 0.5 * (Drhs(i, j) + Drhs(i, j-1))
    #     UFe(i,j)      = 0.5 * (ubar(i+1, j, krhs) + ubar(i, j, krhs))
    #     VFe(i,j)      = 0.5 * (vbar(i, j, krhs) + vbar(i, j+1, krhs))
    #     cff1          = UFx(i,j) - UFx(i-1,j) # line 2215
    #     cff2          = VFx(i,j) - VFx(i-1,j)
    #     cff3          = DUSon(i,j) * cff1
    #     cff4          = DVSon(i,j) * cff2
    #     rhs_ubar(i,j) = rhs_ubar(i,j) + cff3 + cff4
    #
    #     cff1          = UFe(i,j) - UFe(i,j-1)
    #     cff2          = VFe(i,j) - VFe(i,j-1)
    #     cff3          = DUSom(i,j) * cff1
    #     cff4          = DVSom(i,j) * cff2
    #     rhs_vbar(i,j) = rhs_vbar(i,j) + cff3 + cff4

    # # add in bottom stress
    # rhs_ubar[:,:] -= bustr*om_u*on_u
    # rhs_vbar[:,:] -= bvstr*om_v*on_v

    #  Time step 2D momentum equations.

    # compute the water column depth at times "new" and "stp"

    return rhs_ubar, rhs_vbar



def step2dPredictor(compTimes, GRID, OCEAN, BOUNDARY):

    # Aliases
    zeta_t2, zeta_t1, zeta_t0 = (OCEAN.zeta_t2, OCEAN.zeta_t1, OCEAN.zeta_t0)
    ubar_t2, ubar_t1, ubar_t0 = (OCEAN.ubar_t2, OCEAN.ubar_t1, OCEAN.ubar_t0)
    vbar_t2, vbar_t1, vbar_t0 = (OCEAN.vbar_t2, OCEAN.vbar_t1, OCEAN.vbar_t0)
    rzeta_t1, rzeta_t0 = (OCEAN.rzeta_t1, OCEAN.rzeta_t0)
    rubar_t1, rubar_t0 = (OCEAN.rubar_t1, OCEAN.rubar_t0)
    rvbar_t1, rvbar_t0 = (OCEAN.rvbar_t1, OCEAN.rvbar_t0)
    h = GRID.h.ravel()

    Δt = compTimes.get2DTimes()


    # Free-surface equation.
    # =================================

    # During the first time-step, the predictor step is Forward-Euler. Otherwise, the predictor step is Leap-frog.

    rhs_zeta_t1 = computeZetaRHS(zeta_t1, h, ubar_t1, vbar_t1, GRID)


    if compTimes.isFirst2DStep():
        # The first time it performs a simple Euler time step. RHS is computed at tn and time derivatives are
        # (f(tn) - f(tn-1))/dtfast

        zeta_t2[:] = zeta_t1 + Δt*rhs_zeta_t1

        gzeta = 0.5*(zeta_t2 + zeta_t1)


    else:
        # If is not the first time step, the predictor consists of a leapfrog time step where the RHS is computed at tn
        # and time derivatives are centered at tn (f(tn+1) - f(tn-1))/2*dtfast.

        zeta_t2[:] = zeta_t1 + Δt*rhs_zeta_t1

        weight = 4.0/25.0
        gzeta = (1 - weight)*zeta_t1 + weight*0.5*(zeta_t2 + zeta_t0)

        # In the predictor step, save the rhs for future use.
        rzeta_t1[:] = rhs_zeta_t1



    # Apply free-surface lateral BC
    zetabc(zeta_t2, compTimes, BOUNDARY)




    # Momentum equations.
    # ===================


    #compute right-hand-side for the 2D momentum equations
    rhs_ubar, rhs_vbar = computeMomentumRHS(h, gzeta, ubar_t2, vbar_t2)


    # Interpolate depth at points U, V
    D_t2 = zeta_t2 + h
    D_t1U = RtoU(D_t1)
    D_t1V = RtoV(D_t1)
    D_t2U = RtoU(D_t2)
    D_t2V = RtoV(D_t2)


    # During the first time-step, the predictor step is Forward-Euler
    # and the corrector step is Backward-Euler. Otherwise, the predictor
    # step is Leap-frog and the corrector step is Adams-Moulton.


    ubar_t2[:] = (ubar_t1*(D_t1U + Δt*rhs_ubar))/D_t2U
    vbar_t2[:] = (vbar_t1*(D_t1V + Δt*rhs_vbar))/D_t2V

    # In the predictor step, save the rhs for future use.
    rubar_t1[:] = rhs_ubar
    rvbar_t1[:] = rhs_vbar


    # Apply lateral boundary conditions.
    barotropicVelocityBC(OCEAN, BOUNDARY)


def step2dCorrector(compTimes, GRID, OCEAN, BOUNDARY):

    # Aliases
    zeta_t2, zeta_t1, zeta_t0 = (OCEAN.zeta_t2, OCEAN.zeta_t1, OCEAN.zeta_t0)
    ubar_t2, ubar_t1, ubar_t0 = (OCEAN.ubar_t2, OCEAN.ubar_t1, OCEAN.ubar_t0)
    vbar_t2, vbar_t1, vbar_t0 = (OCEAN.vbar_t2, OCEAN.vbar_t1, OCEAN.vbar_t0)
    rzeta_t1, rzeta_t0 = (OCEAN.rzeta_t1, OCEAN.rzeta_t0)
    rubar_t1, rubar_t0 = (OCEAN.rubar_t1, OCEAN.rubar_t0)
    rvbar_t1, rvbar_t0 = (OCEAN.rvbar_t1, OCEAN.rvbar_t0)
    h = GRID.h

    Δt = compTimes.get2DTimes()

    # Free-surface equation.
    # =================================

    # During the first time-step,the corrector step is Backward-Euler. Otherwise, the corrector step is Adams-Moulton.

    rhs_zeta_t2 = computeZetaRHS(zeta_t2, h, ubar_t2, vbar_t2, GRID)

    # Adams-Moulton order 3
    zeta_t2[:] = zeta_t1 + Δt*(AM3_2*rhs_zeta_t2 + AM3_2*rzeta_t1 + AM3_2*rzeta_t0)

    weight = 2.0/5.0
    gzeta  = (1 - weight)*zeta_t2 + weight*zeta_t1




    # Apply free-surface lateral BC
    zetabc(zeta_t2, compTimes, BOUNDARY)


    #compute right-hand-side for the 2D momentum equations
    rhs_ubar, rhs_vbar = computeMomentumRHS(h, gzeta)


    # And interpolate them at points U, V
    D_t1 = zeta_t1 + h
    D_t2 = zeta_t2 + h
    D_t1U = RtoU(D_t1)
    D_t1V = RtoV(D_t1)
    D_t2U = RtoU(D_t2)
    D_t2V = RtoV(D_t2)


    # During the first time-step, the corrector step is Backward-Euler. Otherwise, the corrector step is Adams-Moulton.
    ubar_t2[:] = (ubar_t1*(D_t1U + Δt*(AM3_2*rhs_ubar + AM3_1*rubar_t1 + AM3_0*rubar_t0)))/D_t2U
    vbar_t2[:] = (vbar_t1*(D_t1V + Δt*(AM3_2*rhs_vbar + AM3_1*rvbar_t1 + AM3_0*rvbar_t0)))/D_t2V


    # Apply lateral boundary conditions.
    barotropicVelocityBC(OCEAN, BOUNDARY)


    # Compute integral mass flux across open boundaries and adjust for volume conservation.

    #       IF (ANY(VolCons(:,ng))) THEN
    #         CALL obc_flux (LBi, UBi, LBj, UBj,                         &
    #      &                      IminS, ImaxS, JminS, JmaxS,                 &
    #      &                      knew,                                       &
    # # ifdef MASKING
    #      &                      umask, vmask,                               &
    # # endif
    #      &                      h, om_v, on_u,                              &
    #      &                      ubar, vbar, zeta)
    #       END IF
    # !
    # !-----------------------------------------------------------------------
    # !  Apply momentum transport point sources (like river runoff), if any.
    # !-----------------------------------------------------------------------
    # !
    #       IF (LuvSrc(ng)) THEN
    #         DO is=1,Nsrc(ng)
    #           i=SOURCES(ng)%Isrc(is)
    #           j=SOURCES(ng)%Jsrc(is)
    #           IF (((IstrR.le.i).and.(i.le.IendR)).and.                      &
    #      &        ((JstrR.le.j).and.(j.le.JendR))) THEN
    #             IF (INT(SOURCES(ng)%Dsrc(is)).eq.0) THEN
    #               cff=1.0_r8/(on_u(i,j)*                                    &
    #      &                    0.5_r8*(zeta(i-1,j,knew)+h(i-1,j)+            &
    #      &                            zeta(i  ,j,knew)+h(i  ,j)))
    #               ubar(i,j,knew)=SOURCES(ng)%Qbar(is)*cff
    #             ELSE
    #               cff=1.0_r8/(om_v(i,j)*                                    &
    #      &                    0.5_r8*(zeta(i,j-1,knew)+h(i,j-1)+            &
    #      &                            zeta(i,j  ,knew)+h(i,j  )))
    #               vbar(i,j,knew)=SOURCES(ng)%Qbar(is)*cff
    #             END IF
    #           END IF
    #         END DO
    #       END IF
    # !


