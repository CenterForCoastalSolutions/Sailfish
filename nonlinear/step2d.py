import mod_comptimes
from mod_operators import *
from barotropicVelocityBC import barotropicVelocityBC
from zetabc               import zetabc
from mod_constants import *
from misc          import *
import time




# import rmm
# pool = rmm.mr.PoolMemoryResource(
#     rmm.mr.ManagedMemoryResource(),
#     initial_pool_size=55*(2**30),
#     maximum_pool_size=55*(2**30)
# )
# rmm.mr.set_current_device_resource(pool)
# cp.cuda.set_allocator(rmm.rmm_cupy_allocator)



def step2dPredictor(compTimes, GRID, OCEAN, BOUNDARY):
    from mod_operators import grsz, bksz

    # Aliases
    zeta_t2, zeta_t1, zeta_t0 = (OCEAN.zeta_t2, OCEAN.zeta_t1, OCEAN.zeta_t0)
    ubar_t2, ubar_t1, ubar_t0 = (OCEAN.ubar_t2, OCEAN.ubar_t1, OCEAN.ubar_t0)
    vbar_t2, vbar_t1, vbar_t0 = (OCEAN.vbar_t2, OCEAN.vbar_t1, OCEAN.vbar_t0)
    rzeta_t1, rzeta_t0 = (OCEAN.rzeta_t1, OCEAN.rzeta_t0)
    rubar_t1, rubar_t0 = (OCEAN.rubar_t1, OCEAN.rubar_t0)
    rvbar_t1, rvbar_t0 = (OCEAN.rvbar_t1, OCEAN.rvbar_t0)
    h = GRID.h.ravel()

    Δt = compTimes.dtfast


    # Free-surface equation.
    # =================================

    # During the first time-step, the predictor step is Forward-Euler. Otherwise, the predictor step is Leap-frog.
    computeZetaRHS((grsz[0]*GPUMUL,), (bksz[0]//GPUMUL,), (zeta_t1, h, ubar_t1, vbar_t1, OCEAN.DU_avg1, OCEAN.DV_avg1, OCEAN.Zt_avg1, True, compTimes.weight1, compTimes.weight2, compTimes.iif, rzeta_t1))

    # print('XXXXXX', bksz[0]//GPUMUL, grsz[0]*GPUMUL)

    if compTimes.isFirst2DStep():
        # The first time it performs a simple Euler time step. RHS is computed at tn and time derivatives are
        # (f(tn) - f(tn-1))/dtfast

        zeta_t2[:] = zeta_t1 + Δt*rzeta_t1

        wGzeta = 0.5   # This is to compute gzeta (see computeMomentumRHSPred)

    else:
        # If is not the first time step, the predictor consists of a leapfrog time step where the RHS is computed at tn
        # and time derivatives are centered at tn (f(tn+1) - f(tn-1))/2*dtfast.
        computeZetaPred(grsz, bksz, (Δt, zeta_t0, zeta_t1, zeta_t2, rzeta_t1))

        wGzeta = 2.0 / 5.0  # This is to compute gzeta (see computeMomentumRHSPred)


    # Apply free-surface lateral BC
    zetabc(zeta_t2, compTimes, BOUNDARY)


    # Momentum equations.
    # ===================

    #compute right-hand-side for the 2D momentum equations
    computeMomentumRHSPred(grsz, bksz, (h, rubar_t1, rvbar_t1, zeta_t0, zeta_t2, g, wGzeta))


    # During the first time-step, the predictor step is Forward-Euler
    # and the corrector step is Backward-Euler. Otherwise, the predictor
    # step is Leap-frog and the corrector step is Adams-Moulton.
    computeMomentumPred(grsz, bksz, (Δt, ubar_t1, ubar_t2, vbar_t1, vbar_t2, rubar_t1, rvbar_t1, h, zeta_t1, zeta_t2))


    # Apply lateral boundary conditions.
    barotropicVelocityBC(ubar_t2, vbar_t2, compTimes, BOUNDARY)



def step2dCorrector(compTimes, GRID, OCEAN, BOUNDARY):
    from mod_operators import grsz, bksz

    # Aliases
    zeta_t2,  zeta_t1,  zeta_t0  = (OCEAN.zeta_t2,  OCEAN.zeta_t1,  OCEAN.zeta_t0)
    ubar_t2,  ubar_t1,  ubar_t0  = (OCEAN.ubar_t2,  OCEAN.ubar_t1,  OCEAN.ubar_t0)
    vbar_t2,  vbar_t1,  vbar_t0  = (OCEAN.vbar_t2,  OCEAN.vbar_t1,  OCEAN.vbar_t0)
    rzeta_t2, rzeta_t1, rzeta_t0 = (OCEAN.rzeta_t2, OCEAN.rzeta_t1, OCEAN.rzeta_t0)
    rubar_t2, rubar_t1, rubar_t0 = (OCEAN.rubar_t2, OCEAN.rubar_t1, OCEAN.rubar_t0)
    rvbar_t2, rvbar_t1, rvbar_t0 = (OCEAN.rvbar_t2, OCEAN.rvbar_t1, OCEAN.rvbar_t0)


    h = GRID.h.ravel()

    Δt = compTimes.dtfast

    # Free-surface equation.
    # =================================

    computeZetaRHS((grsz[0]*GPUMUL,), (bksz[0]//GPUMUL,), (zeta_t2, h, ubar_t2, vbar_t2, OCEAN.DU_avg1, OCEAN.DV_avg1, OCEAN.Zt_avg1, False,
                                                 compTimes.weight1, compTimes.weight2, compTimes.iif, rzeta_t1))



    # Adams-Moulton order 3
    AdamsMoultonCorr3rd(grsz, bksz, (Δt, zeta_t2, rzeta_t0, rzeta_t1, rzeta_t2))

    # Apply free-surface lateral BC
    zetabc(zeta_t2, compTimes, BOUNDARY)


    #compute right-hand-side for the 2D momentum equations
    computeMomentumRHSCorr(grsz, bksz, (h, rubar_t2, rvbar_t2, zeta_t0, zeta_t1, zeta_t2, g))


    # During the first time-step, the corrector step is Backward-Euler. Otherwise, the corrector step is Adams-Moulton.
    AdamsMoultonCorr3rd2(grsz, bksz, (Δt, ubar_t1, vbar_t1,
                                      rubar_t0, rubar_t1, rubar_t2,
                                      rvbar_t0, rvbar_t1, rvbar_t2,
                                      h, zeta_t1, zeta_t2))


    # Apply lateral boundary conditions.
    barotropicVelocityBC(ubar_t2, vbar_t2, compTimes, BOUNDARY)



     # After all fast time steps are completed, apply boundary conditions to time averaged fields.






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


