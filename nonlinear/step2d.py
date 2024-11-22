import mod_comptimes
from mod_operators import *
from barotropicVelocityBC import barotropicVelocityBC
from zetabc               import zetabc
from mod_constants import *
from misc          import *
from numba import njit
import time




# import rmm
# pool = rmm.mr.PoolMemoryResource(
#     rmm.mr.ManagedMemoryResource(),
#     initial_pool_size=55*(2**30),
#     maximum_pool_size=55*(2**30)
# )
# rmm.mr.set_current_device_resource(pool)
# cp.cuda.set_allocator(rmm.rmm_cupy_allocator)


@njit
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



    # print('time step: ', compTimes.iic)

    # Free-surface equation.
    # =================================

    # During the first time-step, the predictor step is Forward-Euler. Otherwise, the predictor step is Leap-frog.
    computeZetaRHS((grsz[0]*GPUMUL,), (bksz[0]//GPUMUL,), (zeta_t1, h, ubar_t1, vbar_t1, OCEAN.DU_avg1, OCEAN.DV_avg1, OCEAN.Zt_avg1, True, compTimes.weight1, compTimes.weight2, compTimes.iif, rzeta_t1))

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
    # TODO: I don't think this comment is correct.

    computeMomentumPred(grsz, bksz, (Δt, ubar_t1, ubar_t2, vbar_t1, vbar_t2, rubar_t1, rvbar_t1, h, zeta_t1, zeta_t2))


    # Apply lateral boundary conditions.
    barotropicVelocityBC(ubar_t2, vbar_t2, compTimes, BOUNDARY)



def step2dCorrector(compTimes, GRID, OCEAN, BOUNDARY):
    from mod_operators import grsz, bksz

    # Aliases
    zeta_t2, zeta_t1, zeta_t0    = (OCEAN.zeta_t2, OCEAN.zeta_t1, OCEAN.zeta_t0)
    ubar_t2, ubar_t1, ubar_t0    = (OCEAN.ubar_t2, OCEAN.ubar_t1, OCEAN.ubar_t0)
    vbar_t2, vbar_t1, vbar_t0    = (OCEAN.vbar_t2, OCEAN.vbar_t1, OCEAN.vbar_t0)
    rzeta_t2, rzeta_t1, rzeta_t0 = (OCEAN.rzeta_t2, OCEAN.rzeta_t1, OCEAN.rzeta_t0)
    rubar_t2, rubar_t1, rubar_t0 = (OCEAN.rubar_t2, OCEAN.rubar_t1, OCEAN.rubar_t0)
    rvbar_t2, rvbar_t1, rvbar_t0 = (OCEAN.rvbar_t2, OCEAN.rvbar_t1, OCEAN.rvbar_t0)

    # zeta_t2, zeta_t1, zeta_t0, ubar_t2, ubar_t1, ubar_t0, vbar_t2, vbar_t1, vbar_t0, \
    # rzeta_t1, rzeta_t0, rubar_t1, rubar_t0, rvbar_t1, rvbar_t0 = OCEAN.getVars()
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


    # DUon = ubar(i,j,krhs)*on_u*RtoU(zeta + h)

    # # TODO: REMEMBER XXXXX
    # # cp.cuda.runtime.deviceSynchronize()
    # OCEAN.DUon[:] = OCEAN.ubar_t2*h
    # OCEAN.DVom[:] = OCEAN.vbar_t2*h
    # OCEAN.DU_avg1[:] = OCEAN.DUon
    # OCEAN.DV_avg1[:] = OCEAN.DVom

    # TODO: REMEMBER XXXXX
    # cff1 = weight(1,iif-1)
    # cff2 = (8.0/12.0)*weight(2,iif) - (1.0/12.0)*weight(2,iif+1)
    #
    # Zt_avg1 += cff1*zeta(krhs)
    #
    # DU_avg1 += cff1*DUon(i,j)
    # DU_avg2 += cff2*DUon(i,j)



  #   # Compute time averaged fields over all short time-steps.
  #   IF (PREDICTOR_2D_STEP(ng)) THEN
  #       IF (FIRST_2D_STEP) THEN
  #
  #       # Reset arrays for 2D fields averaged within the short time-steps.
  #       cff2=(-1.0_r8/12.0_r8)*weight(2,iif(ng)+1,ng)
  #             DO j=JstrR,JendR
  #               DO i=IstrR,IendR
  #                 Zt_avg1(i,j)=0.0_r8
  #               END DO
  #               DO i=Istr,IendR
  #                 DU_avg1(i,j)=0.0_r8
  #                 DU_avg2(i,j)=cff2*DUon(i,j)
  #               END DO
  #             END DO
  #
  #             DO j=Jstr,JendR
  #               DO i=IstrR,IendR
  #                 DV_avg1(i,j)=0.0_r8
  #                 DV_avg2(i,j)=cff2*DVom(i,j)
  #               END DO
  #             END DO
  #       ELSE
  #
  # # Accumulate field averages of previous time-step after they are
  # # computed in the previous corrector step, updated their boundaries,
  # # and synchronized.
  #
  #         cff1=weight(1,iif(ng)-1,ng)
  #         cff2=(8.0_r8/12.0_r8)*weight(2,iif(ng)  ,ng)-                 &
  #    &         (1.0_r8/12.0_r8)*weight(2,iif(ng)+1,ng)
  #         DO j=JstrR,JendR
  #           DO i=IstrR,IendR
  #             Zt_avg1(i,j)=Zt_avg1(i,j)+cff1*zeta(i,j,krhs)
  #           END DO
  #           DO i=Istr,IendR
  #             DU_avg1(i,j)=DU_avg1(i,j)+cff1*DUon(i,j)
  #             DU_avg2(i,j)=DU_avg2(i,j)+cff2*DUon(i,j)
  #           END DO
  #         END DO
  #         DO j=Jstr,JendR
  #           DO i=IstrR,IendR
  #             DV_avg1(i,j)=DV_avg1(i,j)+cff1*DVom(i,j)
  #
  #             DV_avg2(i,j)=DV_avg2(i,j)+cff2*DVom(i,j)
  #           END DO
  #         END DO
  #
  #       END IF
  #     ELSE
  #
  #       IF (FIRST_2D_STEP) THEN
  #         cff2=weight(2,iif(ng),ng)
  #       ELSE
  #         cff2=(5.0_r8/12.0_r8)*weight(2,iif(ng),ng)
  #       END IF
  #
  #
  #       DO j=JstrR,JendR
  #         DO i=Istr,IendR
  #           DU_avg2(i,j)=DU_avg2(i,j)+cff2*DUon(i,j)
  #         END DO
  #       END DO
  #       DO j=Jstr,JendR
  #         DO i=IstrR,IendR
  #           DV_avg2(i,j)=DV_avg2(i,j)+cff2*DVom(i,j)
  #         END DO
  #       END DO
  #     END IF


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


