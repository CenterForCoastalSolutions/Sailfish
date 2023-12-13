import mod_comptimes
from mod_operators import *
from barotropicVelocityBC import barotropicVelocityBC
from zetabc               import zetabc
from mod_constants import *
from misc          import *
import time


# filename = os.path.join(exePath, r'nonlinear/step2d_kernels.cpp')
# filename = os.path.join(exePath, r'modules/mod_cppkernels.cpp')
# with open(filename, 'r') as file:
#     code = file.read()
# moduleCPPKernels = cp.RawModule(code=code, options=('-default-device',  '--restrict', '--std=c++17'))
#
# computeMomentumRHSCorr = moduleCPPKernels.get_function('computeMomentumRHSCorr')
# computeMomentumRHSPred = moduleCPPKernels.get_function('computeMomentumRHSPred')
# computeZetaRHS3     = moduleCPPKernels.get_function('computeZetaRHS')
# aaa     = moduleCPPKernels.get_function('aaa')
# bbb     = moduleCPPKernels.get_function('bbb')
# AdamsMoultonCorr3rd = moduleCPPKernels.get_function('AdamsMoultonCorr3rd')
# AdamsMoultonCorr3rd2 = moduleCPPKernels.get_function('AdamsMoultonCorr3rd')
# Pred  = moduleCPPKernels.get_function('Pred')
# Pred2 = moduleCPPKernels.get_function('Pred2')

#
# import rmm
# pool = rmm.mr.PoolMemoryResource(
#     rmm.mr.ManagedMemoryResource(),
#     initial_pool_size=5*(2**30),
#     maximum_pool_size=5*(2**30)
# )
# rmm.mr.set_current_device_resource(pool)
# cp.cuda.set_allocator(rmm.rmm_cupy_allocator)

# In this module, t2, t1 and t0 refer to

# cp.cuda.set_allocator(cp.cuda.MemoryPool().malloc)

# @cp.fuse
# import numba.cuda
# @numba.cuda.jit('float64(float64, float64, float64, float64)', device=True, inline=True)
# def computeZetaRHS(zeta, h, ubar, vbar):
#
#     #Apply mass point sources (volume vertical influx), if any
#     # TODO: Implement XXXX
#     # if LwSrc:
#     #     ij = SOURCES.IJsrc   # This has to be precomputed as i + j*width
#     #     zeta[knew,:,:].ravel()[ij] += SOURCES.Qbar*pm.ravel()[ij]*pn.ravel()[ij]*dtfast
#     # addForcings()
#
#     # compute the water column depth
#     D = zeta + h
#
#     DU = ubar*RtoU(D)
#     DV = vbar*RtoV(D)
#
#     return divUVtoR(DU, DV)


# computeMomentumRHS = cp.ElementwiseKernel(
#     'float64 h, float64 gzeta',
#     'float64 rx, float64 ry',
#     'z = (x - y) * (x - y)',
#     'computeMomentumRHS')




# def computeMomentumRHS2(h, gzeta):
#     rhs_ubar = 0.5*g*(RtoU(h)*DξRtoU(gzeta) + DξRtoU(gzeta*gzeta))
#     rhs_vbar = 0.5*g*(RtoV(h)*DξRtoU(gzeta) + DξRtoU(gzeta*gzeta))


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

    # return rhs_ubar, rhs_vbar



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
    # rhs_zeta_t1 = computeZetaRHS(zeta_t1, h, ubar_t1, vbar_t1)
    computeZetaRHS(grsz, bksz, (zeta_t1, h, ubar_t1, vbar_t1, rzeta_t1))

    if compTimes.isFirst2DStep():
        # The first time it performs a simple Euler time step. RHS is computed at tn and time derivatives are
        # (f(tn) - f(tn-1))/dtfast

        zeta_t2[:] = zeta_t1 + Δt*rzeta_t1

        weight = 0.5   # This is to compute gzeta (see computeMomentumRHSPred)

    else:
        # If is not the first time step, the predictor consists of a leapfrog time step where the RHS is computed at tn
        # and time derivatives are centered at tn (f(tn+1) - f(tn-1))/2*dtfast.
        computeZetaPred(grsz, bksz, (Δt, zeta_t0, zeta_t1, zeta_t2, rzeta_t1))

        weight = 2.0 / 5.0  # This is to compute gzeta (see computeMomentumRHSPred)


    # Apply free-surface lateral BC
    zetabc(zeta_t2, compTimes, BOUNDARY)


    # Momentum equations.
    # ===================


    #compute right-hand-side for the 2D momentum equations
    computeMomentumRHSPred(grsz, bksz, (h, rubar_t1, rvbar_t1, zeta_t1, zeta_t2, g, weight))


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

    computeZetaRHS(grsz, bksz, (zeta_t2, h, ubar_t2, vbar_t2, rzeta_t1))

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

    # TODO: REMEMBER XXXXX
    # YYYY
    # print("uuuuuu", OCEAN.DUon.data, OCEAN.DVom.data, OCEAN.ubar_t2.data, OCEAN.vbar_t2.data, GRID.h.ravel().data)
    cp.cuda.runtime.deviceSynchronize()
    print("uuuuuu", OCEAN.ubar_t2.data, OCEAN.vbar_t2.data)
    OCEAN.DUon[:] = OCEAN.ubar_t2*GRID.h.ravel()
    OCEAN.DVom[:] = OCEAN.vbar_t2*GRID.h.ravel()
    # OCEAN.DU_avg2[:] = OCEAN.DUon
    # OCEAN.DV_avg2[:] = OCEAN.DVom
    OCEAN.DU_avg1[:] = OCEAN.DUon
    OCEAN.DV_avg1[:] = OCEAN.DVom



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


