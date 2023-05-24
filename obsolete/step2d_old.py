import mod_comptimes
from mod_operators import *
from barotropicVelocityBC import barotropicVelocityBC
from zetabc               import zetabc
from mod_constants import *
from misc          import *

def step2d(isPredictorStep, compTimes, GRID, OCEAN):


    # Do not perform the actual time stepping during the auxiliary nfast(ng)+1) time step.
    if compTimes.iif > compTimes.nfast:
        return

    ptsk = 3 - compTimes.kstp   # TODO: Do this in compTimes, check if this variable has another name knew?


    # Aliases/views
    zeta_t2  = OCEAN.zeta [compTimes.i2,:,:].ravel()
    ubar_t2  = OCEAN.ubar [compTimes.i2,:,:].ravel()
    vbar_t2  = OCEAN.vbar [compTimes.i2,:,:].ravel()
    rzeta_t2 = OCEAN.rzeta[compTimes.i2,:,:].ravel()
    rubar_t2 = OCEAN.rubar[compTimes.i2,:,:].ravel()
    rvbar_t2 = OCEAN.rvbar[compTimes.i2,:,:].ravel()
    zeta_t1  = OCEAN.zeta [compTimes.i1,:,:].ravel()
    ubar_t1  = OCEAN.ubar [compTimes.i1,:,:].ravel()
    vbar_t1  = OCEAN.vbar [compTimes.i1,:,:].ravel()
    rzeta_t1 = OCEAN.rzeta[compTimes.i1,:,:].ravel()
    rubar_t1 = OCEAN.rubar[compTimes.i1,:,:].ravel()
    rvbar_t1 = OCEAN.rvbar[compTimes.i1,:,:].ravel()
    rzeta_t0 = OCEAN.rzeta[compTimes.i0,:,:].ravel()
    rubar_t0 = OCEAN.rubar[compTimes.i0,:,:].ravel()
    rvbar_t0 = OCEAN.rvbar[compTimes.i0,:,:].ravel()


    h = GRID.h

    D_i1 = zeta_t1 + h

    # DUon can be read as "DU over n", in other words: D*U/n and D*V/m
    DU_t1 = ubar_t1*RtoU(D_i1)
    DV_t1 = vbar_t1*RtoV(D_i1)


    know, Δt = compTimes.get2DTimes()
    kstp, krhs = compTimes.updateTimes(isPredictorStep) xxxx

    # Time step free-surface equation.
    # =================================

    # During the first time-step, the predictor step is Forward-Euler and the corrector step is Backward-Euler. Otherwise,
    # the predictor step is Leap-frog and the corrector step is Adams-Moulton.

    rhs_zeta = DivUVtoR(DU, DV)

    if compTimes.isFirst2DStep():
        # The first time it performs a simple Euler time step. RHS is computed at tn and time derivatives are (f(tn) - f(tn-1))/dtfast

        zeta_t2[:] = zeta_t1+ Δt*rhs_zeta

        gzeta = 0.5*(zeta_t2 + zeta_t1)


    elif isPredictorStep:
        # The predictor consists of a leapfrog time step where the RHS is computed at tn and time derivatives are centered at tn (f(tn+1) - f(tn-1))/2*dtfast.

        zeta_t2[:] = zeta_t1 + Δt*rhs_zeta

        w = 4.0/25.0
        gzeta = (1 - w)*zeta_t1 + w*0.5*(zeta_tstp + zeta_tnew)

        # In the predictor step, save the rhs for future use.
        rzeta_t2[:] = rhs_zeta


    else:  # Corrector step

        zeta_t2[:] = zeta_t1 + Δt*(AM3_2*rzeta_t2 + AM3_2*rzeta_t1 + AM3_2*rzeta_t0)

        w = 2.0/5.0
        gzeta  = (1 - w)*zeta_t2 + w*zeta_t1



    gzeta2 = gzeta*gzeta



    #Apply mass point sources (volume vertical influx), if any
    # TODO: Implement XXXX
    # if LwSrc:
    #     ij = SOURCES.IJsrc   # This has to be precomputed as i + j*width
    #     zeta[knew,:,:].ravel()[ij] += SOURCES.Qbar*pm.ravel()[ij]*pn.ravel()[ij]*dtfast



    # Apply free-surface lateral BC
    zetabc(OCEAN.zeta, kout, compTimes, BOUNDARY)


    #compute right-hand-side for the 2D momentum equations
    rhs_ubar = 0.5*g*RtoU(h)*DξU(gzeta + gzeta2)
    rhs_vbar = 0.5*g*RtoV(h)*DηV(gzeta + gzeta2)




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
    Dstp = zeta_tstp + h
    Dnew = zeta_tnew + h

    # And interpolate them at points U, V
    DstpU = RtoU(Dstp)
    DstpV = RtoV(Dstp)
    DnewU = RtoU(Dnew)
    DnewV = RtoV(Dnew)


    # Metric factors interpolated at U, V points.
    pmnU = RtoU(pm)*RtoU(pn)
    pmnV = RtoV(pm)*RtoV(pn)


    # During the first time-step, the predictor step is Forward-Euler
    # and the corrector step is Backward-Euler. Otherwise, the predictor
    # step is Leap-frog and the corrector step is Adams-Moulton.


    if compTimes.isFirst2DStep() or isPredictorStep:

        ubar_tnew[:] = (ubar_tstp*(DstpU + Δt*pmnU*rhs_ubar))/DnewU
        vbar_tnew[:] = (vbar_tstp*(DstpV + Δt*pmnV*rhs_vbar))/DnewV

        # In the predictor step, save the rhs for future use.
        rubar_trhs[:] = rhs_ubar
        rvbar_trhs[:] = rhs_vbar


    else:

        ubar_tnew[:] = (ubar_tstp*(DstpU + Δt*pmnU*(AM3a*rhs_ubar + AM3b*rubar_tstp + aM3c*rubar_ptsk)))/DnewU
        vbar_tnew[:] = (vbar_tstp*(DstpV + Δt*pmnV*(AM3a*rhs_vbar + AM3b*rvbar_tstp + aM3c*rvbar_ptsk)))/DnewU






    # Apply lateral boundary conditions.
    barotropicVelocityBC(ubar, vbar, zeta, krhs, kstp, knew)


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


