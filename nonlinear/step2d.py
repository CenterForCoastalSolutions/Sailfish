

def step2d(compTimes, GRID, OCEAN):


    ptsk = 3 - compTimes.kstp
    is2DCorrectorStep = not compTimes.PREDICTOR_2D_STEP

    shp = Slice(JstrVm2-1:Jendp2, IstrUm2-1:Iendp2)


    zeta_trhs = OCEAN.zeta[compTimes.krhs,:,:].ravel()
    ubar_trhs = OCEAN.ubar[compTimes.krhs,:,:].ravel()
    vbar_trhs = OCEAN.vbar[compTimes.krhs,:,:].ravel()
    zeta_tstp = OCEAN.zeta[compTimes.kstp,:,:].ravel()
    ubar_tstp = OCEAN.ubar[compTimes.kstp,:,:].ravel()
    vbar_tstp = OCEAN.vbar[compTimes.kstp,:,:].ravel()
    on_u = GRID.on_u[shp]
    om_v = GRID.om_v[shp]

    Drhs = zeta_rhs[shp] + h[shp]


    DUon[shp] = ubar_rhs[shp]*on_u*RtoU(Drhs[shp])
    DVom[shp] = vbar_rhs[shp]*om_v*RtoV(Drhs[shp])

    know, Δt = compTimes.get2DTimes()
    kstp = compTimes.kstp
    krhs = compTimes.krhs

    if compTimes.isFirst2DStep() :
        # The first time it performs a simple Euler time step. RHS is computed at tn and time derivatives are (f(tn) - f(tn-1))/dtfast


        rhs_zeta[:,:] = DivUVtoR(DUon[shp], DVom[shp], GRID)
        zeta_new[:,:] = zeta_tstp + Δt*rhs_zeta

        #if def masking apply mask (this is completed way earlier)
        Dnew[:,:]     = zeta_new + h
        zwrk[:,:]     = 0.5*zeta_tstp + 0.5*zeta_new[:,:])

        gzeta[:,:]    = zwrk
        gzeta2[:,:]   = zwrk*zwrk

    elif compTimes.predictor2DStep:
        # The predictor consists of a leapfrog time step where the RHS is computed at tn and time derivatives are centered at tn (f(tn+1) - f(tn-1))/2*dtfast.

        cff4          = 4.0/25.0
        cff5          = 1.0 - 2.0*cff4

        rhs_zeta[:,:] = DivUVtoR(DUon[shp], DVom[shp], GRID)
        zeta_new[:,:] = zeta_tstp + Δt*rhs_zeta

        Dnew[:,:]     = zeta_new + h
        zwrk[:,:]     = cff5*zeta_trhs + cff4*(zeta_tstp + zeta_new)

        gzeta[:,:]    = zwrk
        gzeta2[:,:]   = zwrk*zwrk

    elif compTimes.corrector2DStep:

        cff1          = 5.0/12.0
        cff2          = 8.0/12.0
        cff3          = 1.0/12.0
        cff4          = 2.0/5.0
        cff5          = 1.0 - cff4
        rzeta1        = DivUVtoR(DUon[shp], DVom[shp], GRID)
        zeta_new[:,:] = zeta_tstp + Δt*(cff1*rzeta1 + cff2*rzeta[kstp,:,:] - cff3*rzeta[ptsk,:,:])  # REMEMBER WE CHANGED rzeta definition
        Dnew[:,:]     = zeta_new + h
        zwrk[:,:]     = cff5*zeta_new + cff4*zeta_trhs

    else:
        msgError('This code shouldn''t be reachable')


    zeta[knew,:,:] = zeta_new


    if PREDICTOR_2D_STEP:
        rzeta[khrs,:,:] = rhs_zeta

    #Apply mass point sources (volume vertical influx), if any
    if LwSrc:
        # i = SOURCES.Isrx
        # j = SOURCES.Jsrc
        ij = SOURCES.IJsrc   # This has to be precomputed as i + j*width
        zeta[knew,:,:].ravel()[ij] += SOURCES.Qbar*pm.ravel()[ij]*pn.ravel()[ij]*dtfast


    zetabc(krhs, kstp, knew, zeta)

    #compute right-hand-side for the 2D momentum equations
    cff1 = 0.5 * g
    cff2 = 1.0 / 3.0

    #compute pressure gradient terms
    rhs_ubar[:,:] = cff1*on_u*(RtoU(h)*DxU(gzeta + gzeta2))

    rhs_vbar[:,:] = cff1*om_v*(RtoV(h)*DyV(gzeta + gzeta2))

    # Add in horizontal advection of momentum.

    # # second-order, centered differences advection.
    # if UV_C2ADVECTION:
    #     UFξ[:,:] = UtoR(DUon)*UtoR(ubar[krhs,:,:])
    #     UFη[:,:] = UtoR(DVom)*UtoR(ubar[krhs,:,:])
    #     VFξ[:,:] = UtoR(DUon)*VtoR(vbar[krhs,:,:])
    #     VFη[:,:] = UtoR(DVom)*VtoR(vbar[krhs,:,:])

    rhs_ubar[:,:] -= (DxR1(UFξ) + DyR2(UFη))
    rhs_vbar[:,:] -= (DyR1(VFx) + DxR2(VFe))


# ****************************
# ****************************
# ****************************

    # LINE 1664 .... if defined UV_VIS2 || defined UV_VIS4

    if UV_ADV:
        #!---------------------------------------------------------------------------
        #! Contribution of a term corresponding to product of
        #! Stokes and Eulerian Velocity Eqn. 26 and 27.
        #! This removes terms that were unneccessarily added in flux form.
        #!---------------------------------------------------------------------------
        cff         = 0.5 * (Drhs(i-1,j) + Drhs(i,j))
        DUSon(i,j)  = cff * on_u(i, j) * ubar_stokes(i, j)
        DVSon(i, j) = 0.25 * cff * on_u(i, j) * (vbar_stokes(i  ,j  ) + vbar_stokes(i  ,j+1) + vbar_stokes(i-1,j  ) + vbar_stokes(i-1,j+1))

        UFx(i,j)    = 0.5 * (ubar(i, j, krhs) + ubar(i+1, j, krhs))
        VFx(i,j)    = 0.5 * (vbar(i, j, krhs) + vbar(i, j+1, krhs))

        cff         = 0.5 * (Drhs(i,j) + Drhs(i,j-1))
        DUSom(i,j)  = cff * 0.25 * om_v(i,j) * (ubar_stokes(i,j) + ubar_stokes(i+1, j) + ubar_stokes(i  , j-1) + ubar_stokes(i+1, j-1))

        DVSom(i,j)    = cff * om_v(i, j) * vbar_stokes(i, j)
        cff           = 0.5 * (Drhs(i, j) + Drhs(i, j-1))
        UFe(i,j)      = 0.5 * (ubar(i+1, j, krhs) + ubar(i, j, krhs))
        VFe(i,j)      = 0.5 * (vbar(i, j, krhs) + vbar(i, j+1, krhs))
        cff1          = UFx(i,j) - UFx(i-1,j) # line 2215
        cff2          = VFx(i,j) - VFx(i-1,j)
        cff3          = DUSon(i,j) * cff1
        cff4          = DVSon(i,j) * cff2
        rhs_ubar(i,j) = rhs_ubar(i,j) + cff3 + cff4

        cff1          = UFe(i,j) - UFe(i,j-1)
        cff2          = VFe(i,j) - VFe(i,j-1)
        cff3          = DUSom(i,j) * cff1
        cff4          = DVSom(i,j) * cff2
        rhs_vbar(i,j) = rhs_vbar(i,j) + cff3 + cff4

        # ifndef SOLVE3D line 2245
        # add in bottom stress

        fac           = bustr(i,j) * om_u(i,j) * on_u(i,j)
        rhs_ubar(i,j) = rhs_ubar(i,j) - fac

        fac           = bvstr(i,j) * om_v(i,j) * on_v(i,j)
        rhs_vbar(i,j) = rhs_vbar(i,j) - fac

        # !  Time step 2D momentum equations. LINE 2627

        # compute the water column depth
        Dstp(i,j)     = zeta(i,j,kstp) + h(i,j)

        # !  During the first time-step, the predictor step is Forward-Euler
        # !  and the corrector step is Backward-Euler. Otherwise, the predictor
        # !  step is Leap-frog and the corrector step is Adams-Moulton.

        if FIRST_2D_STEP:

            cff1           = 0.5 * dtfast
            cff            = (pm(i,j) + pm(i-1,j)) * (pn(i,j) + pn(i-1,j))
            fac            = 1.0 / (Dnew(i,j) + Dnew(i-1,j))
            ubar(i,j,knew) = (ubar(i,j,kstp) * (Dstp(i,j) + Dstp(i-1,j)) + cff * cff1 * rhs_ubar(i,j)) * fac
            cff            = (pm(i,j) + pm(i,j-1)) * (pn(i,j) + pn(i,j-1))
            fac            = 1.0 / (Dnew(i,j) + Dnew(i,j-1))
            vbar(i,j,knew) = (vbar(i,j,kstp) * (Dstp(i,j) + Dstp(i,j-1)) + cff * cff1 * rhs_vbar(i,j)) * fac

        else if PREDICTOR_2D_STEP:
            cff1           = dtfast
            cff            = (pm(i,j) + pm(i-1,j)) * (pn(i,j) + pn(i-1,j))
            fac            = 1.0 / (Dnew(i,j) + Dnew(i-1,j))
            ubar(i,j,knew) = (ubar(i,j,kstp) * (Dstp(i,j) + Dstp(i-1,j)) + cff * cff1 * rhs_ubar(i,j)) * fac
            cff            = (pm(i,j) + pm(i,j-1)) * (pn(i,j) + pn(i,j-1))
            fac            = 1.0 / (Dnew(i,j) + Dnew(i,j-1))
            vbar(i,j,knew) = (vbar(i,j,kstp) * (Dstp(i,j) + Dstp(i,j-1)) + cff * cff1 * rhs_vbar(i,j)) * fac

        else if CORRECTOR_2D_STEP:
            cff1=0.5*dtfast*5.0/12.0
            cff2=0.5*dtfast*8.0/12.0
            cff3=0.5*dtfast*1.0/12.0
            cff=(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
            fac=1.0_r8/(Dnew(i,j)+Dnew(i-1,j))
            ubar(i,j,knew)=(ubar(i,j,kstp)*(Dstp(i,j)+Dstp(i-1,j))+cff*(cff1*rhs_ubar(i,j)+cff2*rubar(i,j,kstp)-cff3*rubar(i,j,ptsk)))*fac

            cff=(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
            fac=1.0_r8/(Dnew(i,j)+Dnew(i,j-1))
            vbar(i,j,knew)=(vbar(i,j,kstp)*(Dstp(i,j)+Dstp(i,j-1))+cff*(cff1*rhs_vbar(i,j)+cff2*rvbar(i,j,kstp)-cff3*rvbar(i,j,ptsk)))*fac

            # LINE 2982

        #             !  Apply lateral boundary conditions.
        # !-----------------------------------------------------------------------
        # !
        #       CALL u2dbc (LBi, UBi, LBj, UBj,                              &
        #      &                 IminS, ImaxS, JminS, JmaxS,                      &
        #      &                 krhs, kstp, knew,                                &
        #      &                 ubar, vbar, zeta)
        #       CALL v2dbc (LBi, UBi, LBj, UBj,                              &
        #      &                 IminS, ImaxS, JminS, JmaxS,                      &
        #      &                 krhs, kstp, knew,                                &
        #      &                 ubar, vbar, zeta)
        # !
        # !  Compute integral mass flux across open boundaries and adjust
        # !  for volume conservation.
        # !
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


