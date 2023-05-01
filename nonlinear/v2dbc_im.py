
# !  This subroutine sets lateral boundary conditions for vertically     !
# !  integrated V-velocity.                                              !

      USE mod_param
      USE mod_ocean
      USE mod_stepping

      CALL v2dbc_tile (ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
     &                 krhs(ng), kstp(ng), kout,                        &
     &                 OCEAN(ng) % ubar,                                &
     &                 OCEAN(ng) % vbar,                                &
     &                 OCEAN(ng) % zeta)

      RETURN
      END SUBROUTINE v2dbc

      USE mod_param
      USE mod_boundary
      USE mod_clima
      USE mod_forces
      USE mod_grid
      USE mod_ncparam
#ifdef NESTING
      USE mod_nesting
#endif
#ifdef WEC
      USE mod_ocean
#endif
      USE mod_scalars

# !  Local variable declarations.

      integer :: Jmin, Jmax
      integer :: i, j, kNow
#ifdef NESTING
      integer :: Idg, Jdg, cr, dg, m, rg, tnew, told
#endif

      real(r8), parameter :: eps = 1.0E-20_r8

      real(r8) :: Ce, Cx, Ze
      real(r8) :: bry_pgr, bry_cor, bry_str, bry_wec, bry_val
      real(r8) :: cff, cff1, cff2, cff3, dt2d, dVde, dVdt, dVdx
      real(r8) :: obc_in, obc_out, phi, tau
#ifdef WEC_VF
      real(r8), parameter :: Lwave_min = 1.0_r8
      real(r8) :: sigma, osigma, waven, waveny
#endif
#if defined ATM_PRESS && defined PRESS_COMPENSATE
      real(r8) :: OneAtm, fac
#endif

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: grad

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set time-indices
!-----------------------------------------------------------------------
!

    domain   = DOMAIN[ng]
    grid     = GRID[ng]
    boundary = BOUNDARY[ng]
    # lbc      = lbc[ng]
    forces   = FORCES[ng]






    if FIRST_2D_STEP:
        kNow = krhs
        dt2d = dtfast[ng]

    elif PREDICTOR_2D_STEP[ng]:
        kNow=krhs
        dt2d = 2.0*dtfast[ng]

    else:
        kNow = kstp
        dt2d = dtfast[ng]


    if ATM_PRESS and PRESS_COMPENSATE:
        OneAtm = 1013.25
        fac = 100.0/(g*rho0)




# closed boundary condition.
vbar[kout, idxClosedBC] = 0.0



# Clamped boundary condition.
vbar[kout, idxClampedBC] = boundary.***vbar_south[idxBCClamped]


# Gradient boundary condition.
vbar[kout, idxGradientBC] = vbar[kout, idxGradientBC2]



##### Reduced-physics boundary condition.
##### Reduced-physics boundary condition.
##### Reduced-physics boundary condition.


# Computes the gradient at the boundary

# NO IFDEF HERE???
if (LBC[isouth, isFsur, ng].acquire):  # isFreeSurface
    bry_pgr = -g * (zeta[idxReducedPhysicsBC, know] - boundary.zeta_south[i])*0.5*grid.pn[idxReducedPhysicsBC]
else:
    bry_pgr = -g * (zeta(idxReducedPhysicsBC, know) - zeta(idxReducedPhysicsBC2, kNow))*0.5*RtoV(grid.pn, idxReducedPhysicsBC)

# ifdef UV_COR
    bry_cor = UtoV(ubar, idxReducedPhysicsBC) * RtoV(grid.f, idxReducedPhysicsBC)
# else
    bry_cor = 0.0
# endif

    zhR = RtoV(grid.h, idx) + RtoV(grid.zeta[:, :, know], idx)


    bry_str = (forces.svstr[idxReducedPhysicsBC] - forces.bvstr[idxReducedPhysicsBC]) / zhR
    vbar[idxReducedPhysicsBC, kout] = vbar(idxReducedPhysicsBC, kNow) + dt2d * (bry_pgr + bry_cor + bry_str)





########### FLATHER
########### FLATHER
########### FLATHER


        # elif (LBC(isouth,isVbar,ng)%Flather):
        #   DO i=Istr,Iend
        #     IF (LBC_apply(ng)%south(i)) THEN

    if (defined SSH_TIDES && !defined UV_TIDES)
        if (LBC(isouth, isFsur, ng) % acquire):
            bry_pgr = -g * (zeta[idxReducedPhysicsBC, know] - boundary.zeta_south[i]) * 0.5 * grid.pn[idxReducedPhysicsBC]
        else:
            bry_pgr = -g * (zeta(idxReducedPhysicsBC, know) - zeta(idxReducedPhysicsBC2, kNow)) * 0.5 * RtoV(grid.pn,
                                                                                                             idxReducedPhysicsBC)

            # ifdef UV_COR
            bry_cor = UtoV(ubar, idxReducedPhysicsBC) * RtoV(grid.f, idxReducedPhysicsBC)
            # else
            bry_cor = 0.0
            # endif

            zhR = RtoV(grid.h, idx) + RtoV(grid.zeta[:, :, know], idx)
            izhR = 1.0 / zhR

            bry_str = izhR*(forces.svstr[idxReducedPhysicsBC] - forces.bvstr[idxReducedPhysicsBC])

            Ce=sqrt(g*izhR)

        if WEC_VF:
    #             here we reduce pgr from depth 2m to pgr is 0 at 0 depth.
                  bry_pgr=MIN(0.5 * (grid.h[idxFlatherBC] + zeta[idxFlatherBC, kNow]), 1.0)*bry_pgr
                  cff1 = 1.5*pi - forces.Dwave[idxFlatherBC] - grid.angler[idxFlatherBC]

                  waven = 2.0*pi/MAX(forces.Lwave[idxFlatherBC], Lwave_min)
                  waveny = waven*SIN(cff1)
                  sigma = sqrt(MAX(g*waven*TANH(waven*zhR), eps))
                  osigma = 1.0/sigma

                  cff1 = (1.0 - wec_alpha(ng))*osigma*forces.Dissip_break[idxFlatherBC]*waveny
                  bry_wec=cff1

        bry_wec = 0.0
        if WEC_ROLLER:
              bry_wec = bry_wec + osigma*forces.Dissip_roller[idxFlatherBC]*waveny

              cff2 = grid.on_v[idxFlatherBC]*Ce

              bry_val= vbar[idxFlatherBC, kNow] + cff2*(bry_pgr + bry_cor + bry_wec + bry_str)
    else:
        bry_val=boundary.vbar_south(i)

    if defined ATM_PRESS && defined PRESS_COMPENSATE:
          vbar[idxFlatherBC,kout] = bry_val - Ce*RtoV(zeta[:,know], idxFlatherBC) + +
                                            fac * 2*RtoV(forces.Pair, idxFlatherBC) -
             &                              2.0*OneAtm)) -
                                            boundary.zeta_south[i])
    else:
          vbar[idxFlatherBC,kout] = bry_val - Ce*(RtoV(zeta(idxFlatherBC,know) -  boundary.zeta_south(i))








    # Lateral boundary conditions at the southern edge.
    # -------------------------------------------------

    # Southern edge, implicit upstream radiation condition.

    grad[idxImplicitUpstreamRadiationBC]  = Vdiff(vbar[:, kNow], idxImplicitUpstreamRadiationBC)
    grad[idxImplicitUpstreamRadiationBC2] = Vdiff(vbar[:, kNow], idxImplicitUpstreamRadiationBC2)




          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              dVdt = vbar(i, Jstr + 1, kNow) - vbar(i, Jstr + 1, kout)
              dVde = vbar(i, Jstr + 1, kout) - vbar(i, Jstr + 2, kout)

              IF (LBC(isouth,isVbar,ng)%nudging) THEN
                IF (LnudgeM2CLM(ng)) THEN
                  obc_out=RtoV(CLIMA(ng)%M2nudgcof(i,Jstr-1))
                  obc_in =obcfac(ng)*obc_out
                ELSE
                  obc_out=M2obc_out(ng,isouth)
                  obc_in =M2obc_in (ng,isouth)
                END IF

                IF ((dVdt*dVde)<0.0) THEN
                  tau=obc_in
                ELSE
                  tau=obc_out
                END IF

#ifdef IMPLICIT_NUDGING
                IF (tau > 0.0) tau=1.0/tau
#else
                tau = tau*dt2d
#endif
              END IF

              IF ((dVdt*dVde) < 0.0) dVdt=0.0
              IF ((dVdt*(grad(i  ,Jstr+1) + grad(i+1,Jstr+1))) > 0.0) THEN
                dVdx=grad(i  ,Jstr+1)
              ELSE
                dVdx=grad(i+1,Jstr+1)
              END IF

              cff=MAX(dVdx*dVdx + dVde*dVde,eps)
#ifdef RADIATION_2D
              Cx = MIN(cff,MAX(dVdt*dVdx,-cff))
#else
              Cx = 0.0_r8
#endif
              Ce=dVdt*dVde

#if defined CELERITY_WRITE && defined FORWARD_WRITE
              BOUNDARY(ng)%vbar_south_Cx(i)=Cx
              BOUNDARY(ng)%vbar_south_Ce(i)=Ce
              BOUNDARY(ng)%vbar_south_C2(i)=cff
#endif
              vbar(i,Jstr,kout)=(cff * vbar(i, Jstr, kNow) +  Ce * vbar(i,Jstr+1,kout) - &
                                 & MAX(Cx,0.0_r8)*grad(i  ,Jstr) - MIN(Cx,0.0_r8)*grad(i+1,Jstr))/(cff+Ce)

              IF (LBC(isouth,isVbar,ng)%nudging) THEN
#ifdef IMPLICIT_NUDGING
                phi=dt(ng)/(tau+dt(ng))
                vbar(i,Jstr,kout)=(1.0_r8-phi)*vbar(i,Jstr,kout) + phi*BOUNDARY(ng)%vbar_south(i)

#else
                vbar(i,Jstr,kout)=vbar(i,Jstr,kout) + tau*(BOUNDARY(ng) % vbar_south(i) - vbar(i, Jstr, kNow))
#endif
              END IF

#ifdef MASKING
              vbar(i,Jstr,kout)=vbar(i,Jstr,kout)*GRID(ng)%vmask(i,Jstr)
#endif
            END IF
          END DO

# ------------------------------------------------
# ------------------------------------------------
# ------------------------------------------------
#   Southern edge, Shchepetkin boundary condition (Maison et al., 2010).

        ELSE IF (LBC(isouth,isVbar,ng)%Shchepetkin) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i)) THEN
#if defined SSH_TIDES && !defined UV_TIDES
              IF (LBC(isouth,isFsur,ng)%acquire) THEN
                bry_pgr=-g*(zeta(i, Jstr, kNow) - &
                            & BOUNDARY(ng) % zeta_south(i))*                &
     &                  0.5_r8*GRID(ng)%pn(i,Jstr)
              ELSE
                bry_pgr= -g * (zeta(i, Jstr, kNow) - &
                               & zeta(i, Jstr - 1, kNow)) * &
     &                  0.5_r8*(GRID(ng)%pn(i,Jstr-1)+                  &
     &                          GRID(ng)%pn(i,Jstr  ))
              END IF
# ifdef UV_COR
              bry_cor=-0.125_r8 * (ubar(i, Jstr - 1, kNow) + &
                                   & ubar(i + 1, Jstr - 1, kNow) + &
                                   & ubar(i, Jstr, kNow) + &
                                   & ubar(i + 1, Jstr, kNow)) * &
     &                          (GRID(ng)%f(i,Jstr-1)+                  &
     &                           GRID(ng)%f(i,Jstr  ))
# else
              bry_cor=0.0_r8
# endif
              cff1=1.0_r8 / (0.5_r8*(GRID(ng) % h(i,Jstr-1) + &
                      & zeta(i, Jstr - 1, kNow) + &
                      & GRID(ng) % h(i,Jstr  ) + &
                      & zeta(i, Jstr, kNow)))
              bry_str=cff1*(FORCES(ng)%svstr(i,Jstr)-                   &
     &                      FORCES(ng)%bvstr(i,Jstr))
              Ce=1.0_r8/SQRT(g*0.5_r8*(GRID(ng) % h(i,Jstr-1) + &
                                       & zeta(i, Jstr - 1, kNow) + &
                                       & GRID(ng) % h(i,Jstr  ) + &
                                       & zeta(i, Jstr, kNow)))
              cff2=GRID(ng)%on_v(i,Jstr)*Ce
!!            cff2=dt2d
              bry_val= vbar(i, Jstr + 1, kNow) + &
     &                cff2*(bry_pgr+                                    &
     &                      bry_cor+                                    &
     &                      bry_str)
#else
              bry_val=BOUNDARY(ng)%vbar_south(i)
#endif
#ifdef WET_DRY
              cff=0.5_r8*(GRID(ng) % h(i,Jstr-1) + &
                          & zeta(i, Jstr - 1, kNow) + &
                          & GRID(ng) % h(i,Jstr  ) + &
                          & zeta(i, Jstr, kNow))
#else
              cff=0.5_r8*(GRID(ng)%h(i,Jstr-1)+                         &
     &                    GRID(ng)%h(i,Jstr  ))
#endif
              cff1=SQRT(g/cff)
              Ce=dt2d*cff1*cff*0.5_r8*(GRID(ng)%pn(i,Jstr-1)+           &
     &                                 GRID(ng)%pn(i,Jstr  ))
              Ze= (0.5_r8+Ce) * zeta(i, Jstr, kNow) + &
     &           (0.5_r8-Ce)*zeta(i, Jstr - 1, kNow)
              IF (Ce.gt.Co) THEN
                cff2=(1.0_r8-Co/Ce)**2
                cff3=zeta(i,Jstr,kout)+                                 &
     & Ce * zeta(i, Jstr - 1, kNow) - &
     &               (1.0_r8+Ce)*zeta(i, Jstr, kNow)
                Ze=Ze+cff2*cff3
              END IF
              vbar(i,Jstr,kout)=0.5_r8*                                 &
     &                          ((1.0_r8-Ce) * vbar(i, Jstr, kNow) + &
                                 & Ce * vbar(i, Jstr + 1, kNow) + &
                                 & bry_val - &
                                 & cff1 * (Ze-BOUNDARY(ng)%zeta_south(i)))
#ifdef MASKING
              vbar(i,Jstr,kout)=vbar(i,Jstr,kout)*                      &
     &                          GRID(ng)%vmask(i,Jstr)
#endif
            END IF
          END DO


# ------------------------------------------------
# ------------------------------------------------
# ------------------------------------------------

#ifdef NESTING
# !
# !  If refinement grid and southern edge, impose mass flux from donor
# !  coaser grid for volume and mass conservation.
# !
        ELSE IF (LBC(isouth,isVbar,ng)%nested) THEN
          DO cr=1,Ncontact
            dg=Vcontact(cr)%donor_grid
            rg=Vcontact(cr)%receiver_grid
            IF (RefinedGrid(ng).and.                                    &
     &          (rg.eq.ng).and.(DXmax(dg).gt.DXmax(rg))) THEN
              told=3-RollingIndex(cr)
              tnew=RollingIndex(cr)
              DO i=Istr,Iend
                m=BRY_CONTACT(isouth,cr)%C2Bindex(i)
                Idg=Vcontact(cr)%Idg(m)             ! for debugging
                Jdg=Vcontact(cr)%Jdg(m)             ! purposes
                cff=0.5_r8*GRID(ng)%om_v(i,Jstr)*                       &
     &              (GRID(ng)%h(i,Jstr-1)+zeta(i,Jstr-1,kout)+          &
     &               GRID(ng)%h(i,Jstr  )+zeta(i,Jstr  ,kout))
                cff1=GRID(ng)%om_v(i,Jstr)/REFINED(cr)%om_v(m)
                bry_val=cff1*REFINED(cr)%DV_avg2(1,m,tnew)/cff
# ifdef WEC
                bry_val=bry_val-OCEAN(ng)%vbar_stokes(i,Jstr)
# endif
# ifdef MASKING
                bry_val=bry_val*GRID(ng)%vmask(i,Jstr)
# endif
# ifdef NESTING_DEBUG
                BRY_CONTACT(isouth,cr)%Mflux(i)=cff*bry_val
# endif
                vbar(i,Jstr,kout)=bry_val
              END DO
            END IF
          END DO
#endif
        END IF
      END IF

