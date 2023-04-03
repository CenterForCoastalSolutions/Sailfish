
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
if (LBC(isouth, isFsur, ng) % acquire):
    bry_pgr = -g * (zeta[idxReducedPhysicsBC, know] - boundary.zeta_south[i]) * 0.5 * grid.pn[idxReducedPhysicsBC]
else:
    bry_pgr = -g * (zeta(idxReducedPhysicsBC, know) - zeta(idxReducedPhysicsBC2, kNow)) * 0.5 * RtoV(grid.pn,
                                                                                                     idxReducedPhysicsBC) )

    # ifdef UV_COR
    # ubar interpolated at at V point.
    #        bry_cor=-0.125_r8 * (ubar(i,     Jstr - 1, kNow) +
    #                             ubar(i + 1, Jstr - 1, kNow) +
    #                             ubar(i,     Jstr,     kNow) +
    #                             ubar(i + 1, Jstr,     kNow)) *
    # &                          (grid.f(i,Jstr-1) + grid.f(i,Jstr  ))

    bry_cor = UtoV(ubar, idxReducedPhysicsBC) * RtoV(grid.f, idxReducedPhysicsBC)
    # else
    bry_cor = 0.0
    # endif

    zhR = RtoV(grid.h, idx) + RtoV(grid.zeta[:, :, know], idx)
    # zhR = 0.5 * (grid.h(i, Jstr - 1) + zeta(i, Jstr - 1, kNow) +
    #              grid.h(i, Jstr) + zeta(i, Jstr, kNow))
    bry_str = (forces.svstr(i, Jstr) - forces.bvstr(i, Jstr)) / zhR
    vbar(i, Jstr, kout) = vbar(i, Jstr, kNow) + dt2d * (bry_pgr + bry_cor + bry_str)

    END
    IF
    END
    DO

# C++

# Computes the gradient at the boundary

# NO IFDEF HERE???
if (LBC(isouth, isFsur, ng) % acquire):
    bry_pgr = -g*(zeta_Know - boundary.zeta_south[i])*0.5*pn
else:
    bry_pgr = -g * Diff(zeta_Know - zeta(idxReducedPhysicsBC2, kNow)) * 0.5 * RtoV(grid.pn,
                                                                                                     idxReducedPhysicsBC) )

    # ifdef UV_COR
    # ubar interpolated at at V point.
    #        bry_cor=-0.125_r8 * (ubar(i,     Jstr - 1, kNow) +
    #                             ubar(i + 1, Jstr - 1, kNow) +
    #                             ubar(i,     Jstr,     kNow) +
    #                             ubar(i + 1, Jstr,     kNow)) *
    # &                          (grid.f(i,Jstr-1) + grid.f(i,Jstr  ))

    bry_cor = UtoV(ubar, idxReducedPhysicsBC) * RtoV(grid.f, idxReducedPhysicsBC)
    # else
    bry_cor = 0.0
    # endif

    zhR = RtoV(grid.h, idx) + RtoV(grid.zeta[:, :, know], idx)
    # zhR = 0.5 * (grid.h(i, Jstr - 1) + zeta(i, Jstr - 1, kNow) +
    #              grid.h(i, Jstr) + zeta(i, Jstr, kNow))
    bry_str = (forces.svstr(i, Jstr) - forces.bvstr(i, Jstr)) / zhR
    vbar(i, Jstr, kout) = vbar(i, Jstr, kNow) + dt2d * (bry_pgr + bry_cor + bry_str)

    END
    IF
    END
    DO



########### FLATHER
########### FLATHER
########### FLATHER


        # elif (LBC(isouth,isVbar,ng)%Flather):
        #   DO i=Istr,Iend
        #     IF (LBC_apply(ng)%south(i)) THEN


#if (defined SSH_TIDES && !defined UV_TIDES)
              # Computes the gradient at the boundary
              if (LBC(isouth,isFsur,ng)%acquire):
                bry_pgr= -g*(zeta(i, Jstr, kNow) -  boundary.zeta_south(i))*0.5*grid.pn(i,Jstr)
              else:
                bry_pgr= -g*(zeta(i, Jstr, kNow) - zeta(i, Jstr - 1, kNow))*0.5*(grid.pn(i,Jstr-1) + grid.pn(i,Jstr))

# ifdef UV_COR
              bry_cor = UtoV(ubar, idxFlatherBC)*RtoV(grid.f, idxFlatherBC)

# else
              bry_cor = 0.0
# endif

              zhR = RtoV(grid.h[idxFlatherBC] + zeta[idxFlatherBC])
              izhR = 1.0/zhR
              bry_str=izhR*(forces.svstr(i,Jstr) - forces.bvstr(i,Jstr))

              Ce=sqrt(g*izhR)
# ifdef WEC_VF
#             here we reduce pgr from depth 2m to pgr is 0 at 0 depth.
              bry_pgr=MIN(0.5 * (grid.h[idxFlatherBC] + zeta[idxFlatherBC, kNow]), 1.0)*bry_pgr
              cff = 1.0 / zhR
              cff1 = 1.5*pi - forces.Dwave[idxFlatherBC] - grid.angler[idxFlatherBC]
              waven = 2.0*pi/MAX(forces.Lwave[idxFlatherBC], Lwave_min)
              waveny = waven*SIN(cff1)
              sigma = sqrt(MAX(g*waven*TANH(waven*zhR), eps))
              osigma = 1.0/sigma

              cff1 = (1.0 - wec_alpha(ng))*osigma*forces.Dissip_break[idxFlatherBC]*waveny
              bry_wec=cff1
#  ifdef WEC_ROLLER
              bry_wec = bry_wec + osigma*forces.Dissip_roller[idxFlatherBC]*waveny
#  endif
# else
              bry_wec=0.0_r8
# endif
              cff2=grid.on_v(i,Jstr)*Ce

              bry_val= vbar[idxFlatherBC, kNow] + cff2*(bry_pgr + bry_cor + bry_wec + bry_str)
#else
              bry_val=boundary.vbar_south(i)
#endif
#if defined ATM_PRESS && defined PRESS_COMPENSATE
              vbar[i,Jstr,kout] = bry_val - Ce*(0.5*(zeta(i,Jstr-1,know) +  zeta(i, Jstr, kNow) +
                                                fac * (forces.Pair(i,Jstr-1) + forces.Pair(i,Jstr  ) -
                 &                              2.0*OneAtm)) -
                                                boundary.zeta_south(i))
#else
              vbar[idxFlatherBC,kout] = bry_val - Ce*(RtoV(zeta(idxFlatherBC,know) -  boundary.zeta_south(i))
#endif

            END IF
          END DO







    # Lateral boundary conditions at the southern edge.
    # -------------------------------------------------

    if domain.Southern_Edge:

        # ------------------------------------------------
        # ------------------------------------------------
        # ------------------------------------------------

        # Southern edge, implicit upstream radiation condition.
        if LBC[isouth,isVbar,ng].radiation:
            grad[Istr:Iend+1, Jstr  ] = vbar[Istr:Iend+1, Jstr,   kNow] - vbar[Istr - 1:Iend, Jstr    , kNow]
            grad[Istr:Iend+1, Jstr+1] = vbar[Istr:Iend+1, Jstr+1, kNow] - vbar[Istr - 1:Iend, Jstr + 1, kNow]

     #      DO i=Istr,Iend+1
     #        grad(i,Jstr  )= vbar(i, Jstr, kNow) - &
     # &                     vbar(i - 1, Jstr, kNow)
     #        grad(i,Jstr+1)= vbar(i, Jstr + 1, kNow) - &
     # &                     vbar(i - 1, Jstr + 1, kNow)
     #      END DO


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
#   Southern edge, Flather boundary condition.

        ELSE IF (LBC(isouth,isVbar,ng)%Flather) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i)) THEN
!#if (defined SSH_TIDES && !defined UV_TIDES) || defined INWAVE_MODEL
#if (defined SSH_TIDES && !defined UV_TIDES)
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
# ifdef WEC_VF
!             here we reduce pgr from depth 2m to pgr is 0 at 0 depth.
              bry_pgr=MIN(0.5_r8 * (GRID(ng) % h(i,Jstr) + &
                                    & zeta(i, Jstr, kNow)), 1.0_r8)*bry_pgr
              cff=1.0_r8 / (0.5_r8*(GRID(ng) % h(i,Jstr-1) + &
                     & zeta(i, Jstr - 1, kNow) + &
                     & GRID(ng) % h(i,Jstr  ) + &
                     & zeta(i, Jstr, kNow)))
              cff1=1.5_r8*pi-FORCES(ng)%Dwave(i,Jstr)-                  &
     &             GRID(ng)%angler(i,Jstr)
              waven=2.0_r8*pi/MAX(FORCES(ng)%Lwave(i,Jstr),Lwave_min)
              waveny=waven*SIN(cff1)
              sigma=SQRT(MAX(g*waven*TANH(waven/cff),eps))
              osigma=1.0_r8/sigma
              cff1=(1.0_r8-wec_alpha(ng))*osigma*                       &
     &             FORCES(ng)%Dissip_break(i,Jstr)*waveny
              bry_wec=cff1
#  ifdef WEC_ROLLER
              bry_wec=bry_wec+osigma*FORCES(ng)%Dissip_roller(i,Jstr)*  &
     &                        waveny
#  endif
!             bry_wec=0.0_r8   !for now
# else
              bry_wec=0.0_r8
# endif
              cff2=GRID(ng)%on_v(i,Jstr)*Ce
!!            cff2=dt2d
              bry_val= vbar(i, Jstr + 1, kNow) + &
     &                cff2*(bry_pgr+                                    &
     &                      bry_cor+                                    &
     &                      bry_wec+                                    &
     &                      bry_str)
#else
              bry_val=BOUNDARY(ng)%vbar_south(i)
#endif
              cff=1.0_r8 / (0.5_r8*(GRID(ng) % h(i,Jstr-1) + &
                     & zeta(i, Jstr - 1, kNow) + &
                     & GRID(ng) % h(i,Jstr  ) + &
                     & zeta(i, Jstr, kNow)))
              Ce=SQRT(g*cff)
#if defined ATM_PRESS && defined PRESS_COMPENSATE
              vbar(i,Jstr,kout)=bry_val-                                &
     & Ce * (0.5_r8*                             &
     &                              (zeta(i,Jstr-1,know) + &
       & zeta(i, Jstr, kNow) + &
       & fac * (FORCES(ng)%Pair(i,Jstr-1)+    &
     &                                    FORCES(ng)%Pair(i,Jstr  )-    &
     &                                    2.0_r8*OneAtm))-              &
     &                              BOUNDARY(ng)%zeta_south(i))
#else
              vbar(i,Jstr,kout)=bry_val-                                &
     &                          Ce*(0.5_r8*(zeta(i,Jstr-1,know)+        &
     &                                      zeta(i, Jstr, kNow))- &
                                                                & BOUNDARY(ng) % zeta_south(i))
#endif
#ifdef MASKING
              vbar(i,Jstr,kout)=vbar(i,Jstr,kout)*                      &
     &                          GRID(ng)%vmask(i,Jstr)
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
#   Southern edge, reduced-physics boundary condition.

        ELSE IF (LBC(isouth,isVbar,ng)%reduced) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i)) THEN
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
#ifdef UV_COR
              bry_cor=-0.125_r8 * (ubar(i, Jstr - 1, kNow) + &
                                   & ubar(i + 1, Jstr - 1, kNow) + &
                                   & ubar(i, Jstr, kNow) + &
                                   & ubar(i + 1, Jstr, kNow)) * &
     &                          (GRID(ng)%f(i,Jstr-1)+                  &
     &                           GRID(ng)%f(i,Jstr  ))
#else
              bry_cor=0.0_r8
#endif
              cff=1.0_r8 / (0.5_r8*(GRID(ng) % h(i,Jstr-1) + &
                     & zeta(i, Jstr - 1, kNow) + &
                     & GRID(ng) % h(i,Jstr  ) + &
                     & zeta(i, Jstr, kNow)))
              bry_str=cff*(FORCES(ng)%svstr(i,Jstr)-                    &
     &                     FORCES(ng)%bvstr(i,Jstr))
              vbar(i,Jstr,kout)= vbar(i, Jstr, kNow) + &
     &                          dt2d*(bry_pgr+                          &
     &                                bry_cor+                          &
     &                                bry_str)
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


!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the northern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
!
!  Northern edge, implicit upstream radiation condition.
!
        IF (LBC(inorth,isVbar,ng)%radiation) THEN
          DO i=Istr,Iend+1
            grad(i,Jend  )= vbar(i, Jend, kNow) - &
     &                     vbar(i - 1, Jend, kNow)
            grad(i,Jend+1)= vbar(i, Jend + 1, kNow) - &
     &                     vbar(i - 1, Jend + 1, kNow)
          END DO
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              dVdt= vbar(i, Jend, kNow) - vbar(i, Jend, kout)
              dVde=vbar(i,Jend,kout)-vbar(i,Jend-1,kout)

              IF (LBC(inorth,isVbar,ng)%nudging) THEN
                IF (LnudgeM2CLM(ng)) THEN
                  obc_out=0.5_r8*                                       &
     &                    (CLIMA(ng)%M2nudgcof(i,Jend  )+               &
     &                     CLIMA(ng)%M2nudgcof(i,Jend+1))
                  obc_in =obcfac(ng)*obc_out
                ELSE
                  obc_out=M2obc_out(ng,inorth)
                  obc_in =M2obc_in (ng,inorth)
                END IF
                IF ((dVdt*dVde).lt.0.0_r8) THEN
                  tau=obc_in
                ELSE
                  tau=obc_out
                END IF
#ifdef IMPLICIT_NUDGING
                IF (tau.gt.0.0_r8) tau=1.0_r8/tau
#else
                tau=tau*dt2d
#endif
              END IF

              IF ((dVdt*dVde).lt.0.0_r8) dVdt=0.0_r8
              IF ((dVdt*(grad(i  ,Jend)+                                &
     &                   grad(i+1,Jend))).gt.0.0_r8) THEN
                dVdx=grad(i  ,Jend)
              ELSE
                dVdx=grad(i+1,Jend)
              END IF
              cff=MAX(dVdx*dVdx+dVde*dVde,eps)
#ifdef RADIATION_2D
              Cx=MIN(cff,MAX(dVdt*dVdx,-cff))
#else
              Cx=0.0_r8
#endif
              Ce=dVdt*dVde
#if defined CELERITY_WRITE && defined FORWARD_WRITE
              BOUNDARY(ng)%vbar_north_Cx(i)=Cx
              BOUNDARY(ng)%vbar_north_Ce(i)=Ce
              BOUNDARY(ng)%vbar_north_C2(i)=cff
#endif
              vbar(i,Jend+1,kout)=(cff * vbar(i, Jend + 1, kNow) + &
                                   & Ce * vbar(i,Jend  ,kout) - &
                                   & MAX(Cx,0.0_r8)*grad(i  ,Jend+1)-     &
     &                             MIN(Cx,0.0_r8)*grad(i+1,Jend+1))/    &
     &                            (cff+Ce)

              IF (LBC(inorth,isVbar,ng)%nudging) THEN
#ifdef IMPLICIT_NUDGING
                phi=dt(ng)/(tau+dt(ng))
                vbar(i,Jend+1,kout)=(1.0_r8-phi)*vbar(i,Jend+1,kout)+   &
     &                              phi*BOUNDARY(ng)%vbar_north(i)

#else
                vbar(i,Jend+1,kout)=vbar(i,Jend+1,kout)+                &
     &                              tau*(BOUNDARY(ng) % vbar_north(i) - &
                                         & vbar(i, Jend + 1, kNow))
#endif
              END IF
#ifdef MASKING
              vbar(i,Jend+1,kout)=vbar(i,Jend+1,kout)*                  &
     &                            GRID(ng)%vmask(i,Jend+1)
#endif
            END IF
          END DO
!
!  Northern edge, Flather boundary condition.
!
        ELSE IF (LBC(inorth,isVbar,ng)%Flather) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
!#if (defined SSH_TIDES && !defined UV_TIDES) || defined INWAVE_MODEL
#if (defined SSH_TIDES && !defined UV_TIDES)
              IF (LBC(inorth,isFsur,ng)%acquire) THEN
                bry_pgr= -g * (BOUNDARY(ng) % zeta_north(i) - &
                               & zeta(i, Jend, kNow)) * &
     &                  0.5_r8*GRID(ng)%pn(i,Jend)
              ELSE
                bry_pgr= -g * (zeta(i, Jend + 1, kNow) - &
                               & zeta(i, Jend, kNow)) * &
     &                  0.5_r8*(GRID(ng)%pn(i,Jend  )+                  &
     &                          GRID(ng)%pn(i,Jend+1))
              END IF
# ifdef UV_COR
              bry_cor=-0.125_r8 * (ubar(i, Jend, kNow) + &
                                   & ubar(i + 1, Jend, kNow) + &
                                   & ubar(i, Jend + 1, kNow) + &
                                   & ubar(i + 1, Jend + 1, kNow)) * &
     &                          (GRID(ng)%f(i,Jend  )+                  &
     &                           GRID(ng)%f(i,Jend+1))
# else
              bry_cor=0.0_r8
# endif
              cff1=1.0_r8 / (0.5_r8*(GRID(ng) % h(i,Jend  ) + &
                      & zeta(i, Jend, kNow) + &
                      & GRID(ng) % h(i,Jend+1) + &
                      & zeta(i, Jend + 1, kNow)))
              bry_str=cff1*(FORCES(ng)%svstr(i,Jend+1)-                 &
     &                      FORCES(ng)%bvstr(i,Jend+1))
              Ce=1.0_r8/SQRT(g*0.5_r8*(GRID(ng) % h(i,Jend+1) + &
                                       & zeta(i, Jend + 1, kNow) + &
                                       & GRID(ng) % h(i,Jend  ) + &
                                       & zeta(i, Jend, kNow)))
# ifdef WEC_VF
!             here we reduce pgr from depth 2m to pgr is 0 at 0 depth.
              bry_pgr=MIN(0.5_r8 * (GRID(ng) % h(i,Jend) + &
                                    & zeta(i, Jend, kNow)), 1.0_r8)*bry_pgr
              cff=1.0_r8 / (0.5_r8*(GRID(ng) % h(i,Jend  ) + &
                     & zeta(i, Jend, kNow) + &
                     & GRID(ng) % h(i,Jend+1) + &
                     & zeta(i, Jend + 1, kNow)))
              cff1=1.5_r8*pi-FORCES(ng)%Dwave(i,Jend)-                  &
     &             GRID(ng)%angler(i,Jend)
              waven=2.0_r8*pi/MAX(FORCES(ng)%Lwave(i,Jend),Lwave_min)
              waveny=waven*SIN(cff1)
              sigma=SQRT(MAX(g*waven*TANH(waven/cff),eps))
              osigma=1.0_r8/sigma
              cff1=(1.0_r8-wec_alpha(ng))*osigma*                       &
     &             FORCES(ng)%Dissip_break(i,Jend)*waveny
              bry_wec=cff1
#  ifdef WEC_ROLLER
              bry_wec=bry_wec+osigma*FORCES(ng)%Dissip_roller(i,Jend)*  &
     &                        waveny
#  endif
!             bry_wec=0.0_r8   ! for now
# else
              bry_wec=0.0_r8
# endif
              cff2=GRID(ng)%on_v(i,Jend+1)*Ce
!!            cff2=dt2d
              bry_val= vbar(i, Jend, kNow) + &
     &                cff2*(bry_pgr+                                    &
     &                      bry_cor+                                    &
     &                      bry_wec+                                    &
     &                      bry_str)
#else
              bry_val=BOUNDARY(ng)%vbar_north(i)
#endif
              cff=1.0_r8 / (0.5_r8*(GRID(ng) % h(i,Jend  ) + &
                     & zeta(i, Jend, kNow) + &
                     & GRID(ng) % h(i,Jend+1) + &
                     & zeta(i, Jend + 1, kNow)))
              Ce=SQRT(g*cff)
#if defined ATM_PRESS && defined PRESS_COMPENSATE
            vbar(i,Jend+1,kout)=bry_val+                                &
     & Ce * (0.5_r8*                             &
     &                              (zeta(i,Jend  ,know) + &
       & zeta(i, Jend + 1, kNow) + &
       & fac * (FORCES(ng)%Pair(i,Jend  )+    &
     &                                    FORCES(ng)%Pair(i,Jend+1)-    &
     &                                    2.0_r8*OneAtm))-              &
     &                              BOUNDARY(ng)%zeta_north(i))
#else
              vbar(i,Jend+1,kout)=bry_val+                              &
     &                            Ce*(0.5_r8*(zeta(i,Jend  ,know)+      &
     &                                        zeta(i, Jend + 1, kNow))- &
                                                                      & BOUNDARY(ng) % zeta_north(i))
#endif
#ifdef MASKING
              vbar(i,Jend+1,kout)=vbar(i,Jend+1,kout)*                  &
     &                            GRID(ng)%vmask(i,Jend+1)
#endif
            END IF
          END DO
!
!  Northern edge, Shchepetkin boundary condition (Maison et al., 2010).
!
        ELSE IF (LBC(inorth,isVbar,ng)%Shchepetkin) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
#if defined SSH_TIDES && !defined UV_TIDES
              IF (LBC(inorth,isFsur,ng)%acquire) THEN
                bry_pgr= -g * (BOUNDARY(ng) % zeta_north(i) - &
                               & zeta(i, Jend, kNow)) * &
     &                  0.5_r8*GRID(ng)%pn(i,Jend)
              ELSE
                bry_pgr= -g * (zeta(i, Jend + 1, kNow) - &
                               & zeta(i, Jend, kNow)) * &
     &                  0.5_r8*(GRID(ng)%pn(i,Jend  )+                  &
     &                          GRID(ng)%pn(i,Jend+1))
              END IF
# ifdef UV_COR
              bry_cor=-0.125_r8 * (ubar(i, Jend, kNow) + &
                                   & ubar(i + 1, Jend, kNow) + &
                                   & ubar(i, Jend + 1, kNow) + &
                                   & ubar(i + 1, Jend + 1, kNow)) * &
     &                          (GRID(ng)%f(i,Jend  )+                  &
     &                           GRID(ng)%f(i,Jend+1))
# else
              bry_cor=0.0_r8
# endif
              cff1=1.0_r8 / (0.5_r8*(GRID(ng) % h(i,Jend  ) + &
                      & zeta(i, Jend, kNow) + &
                      & GRID(ng) % h(i,Jend+1) + &
                      & zeta(i, Jend + 1, kNow)))
              bry_str=cff1*(FORCES(ng)%svstr(i,Jend+1)-                 &
     &                      FORCES(ng)%bvstr(i,Jend+1))
              Ce=1.0_r8/SQRT(g*0.5_r8*(GRID(ng) % h(i,Jend+1) + &
                                       & zeta(i, Jend + 1, kNow) + &
                                       & GRID(ng) % h(i,Jend  ) + &
                                       & zeta(i, Jend, kNow)))
              cff2=GRID(ng)%on_v(i,Jend+1)*Ce
!!            cff2=dt2d
              bry_val= vbar(i, Jend, kNow) + &
     &                cff2*(bry_pgr+                                    &
     &                      bry_cor+                                    &
     &                      bry_str)
#else
              bry_val=BOUNDARY(ng)%vbar_north(i)
#endif
#ifdef WET_DRY
              cff=0.5_r8*(GRID(ng) % h(i,Jend  ) + &
                          & zeta(i, Jend, kNow) + &
                          & GRID(ng) % h(i,Jend+1) + &
                          & zeta(i, Jend + 1, kNow))
#else
              cff=0.5_r8*(GRID(ng)%h(i,Jend  )+                         &
     &                    GRID(ng)%h(i,Jend+1))
#endif
              cff1=SQRT(g/cff)
              Ce=dt2d*cff1*cff*0.5_r8*(GRID(ng)%pn(i,Jend  )+           &
     &                                 GRID(ng)%pn(i,Jend+1))
              Ze= (0.5_r8+Ce) * zeta(i, Jend, kNow) + &
     &           (0.5_r8-Ce)*zeta(i, Jend + 1, kNow)
              IF (Ce.gt.Co) THEN
                cff2=(1.0_r8-Co/Ce)**2
                cff3=zeta(i,Jend,kout)+                                 &
     & Ce * zeta(i, Jend + 1, kNow) - &
     &               (1.0_r8+Ce)*zeta(i, Jend, kNow)
                Ze=Ze+cff2*cff3
              END IF
              vbar(i,Jend+1,kout)=0.5_r8*                               &
     &                            ((1.0_r8-Ce) * vbar(i, Jend + 1, kNow) + &
                                   & Ce * vbar(i, Jend, kNow) + &
                                   & bry_val + &
                                   & cff1 * (Ze-BOUNDARY(ng)%zeta_north(i)))
#ifdef MASKING
              vbar(i,Jend+1,kout)=vbar(i,Jend+1,kout)*                  &
     &                            GRID(ng)%vmask(i,Jend+1)
#endif
            END IF
          END DO
!
!  Northern edge, clamped boundary condition.
!
        ELSE IF (LBC(inorth,isVbar,ng)%clamped) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              vbar(i,Jend+1,kout)=BOUNDARY(ng)%vbar_north(i)
#ifdef MASKING
              vbar(i,Jend+1,kout)=vbar(i,Jend+1,kout)*                  &
     &                            GRID(ng)%vmask(i,Jend+1)
#endif
            END IF
          END DO
!
!  Northern edge, gradient boundary condition.
!
        ELSE IF (LBC(inorth,isVbar,ng)%gradient) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              vbar(i,Jend+1,kout)=vbar(i,Jend,kout)
#ifdef MASKING
              vbar(i,Jend+1,kout)=vbar(i,Jend+1,kout)*                  &
     &                            GRID(ng)%vmask(i,Jend+1)
#endif
            END IF
          END DO
!
!  Northern edge, reduced-physics boundary condition.
!
        ELSE IF (LBC(inorth,isVbar,ng)%reduced) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              IF (LBC(inorth,isFsur,ng)%acquire) THEN
                bry_pgr= -g * (BOUNDARY(ng) % zeta_north(i) - &
                               & zeta(i, Jend, kNow)) * &
     &                  0.5_r8*GRID(ng)%pn(i,Jend)
              ELSE
                bry_pgr= -g * (zeta(i, Jend + 1, kNow) - &
                               & zeta(i, Jend, kNow)) * &
     &                  0.5_r8*(GRID(ng)%pn(i,Jend  )+                  &
     &                          GRID(ng)%pn(i,Jend+1))
              END IF
#ifdef UV_COR
              bry_cor=-0.125_r8 * (ubar(i, Jend, kNow) + &
                                   & ubar(i + 1, Jend, kNow) + &
                                   & ubar(i, Jend + 1, kNow) + &
                                   & ubar(i + 1, Jend + 1, kNow)) * &
     &                          (GRID(ng)%f(i,Jend  )+                  &
     &                           GRID(ng)%f(i,Jend+1))
#else
              bry_cor=0.0_r8
#endif
              cff=1.0_r8 / (0.5_r8*(GRID(ng) % h(i,Jend  ) + &
                     & zeta(i, Jend, kNow) + &
                     & GRID(ng) % h(i,Jend+1) + &
                     & zeta(i, Jend + 1, kNow)))
              bry_str=cff*(FORCES(ng)%svstr(i,Jend+1)-                  &
     &                     FORCES(ng)%bvstr(i,Jend+1))
              vbar(i,Jend+1,kout)= vbar(i, Jend + 1, kNow) + &
     &                            dt2d*(bry_pgr+                        &
     &                                  bry_cor+                        &
     &                                  bry_str)
#ifdef MASKING
              vbar(i,Jend+1,kout)=vbar(i,Jend+1,kout)*                  &
     &                            GRID(ng)%vmask(i,Jend+1)
#endif
            END IF
          END DO
!
!  Northern edge, closed boundary condition.
!
        ELSE IF (LBC(inorth,isVbar,ng)%closed) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              vbar(i,Jend+1,kout)=0.0_r8
            END IF
          END DO

#ifdef NESTING
!
!  If refinement grid and northern edge, impose mass flux from donor
!  coaser grid for volume and mass conservation.
!
        ELSE IF (LBC(inorth,isVbar,ng)%nested) THEN
          DO cr=1,Ncontact
            dg=Vcontact(cr)%donor_grid
            rg=Vcontact(cr)%receiver_grid
            IF (RefinedGrid(ng).and.                                    &
     &          (rg.eq.ng).and.(DXmax(dg).gt.DXmax(rg))) THEN
              told=3-RollingIndex(cr)
              tnew=RollingIndex(cr)
              DO i=Istr,Iend
                m=BRY_CONTACT(inorth,cr)%C2Bindex(i)
                Idg=Vcontact(cr)%Idg(m)             ! for debugging
                Jdg=Vcontact(cr)%Jdg(m)             ! purposes
                cff=0.5_r8*GRID(ng)%om_v(i,Jend+1)*                     &
     &              (GRID(ng)%h(i,Jend+1)+zeta(i,Jend+1,kout)+          &
     &               GRID(ng)%h(i,Jend  )+zeta(i,Jend  ,kout))
                cff1=GRID(ng)%om_v(i,Jend+1)/REFINED(cr)%om_v(m)
                bry_val=cff1*REFINED(cr)%DV_avg2(1,m,tnew)/cff
# ifdef WEC
                bry_val=bry_val-OCEAN(ng)%vbar_stokes(i,Jend+1)
# endif
# ifdef MASKING
                bry_val=bry_val*GRID(ng)%vmask(i,Jend+1)
# endif
# ifdef NESTING_DEBUG
                BRY_CONTACT(inorth,cr)%Mflux(i)=cff*bry_val
# endif
                vbar(i,Jend+1,kout)=bry_val
              END DO
            END IF
          END DO
#endif
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
!
!  Western edge, implicit upstream radiation condition.
!
        IF (LBC(iwest,isVbar,ng)%radiation) THEN
          DO j=JstrV-1,Jend
            grad(Istr-1,j)= vbar(Istr - 1, j + 1, kNow) - &
     &                     vbar(Istr - 1, j, kNow)
            grad(Istr  ,j)= vbar(Istr, j + 1, kNow) - &
     &                     vbar(Istr, j, kNow)
          END DO
          DO j=JstrV,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              dVdt= vbar(Istr, j, kNow) - vbar(Istr, j, kout)
              dVdx=vbar(Istr,j,kout)-vbar(Istr+1,j,kout)

              IF (LBC(iwest,isVbar,ng)%nudging) THEN
                IF (LnudgeM2CLM(ng)) THEN
                  obc_out=0.5_r8*                                       &
     &                    (CLIMA(ng)%M2nudgcof(Istr-1,j-1)+             &
     &                     CLIMA(ng)%M2nudgcof(Istr-1,j  ))
                  obc_in =obcfac(ng)*obc_out
                ELSE
                  obc_out=M2obc_out(ng,iwest)
                  obc_in =M2obc_in (ng,iwest)
                END IF
                IF ((dVdt*dVdx).lt.0.0_r8) THEN
                  tau=obc_in
                ELSE
                  tau=obc_out
                END IF
#ifdef IMPLICIT_NUDGING
                IF (tau.gt.0.0_r8) tau=1.0_r8/tau
#else
                tau=tau*dt2d
#endif
              END IF

              IF ((dVdt*dVdx).lt.0.0_r8) dVdt=0.0_r8
              IF ((dVdt*(grad(Istr,j-1)+                                &
     &                   grad(Istr,j  ))).gt.0.0_r8) THEN
                dVde=grad(Istr,j-1)
              ELSE
                dVde=grad(Istr,j  )
              END IF

              cff=MAX(dVdx*dVdx+dVde*dVde,eps)
              Cx=dVdt*dVdx
#ifdef RADIATION_2D
              Ce=MIN(cff,MAX(dVdt*dVde,-cff))
#else
              Ce=0.0_r8
#endif
#if defined CELERITY_WRITE && defined FORWARD_WRITE
              BOUNDARY(ng)%vbar_west_Cx(j)=Cx
              BOUNDARY(ng)%vbar_west_Ce(j)=Ce
              BOUNDARY(ng)%vbar_west_C2(j)=cff
#endif
              vbar(Istr-1,j,kout)=(cff * vbar(Istr - 1, j, kNow) + &
                                   & Cx * vbar(Istr  ,j,kout) - &
                                   & MAX(Ce,0.0_r8)*grad(Istr-1,j-1)-     &
     &                             MIN(Ce,0.0_r8)*grad(Istr-1,j  ))/    &
     &                            (cff+Cx)

              IF (LBC(iwest,isVbar,ng)%nudging) THEN
#ifdef IMPLICIT_NUDGING
                phi=dt(ng)/(tau+dt(ng))
                vbar(Istr-1,j,kout)=(1.0_r8-phi)*vbar(Istr-1,j,kout)+   &
     &                              phi*BOUNDARY(ng)%vbar_west(j)
#else
                vbar(Istr-1,j,kout)=vbar(Istr-1,j,kout)+                &
     &                              tau*(BOUNDARY(ng) % vbar_west(j) - &
                                         & vbar(Istr - 1, j, kNow))
#endif
              END IF
#ifdef MASKING
              vbar(Istr-1,j,kout)=vbar(Istr-1,j,kout)*                  &
     &                            GRID(ng)%vmask(Istr-1,j)
#endif
            END IF
          END DO
!
!  Western edge, Chapman boundary condition.
!
        ELSE IF (LBC(iwest,isVbar,ng)%Flather.or.                       &
     &           LBC(iwest,isVbar,ng)%reduced.or.                       &
     &           LBC(iwest,isVbar,ng)%Shchepetkin) THEN
          DO j=JstrV,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              cff=dt2d*0.5_r8*(GRID(ng)%pm(Istr,j-1)+                   &
     &                         GRID(ng)%pm(Istr,j  ))
              cff1=SQRT(g*0.5_r8*(GRID(ng) % h(Istr,j-1) + &
                                  & zeta(Istr, j - 1, kNow) + &
                                  & GRID(ng) % h(Istr,j  ) + &
                                  & zeta(Istr, j, kNow)))
              Cx=cff*cff1
              cff2=1.0_r8/(1.0_r8+Cx)
              vbar(Istr-1,j,kout)=cff2*(vbar(Istr - 1, j, kNow) + &
                                        & Cx * vbar(Istr,j,kout))
#ifdef MASKING
              vbar(Istr-1,j,kout)=vbar(Istr-1,j,kout)*                  &
     &                            GRID(ng)%vmask(Istr-1,j)
#endif
            END IF
          END DO
!
!  Western edge, clamped boundary condition.
!
        ELSE IF (LBC(iwest,isVbar,ng)%clamped) THEN
          DO j=JstrV,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              vbar(Istr-1,j,kout)=BOUNDARY(ng)%vbar_west(j)
#ifdef MASKING
              vbar(Istr-1,j,kout)=vbar(Istr-1,j,kout)*                  &
     &                            GRID(ng)%vmask(Istr-1,j)
#endif
            END IF
          END DO
!
!  Western edge, gradient boundary condition.
!
        ELSE IF (LBC(iwest,isVbar,ng)%gradient) THEN
          DO j=JstrV,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              vbar(Istr-1,j,kout)=vbar(Istr,j,kout)
#ifdef MASKING
              vbar(Istr-1,j,kout)=vbar(Istr-1,j,kout)*                  &
     &                            GRID(ng)%vmask(Istr-1,j)
#endif
            END IF
          END DO
!
!  Western edge, closed boundary condition: free slip (gamma2=1)  or
!                                           no   slip (gamma2=-1).
!
        ELSE IF (LBC(iwest,isVbar,ng)%closed) THEN
          IF (NSperiodic(ng)) THEN
            Jmin=JstrV
            Jmax=Jend
          ELSE
            Jmin=Jstr
            Jmax=JendR
          END IF
          DO j=Jmin,Jmax
            IF (LBC_apply(ng)%west(j)) THEN
              vbar(Istr-1,j,kout)=gamma2(ng)*vbar(Istr,j,kout)
#ifdef MASKING
              vbar(Istr-1,j,kout)=vbar(Istr-1,j,kout)*                  &
     &                            GRID(ng)%vmask(Istr-1,j)
#endif
            END IF
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the eastern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
!
!  Eastern edge, implicit upstream radiation condition.
!
        IF (LBC(ieast,isVbar,ng)%radiation) THEN
          DO j=JstrV-1,Jend
            grad(Iend  ,j)= vbar(Iend, j + 1, kNow) - &
     &                     vbar(Iend, j, kNow)
            grad(Iend+1,j)= vbar(Iend + 1, j + 1, kNow) - &
     &                     vbar(Iend + 1, j, kNow)
          END DO
          DO j=JstrV,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              dVdt= vbar(Iend, j, kNow) - vbar(Iend, j, kout)
              dVdx=vbar(Iend,j,kout)-vbar(Iend-1,j,kout)

              IF (LBC(ieast,isVbar,ng)%nudging) THEN
                IF (LnudgeM2CLM(ng)) THEN
                  obc_out=0.5_r8*                                       &
     &                    (CLIMA(ng)%M2nudgcof(Iend+1,j-1)+             &
     &                     CLIMA(ng)%M2nudgcof(Iend+1,j  ))
                  obc_in =obcfac(ng)*obc_out
                ELSE
                  obc_out=M2obc_out(ng,ieast)
                  obc_in =M2obc_in (ng,ieast)
                END IF
                IF ((dVdt*dVdx).lt.0.0_r8) THEN
                  tau=obc_in
                ELSE
                  tau=obc_out
                END IF
                tau=tau*dt2d
              END IF

              IF ((dVdt*dVdx).lt.0.0_r8) dVdt=0.0_r8
              IF ((dVdt*(grad(Iend,j-1)+                                &
     &                   grad(Iend,j  ))).gt.0.0_r8) THEN
                dVde=grad(Iend,j-1)
              ELSE
                dVde=grad(Iend,j  )
              END IF
              cff=MAX(dVdx*dVdx+dVde*dVde,eps)
              Cx=dVdt*dVdx
#ifdef RADIATION_2D
              Ce=MIN(cff,MAX(dVdt*dVde,-cff))
#else
              Ce=0.0_r8
#endif
#if defined CELERITY_WRITE && defined FORWARD_WRITE
              BOUNDARY(ng)%vbar_east_Cx(j)=Cx
              BOUNDARY(ng)%vbar_east_Ce(j)=Ce
              BOUNDARY(ng)%vbar_east_C2(j)=cff
#endif
              vbar(Iend+1,j,kout)=(cff * vbar(Iend + 1, j, kNow) + &
                                   & Cx * vbar(Iend  ,j,kout) - &
                                   & MAX(Ce,0.0_r8)*grad(Iend+1,j-1)-     &
     &                             MIN(Ce,0.0_r8)*grad(Iend+1,j  ))/    &
     &                            (cff+Cx)

              IF (LBC(ieast,isVbar,ng)%nudging) THEN
                vbar(Iend+1,j,kout)=vbar(Iend+1,j,kout)+                &
     &                              tau*(BOUNDARY(ng) % vbar_east(j) - &
                                         & vbar(Iend + 1, j, kNow))
              END IF
#ifdef MASKING
              vbar(Iend+1,j,kout)=vbar(Iend+1,j,kout)*                  &
     &                            GRID(ng)%vmask(Iend+1,j)
#endif
            END IF
          END DO
!
!  Eastern edge, Chapman boundary condition.
!
        ELSE IF (LBC(ieast,isVbar,ng)%Flather.or.                       &
     &           LBC(ieast,isVbar,ng)%reduced.or.                       &
     &           LBC(ieast,isVbar,ng)%Shchepetkin) THEN
          DO j=JstrV,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              cff=dt2d*0.5_r8*(GRID(ng)%pm(Iend,j-1)+                   &
     &                         GRID(ng)%pm(Iend,j  ))
              cff1=SQRT(g*0.5_r8*(GRID(ng) % h(Iend,j-1) + &
                                  & zeta(Iend, j - 1, kNow) + &
                                  & GRID(ng) % h(Iend,j  ) + &
                                  & zeta(Iend, j, kNow)))
              Cx=cff*cff1
              cff2=1.0_r8/(1.0_r8+Cx)
              vbar(Iend+1,j,kout)=cff2*(vbar(Iend + 1, j, kNow) + &
                                        & Cx * vbar(Iend,j,kout))
#ifdef MASKING
              vbar(Iend+1,j,kout)=vbar(Iend+1,j,kout)*                  &
     &                            GRID(ng)%vmask(Iend+1,j)
#endif
            END IF
          END DO
!
!  Eastern edge, clamped boundary condition.
!
        ELSE IF (LBC(ieast,isVbar,ng)%clamped) THEN
          DO j=JstrV,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              vbar(Iend+1,j,kout)=BOUNDARY(ng)%vbar_east(j)
#ifdef MASKING
              vbar(Iend+1,j,kout)=vbar(Iend+1,j,kout)*                  &
     &                            GRID(ng)%vmask(Iend+1,j)
#endif
            END IF
          END DO
!
!  Eastern edge, gradient boundary condition.
!
        ELSE IF (LBC(ieast,isVbar,ng)%gradient) THEN
          DO j=JstrV,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              vbar(Iend+1,j,kout)=vbar(Iend,j,kout)
#ifdef MASKING
              vbar(Iend+1,j,kout)=vbar(Iend+1,j,kout)*                  &
     &                            GRID(ng)%vmask(Iend+1,j)
#endif
            END IF
          END DO
!
!  Eastern edge, closed boundary condition: free slip (gamma2=1)  or
!                                           no   slip (gamma2=-1).
!
        ELSE IF (LBC(ieast,isVbar,ng)%closed) THEN
          IF (NSperiodic(ng)) THEN
            Jmin=JstrV
            Jmax=Jend
          ELSE
            Jmin=Jstr
            Jmax=JendR
          END IF
          DO j=Jmin,Jmax
            IF (LBC_apply(ng)%east(j)) THEN
              vbar(Iend+1,j,kout)=gamma2(ng)*vbar(Iend,j,kout)
#ifdef MASKING
              vbar(Iend+1,j,kout)=vbar(Iend+1,j,kout)*                  &
     &                            GRID(ng)%vmask(Iend+1,j)
#endif
            END IF
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (.not.(EWperiodic(ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jstr  )) THEN
            vbar(Istr-1,Jstr,kout)=0.5_r8*(vbar(Istr  ,Jstr  ,kout)+    &
     &                                     vbar(Istr-1,Jstr+1,kout))
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr  )) THEN
            vbar(Iend+1,Jstr,kout)=0.5_r8*(vbar(Iend  ,Jstr  ,kout)+    &
     &                                     vbar(Iend+1,Jstr+1,kout))
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
            vbar(Istr-1,Jend+1,kout)=0.5_r8*(vbar(Istr-1,Jend  ,kout)+  &
     &                                       vbar(Istr  ,Jend+1,kout))
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
            vbar(Iend+1,Jend+1,kout)=0.5_r8*(vbar(Iend+1,Jend  ,kout)+  &
     &                                      vbar(Iend  ,Jend+1,kout))
          END IF
        END IF
      END IF

#if defined WET_DRY
!
!-----------------------------------------------------------------------
!  Impose wetting and drying conditions.
!-----------------------------------------------------------------------
!
      IF (.not.EWperiodic(ng)) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=JstrV,Jend
            IF (LBC_apply(ng)%west(j).or.                               &
     &          LBC(iwest,isVbar,ng)%nested) THEN
              cff1=ABS(ABS(GRID(ng)%vmask_wet(Istr-1,j))-1.0_r8)
              cff2=0.5_r8+DSIGN(0.5_r8,vbar(Istr-1,j,kout))*            &
     &                    GRID(ng)%vmask_wet(Istr-1,j)
              cff=0.5_r8*GRID(ng)%vmask_wet(Istr-1,j)*cff1+             &
     &            cff2*(1.0_r8-cff1)
              vbar(Istr,j,kout)=vbar(Istr,j,kout)*cff
            END IF
          END DO
        END IF
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=JstrV,Jend
            IF (LBC_apply(ng)%east(j).or.                               &
     &          LBC(ieast,isVbar,ng)%nested) THEN
              cff1=ABS(ABS(GRID(ng)%vmask_wet(Iend+1,j))-1.0_r8)
              cff2=0.5_r8+DSIGN(0.5_r8,vbar(Iend+1,j,kout))*            &
     &                    GRID(ng)%vmask_wet(Iend+1,j)
              cff=0.5_r8*GRID(ng)%vmask_wet(Iend+1,j)*cff1+             &
     &            cff2*(1.0_r8-cff1)
              vbar(Iend+1,j,kout)=vbar(Iend+1,j,kout)*cff
            END IF
          END DO
        END IF
      END IF

      IF (.not.NSperiodic(ng)) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i).or.                              &
     &          LBC(isouth,isVbar,ng)%nested) THEN
              cff1=ABS(ABS(GRID(ng)%vmask_wet(i,Jstr))-1.0_r8)
              cff2=0.5_r8+DSIGN(0.5_r8,vbar(i,Jstr,kout))*              &
     &                    GRID(ng)%vmask_wet(i,Jstr)
              cff=0.5_r8*GRID(ng)%vmask_wet(i,Jstr)*cff1+               &
     &            cff2*(1.0_r8-cff1)
              vbar(i,Jstr,kout)=vbar(i,Jstr,kout)*cff
            END IF
          END DO
        END IF
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i).or.                              &
     &          LBC(inorth,isVbar,ng)%nested) THEN
              cff1=ABS(ABS(GRID(ng)%vmask_wet(i,Jend+1))-1.0_r8)
              cff2=0.5_r8+DSIGN(0.5_r8,vbar(i,Jend+1,kout))*            &
     &                    GRID(ng)%vmask_wet(i,Jend+1)
              cff=0.5_r8*GRID(ng)%vmask_wet(i,Jend+1)*cff1+             &
     &            cff2*(1.0_r8-cff1)
              vbar(i,Jend+1,kout)=vbar(i,Jend+1,kout)*cff
            END IF
          END DO
        END IF
      END IF
!
      IF (.not.(EWperiodic(ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          IF ((LBC_apply(ng)%south(Istr-1).and.                         &
     &         LBC_apply(ng)%west (Jstr  )).or.                         &
     &        (LBC(iwest,isVbar,ng)%nested.and.                         &
     &         LBC(isouth,isVbar,ng)%nested)) THEN
            cff1=ABS(ABS(GRID(ng)%vmask_wet(Istr-1,Jstr))-1.0_r8)
            cff2=0.5_r8+DSIGN(0.5_r8,vbar(Istr-1,Jstr,kout))*           &
     &                  GRID(ng)%vmask_wet(Istr-1,Jstr)
            cff=0.5_r8*GRID(ng)%vmask_wet(Istr-1,Jstr)*cff1+            &
     &          cff2*(1.0_r8-cff1)
            vbar(Istr-1,Jstr,kout)=vbar(Istr-1,Jstr,kout)*cff
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF ((LBC_apply(ng)%south(Iend+1).and.                         &
     &         LBC_apply(ng)%east (Jstr  )).or.                         &
     &        (LBC(ieast,isVbar,ng)%nested.and.                         &
     &         LBC(isouth,isVbar,ng)%nested)) THEN
            cff1=ABS(ABS(GRID(ng)%vmask_wet(Iend+1,Jstr))-1.0_r8)
            cff2=0.5_r8+DSIGN(0.5_r8,vbar(Iend+1,Jstr,kout))*           &
     &                  GRID(ng)%vmask_wet(Iend+1,Jstr)
            cff=0.5_r8*GRID(ng)%vmask_wet(Iend+1,Jstr)*cff1+            &
     &          cff2*(1.0_r8-cff1)
            vbar(Iend+1,Jstr,kout)=vbar(Iend+1,Jstr,kout)*cff
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF ((LBC_apply(ng)%north(Istr-1).and.                         &
     &         LBC_apply(ng)%west (Jend+1)).or.                         &
     &        (LBC(iwest,isVbar,ng)%nested.and.                         &
     &         LBC(inorth,isVbar,ng)%nested)) THEN
            cff1=ABS(ABS(GRID(ng)%vmask_wet(Istr-1,Jend+1))-1.0_r8)
            cff2=0.5_r8+DSIGN(0.5_r8,vbar(Istr-1,Jend+1,kout))*         &
     &                  GRID(ng)%vmask_wet(Istr-1,Jend+1)
            cff=0.5_r8*GRID(ng)%vmask_wet(Istr-1,Jend+1)*cff1+          &
     &          cff2*(1.0_r8-cff1)
            vbar(Istr-1,Jend+1,kout)=vbar(Istr-1,Jend+1,kout)*cff
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF ((LBC_apply(ng)%north(Iend+1).and.                         &
     &         LBC_apply(ng)%east (Jend+1)).or.                         &
     &        (LBC(ieast,isVbar,ng)%nested.and.                         &
     &         LBC(inorth,isVbar,ng)%nested)) THEN
            cff1=ABS(ABS(GRID(ng)%vmask_wet(Iend+1,Jend+1))-1.0_r8)
            cff2=0.5_r8+DSIGN(0.5_r8,vbar(Iend+1,Jend+1,kout))*         &
     &                  GRID(ng)%vmask_wet(Iend+1,Jend+1)
            cff=0.5_r8*GRID(ng)%vmask_wet(Iend+1,Jend+1)*cff1+          &
     &          cff2*(1.0_r8-cff1)
            vbar(Iend+1,Jend+1,kout)=vbar(Iend+1,Jend+1,kout)*cff
          END IF
        END IF
      END IF
#endif

      RETURN

      END SUBROUTINE v2dbc_tile
      END MODULE v2dbc_mod




