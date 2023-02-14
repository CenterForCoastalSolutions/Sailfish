
if SOLVE3D:

# !=======================================================================
# !  Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
# !    Licensed under a MIT/X style license                              !
# !    See License_ROMS.txt                           Hernan G. Arango   !
# !========================================== Alexander F. Shchepetkin ===
# !                                                                      !
# !  This subroutines computes vertical velocity (m/s) at W-points       !
# !  from the vertical mass flux (omega*hz/m*n).  This computation       !
# !  is done solely for output purposes.                                 !
# !                                                                      !
# !=======================================================================
!
!***********************************************************************
      SUBROUTINE wvelocity (ng, tile, Ninp)
!***********************************************************************
!
      USE mod_param
      USE mod_coupling
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Ninp
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL wvelocity_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     Ninp,                                        &
     &                     GRID(ng) % pm,                               &
     &                     GRID(ng) % pn,                               &
     &                     GRID(ng) % z_r,                              &
     &                     GRID(ng) % z_w,                              &
     &                     COUPLING(ng) % DU_avg1,                      &
     &                     COUPLING(ng) % DV_avg1,                      &
     &                     OCEAN(ng) % u,                               &
     &                     OCEAN(ng) % v,                               &
     &                     OCEAN(ng) % W,                               &
# if defined OMEGA_IMPLICIT
     &                     OCEAN(ng) % Wi,                              &
# endif
     &                     OCEAN(ng) % wvel)

      #
      # real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: vert
      #
      # real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wrk

# !-----------------------------------------------------------------------
# !  Compute "true" vertical velocity (m/s).
# !-----------------------------------------------------------------------


# !  In ROMS, the terrain-following vertical velocity, omega, is given by:
# !
# !         Hz * omega = w - d(z)/d(t) - div(z)
# !
# !  where w is the "true" vertical velocity and
# !
# !         div(z) = pm * u * d(z)/d(xi) + pn * v * d(z)/d(eta)
# !
# !  The vertical coordinate is a function of several parameter but only
# !  the free-surface is time dependent. However, in sediment applications
# !  with stratigraphy, the bathymetry (h) also evolves in time.
# !
# !  Exchange time-averaged fields.

    exchange_u2d_tile (ng, tile, LBi, UBi, LBj, UBj, DU_avg1)
    exchange_v2d_tile (ng, tile, LBi, UBi, LBj, UBj, DV_avg1)

# !  Compute contribution due to quasi-horizontal motions along
# !  S-coordinate surfaces:  (Ui + Vj)*GRADs(z).

    vert = u[:,:,:,Ninp]*dxRtoU(z_r, pm) + v[:,:,:,Ninp]*dyVtoR(z_r, pm)

# !  Compute contribution due to time tendency of the free-surface,
# !  d(zeta)/d(t), which is the vertical velocity at the free-surface
# !  and it is expressed in terms of barotropic mass flux divergence.
# !  Notice that it is divided by the total depth of the water column.
# !  This is needed because this contribution is linearly distributed
# !  throughout the water column by multiplying it by the distance from
# !  the bottom to the depth at which the vertical velocity is computed.


    cff1=3.0_r8/8.0_r8
    cff2=3.0_r8/4.0_r8
    cff3=1.0_r8/8.0_r8
    cff4=9.0_r8/16.0_r8
    cff5=1.0_r8/16.0_r8

    wrk(i, j) = -dxUtoR(DU_avg1(i, j) - DU_avg1(i + 1, j) + &
                 & DV_avg1(i, j) - DV_avg1(i, j + 1)) / &
    & (z_w(i, j, N(ng)) - z_w(i, j, 0))

      J_LOOP : DO j=Jstr,Jend
        DO i=Istr,Iend
          wrk(i,j)=(DU_avg1(i,j)-DU_avg1(i+1,j)+                        &
     &              DV_avg1(i,j)-DV_avg1(i,j+1))/                       &
     &             (z_w(i,j,N(ng))-z_w(i,j,0))
        END DO
!
!  Notice that a cubic interpolation is used to shift the "vert"
!  contribution from vertical RHO- to W-points.
!
        DO i=Istr,Iend
          slope=(z_r(i,j,1)-z_w(i,j,0))/                                &
     &          (z_r(i,j,2)-z_r(i,j,1))            ! extrapolation slope
          wvel(i,j,0)=cff1*(vert(i,j,1)-                                &
     &                      slope*(vert(i,j,2)-                         &
     &                             vert(i,j,1)))+                       &
     &                cff2*vert(i,j,1)-                                 &
     &                cff3*vert(i,j,2)
          wvel(i,j,1)=pm(i,j)*pn(i,j)*                                  &
     &                (W(i,j,1)+                                        &
# if defined OMEGA_IMPLICIT
     &                 Wi(i,j,1)+                                       &
# endif
     &                 wrk(i,j)*(z_w(i,j,1)-z_w(i,j,0)))+               &
     &                cff1*vert(i,j,1)+                                 &
     &                cff2*vert(i,j,2)-                                 &
     &                cff3*vert(i,j,3)
        END DO
        DO k=2,N(ng)-2
          DO i=Istr,Iend
            wvel(i,j,k)=pm(i,j)*pn(i,j)*                                &
     &                  (W(i,j,k)+                                      &
# if defined OMEGA_IMPLICIT
     &                   Wi(i,j,k)+                                     &
# endif
     &                   wrk(i,j)*(z_w(i,j,k)-z_w(i,j,0)))+             &
     &                  cff4*(vert(i,j,k  )+vert(i,j,k+1))-             &
     &                  cff5*(vert(i,j,k-1)+vert(i,j,k+2))
          END DO
        END DO
        DO i=Istr,Iend
          slope=(z_w(i,j,N(ng))-z_r(i,j,N(ng)  ))/                      &
     &          (z_r(i,j,N(ng))-z_r(i,j,N(ng)-1))  ! extrapolation slope
          wvel(i,j,N(ng))=pm(i,j)*pn(i,j)*                              &
     &                    wrk(i,j)*(z_w(i,j,N(ng))-z_w(i,j,0))+         &
     &                    cff1*(vert(i,j,N(ng))+                        &
     &                          slope*(vert(i,j,N(ng)  )-               &
     &                                 vert(i,j,N(ng)-1)))+             &
     &                    cff2*vert(i,j,N(ng)  )-                       &
     &                    cff3*vert(i,j,N(ng)-1)
          wvel(i,j,N(ng)-1)=pm(i,j)*pn(i,j)*                            &
     &                      (W(i,j,N(ng)-1)+                            &
# if defined OMEGA_IMPLICIT
     &                       Wi(i,j,N(ng)-1)+                           &
# endif
     &                       wrk(i,j)*(z_w(i,j,N(ng)-1)-z_w(i,j,0)))+   &
     &                      cff1*vert(i,j,N(ng)  )+                     &
     &                      cff2*vert(i,j,N(ng)-1)-                     &
     &                      cff3*vert(i,j,N(ng)-2)
        END DO
      END DO J_LOOP
!
!  Set lateral boundary conditions.
!
      CALL bc_w3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 0, N(ng),                   &
     &                  wvel)
# ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    wvel)
# endif

      RETURN
      END SUBROUTINE wvelocity_tile
#endif
      END MODULE wvelocity_mod