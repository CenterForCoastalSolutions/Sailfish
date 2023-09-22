
# This routine computes S-coordinate vertical velocity (m^3/s),
#
#                 W = [Hz/(m*n)]*omega,
#
# diagnostically at horizontal RHO-points and vertical W-points.
#
# Added implicit vertical adveciton froman adaptive, Courant-number-dependent implicit scheme for vertical
# advection in oceanic modeling, Alexander F. Shchepetkin, pp 38-69.





    # Vertically integrate horizontal mass flux divergence.
    # ------------------------------------------------------------------------

    # Starting with zero vertical velocity at the bottom, integrate
    # from the bottom (k=0) to the free-surface (k=N).  The w(:,:,N(ng))
    # contains the vertical velocity at the free-surface, d(zeta)/d(t).
    # Notice that barotropic mass flux divergence is not used directly.

    W[0,:,:] = 0.0
    for k = in range(1,N):
        W[k,:,:] = W[k-1,:,:] - Dx(Huon) + Dy(Hvom(i,j+1,k))

      # DO j=Jstr,Jend
      #   DO i=Istr,Iend
      #     W(i,j,0)=0.0_r8
      #
      #   END DO
      #
      #   DO k=1,N(ng)
      #     DO i=Istr,Iend
      #       W(i,j,k) = W(i,j,k-1) - (Huon(i+1,j,k)-Huon(i,j,k) + Hvom(i,j+1,k)-Hvom(i,j,k))
      #
      #     END DO
      #   END DO

    pass

        #  Apply mass point sources (volume vertical influx), if any.
        #
        #  Overwrite W(Isrc,Jsrc,k) with the same divergence of Huon,Hvom as
        #  above but add in point source Qsrc(k) and reaccumulate the vertical
        #  sum to obtain the correct net Qbar given in user input - J. Levin
        #  (Jupiter Intelligence Inc.) and J. Wilkin
        #
        #    Dsrc(is) = 2,  flow across grid cell w-face (positive or negative)

        IF (LwSrc(ng)) THEN
          DO is=1,Nsrc(ng)
            IF (INT(SOURCES(ng)%Dsrc(is)).eq.2) THEN
              ii=SOURCES(ng)%Isrc(is)
              jj=SOURCES(ng)%Jsrc(is)
              IF (((IstrR.le.ii).and.(ii.le.IendR)).and.                &
     &            ((JstrR.le.jj).and.(jj.le.JendR)).and.                &
     &            (j.eq.jj)) THEN

                DO k=1,N(ng)
                  W(ii,jj,k)=W(ii,jj,k-1)-                              &


     &                       (Huon(ii+1,jj,k)-Huon(ii,jj,k)+            &
     &                        Hvom(ii,jj+1,k)-Hvom(ii,jj,k))+           &
     &                       SOURCES(ng)%Qsrc(is,k)
                END DO
              END IF
            END IF
          END DO
        END IF

        DO i=Istr,Iend
          wrk(i)=W(i,j,N(ng))/(z_w(i,j,N(ng))-z_w(i,j,0))

        END DO


        # In order to ensure zero vertical velocity at the free-surface,
        # subtract the vertical velocities of the moving S-coordinates
        # isosurfaces. These isosurfaces are proportional to d(zeta)/d(t).
        # The proportionally coefficients are a linear function of the
        # S-coordinate with zero value at the bottom (k=0) and unity at
        # the free-surface (k=N).

        DO k=N(ng)-1,1,-1
          DO i=Istr,Iend
            W(i,j,k) = W(i,j,k) - wrk(i)*(z_w(i,j,k) - z_w(i,j,0))
          END DO
        END DO
        DO i=Istr,Iend
          W(i,j,N(ng))=0.0_r8
        END DO
      END DO

      # Set lateral boundary conditions.
      CALL bc_w3d_tile (LBi, UBi, LBj, UBj, 0, N(ng), W)




def scale_omega ( LBk, UBk,pm, pn, W, Wscl):

    # Scale omega vertical velocity to m/s.

    DO k=LBk,UBk
      DO j=JstrR,JendR
        DO i=IstrR,IendR
            Wscl(i,j,k) = W(i,j,k)*pm(i,j)*pn(i,j)
        END DO
      END DO
    END DO
