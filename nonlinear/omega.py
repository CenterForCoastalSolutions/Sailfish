
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

# ifdef EMINUSP_SHIMA
        fac2 = 1.0/N
# endif

# if defined OMEGA_IMPLICIT
      cmnx_ratio = amin/amax
      cutoff = 2.0 - amin/amax
      r4cmx = 0.25/(1.0 - amin/amax)
# endif
      DO j=Jstr,Jend
        DO i=Istr,Iend
          W(i,j,0)=0.0_r8

        END DO

        DO k=1,N(ng)
          DO i=Istr,Iend
            W(i,j,k) = W(i,j,k-1) -                                        &

# ifdef EMINUSP_SHIMA
     &                stflux(i,j,isalt)*omn(i,j)*fac2 -                  &
# endif

     &               (Huon(i+1,j,k)-Huon(i,j,k) + Hvom(i,j+1,k)-Hvom(i,j,k))
# if defined OMEGA_IMPLICIT


        # Compute the horizontal Courant number
        # ---------------------------------------------------------

            Cu_adv(i,k) =
                     MAX(Huon(i+1,j,k),0.0) - MIN(Huon(i,j,k),0.0) +
                     MAX(Hvom(i,j+1,k),0.0) - MIN(Hvom(i,j,k),0.0)
# endif
          END DO
        END DO

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

# ifdef EMINUSP_SHIMA
     &                       stflux(ii,jj,isalt)*omn(ii,jj)*fac2-       &
# endif
     &                       (Huon(ii+1,jj,k)-Huon(ii,jj,k)+            &
     &                        Hvom(ii,jj+1,k)-Hvom(ii,jj,k))+           &
     &                       SOURCES(ng)%Qsrc(is,k)
                END DO
              END IF
            END IF
          END DO
        END IF
!
        DO i=Istr,Iend
          wrk(i)=W(i,j,N(ng))/(z_w(i,j,N(ng))-z_w(i,j,0))
# if defined OMEGA_IMPLICIT
          Cu_adv(i,0)=dt(ng)*pm(i,j)*pn(i,j)
# endif
        END DO

        # In order to insure zero vertical velocity at the free-surface,
        # subtract the vertical velocities of the moving S-coordinates
        # isosurfaces. These isosurfaces are proportional to d(zeta)/d(t).
        # The proportionally coefficients are a linear function of the
        # S-coordinate with zero value at the bottom (k=0) and unity at
        # the free-surface (k=N).

        DO k=N(ng)-1,1,-1
          DO i=Istr,Iend
            W(i,j,k)=W(i,j,k)-                                          &

     &               wrk(i)*(z_w(i,j,k)-z_w(i,j,0))



# if defined OMEGA_IMPLICIT

        # Determine implicit part Wi of vertical advection.
        # W  becomes the explicit part We.

            Wi(i,j,k)=W(i,j,k)
            IF (Wi(i,j,k).ge.0.0_r8) THEN        ! Three different variants
              c2d=Cu_adv(i,k)			 ! for computing 2D Courant
              dh=z_w(i,j,k)-z_w(i,j,k-1)	 ! number at the interface:
            ELSE				 ! (1) use value from the
              c2d=Cu_adv(i,k+1)			 !     grid box upstream in
              dh=z_w(i,j,k+1)-z_w(i,j,k)	 !     vertical direction;
            END IF


            cw_max=amax*dh-c2d*Cu_adv(i,0)  ! compare vertical displacement
            IF (cw_max.ge.0.0_r8) THEN      ! to dz*amax. Partition W into
              cw_max2=cw_max*cw_max         ! Wi and We.
              cw_min=cw_max*cmnx_ratio
              cw=ABS(Wi(i,j,k))*Cu_adv(i,0)
              IF (cw.le.cw_min) THEN
                cff=cw_max2
              ELSE IF (cw.le.cutoff*cw_max) THEN
                cff=cw_max2+r4cmx*(cw-cw_min)**2
              ELSE
                cff=cw_max*cw
              END IF

              W(i,j,k)=cw_max2*Wi(i,j,k)/cff
              Wi(i,j,k)=Wi(i,j,k)-W(i,j,k)
            ELSE                            ! All the displacement is
              W(i,j,k)=0.0_r8               ! greater than amax*dz, so
            END IF                          ! keep it all into Wi.
# endif
          END DO
        END DO
        DO i=Istr,Iend
          W(i,j,N(ng))=0.0_r8
        END DO
      END DO

    # Set lateral boundary conditions.

      CALL bc_w3d_tile (LBi, UBi, LBj, UBj, 0, N(ng), W)
# if defined OMEGA_IMPLICIT
      CALL bc_w3d_tile (LBi, UBi, LBj, UBj, 0, N(ng), Wi)
# endif



# endif

def scale_omega ( LBk, UBk,pm, pn, W, Wscl):

    # Scale omega vertical velocity to m/s.

    DO k=LBk,UBk
      DO j=JstrR,JendR
        DO i=IstrR,IendR
            Wscl(i,j,k) = W(i,j,k)*pm(i,j)*pn(i,j)
        END DO
      END DO
    END DO
