                                               !
# This subroutine time-steps the nonlinear  horizontal  momentum
# equations. The vertical viscosity terms are time-stepped using
# an implicit algorithm.


    # Time step momentum equation in the XI-direction.
    # -----------------------------------------------------------------------
    DO j=Jstr,Jend
        DO i=IstrU,Iend

          AK(i,0) = 0.5*(Akv(i-1,j,0)+ Akv(i  ,j,0))
          DO k=1,N(ng)
            AK (i,k) = 0.5*(Akv(i-1,j,k) + Akv(i,j,k))
            Hzk(i,k) = 0.5*(Hz (i-1,j,k) + Hz (i,j,k))

# if defined SPLINES_VVISC || defined DIAGNOSTICS_UV
            oHz(i,k) = 1.0/Hzk(i,k)
# endif
          END DO
        END DO


    # Time step right-hand-side terms.

        IF (iic(ng).eq.ntfirst(ng)) THEN
          cff=0.25*dt
        ELSE IF (iic(ng).eq.(ntfirst(ng)+1)) THEN
          cff=0.25*dt*3.0/2.0
        ELSE
          cff=0.25*dt*23.0/12.0
        END IF

        DO i=IstrU,Iend
          DC(i,0)=cff*(pm(i,j) + pm(i-1,j))*(pn(i,j) + pn(i-1,j))
        END DO

        DO k=1,N
          DO i=IstrU,Iend
            u(i,j,k,nnew)=u(i,j,k,nnew) + DC(i,0)*ru(i,j,k,nrhs)
          END DO
        END DO



    # Compute off-diagonal coefficients [lambda*dt*Akv/Hz] for the
    # implicit vertical viscosity term at horizontal U-points and
    # vertical W-points.

        cff = -lambda*dt/0.5
        DO k=1,N(ng)-1
          DO i=IstrU,Iend
            cff1 = 1.0/(z_r(i,j,k+1) + z_r(i-1,j,k+1) - z_r(i,j,k  ) - z_r(i-1,j,k  ))
            FC(i,k)=cff*cff1*AK(i,k)
          END DO
        END DO

        DO i=IstrU,Iend
          FC(i,0)=0.0
          FC(i,N(ng))=0.0
        END DO

    # Solve the tridiagonal system.

        DO k=1,N(ng)
          DO i=IstrU,Iend
            DC(i,k)=u(i,j,k,nnew)
            BC(i,k)=Hzk(i,k)-FC(i,k)-FC(i,k-1)
          END DO
        END DO
        DO i=IstrU,Iend
          cff=1.0/BC(i,1)
          CF(i,1)=cff*FC(i,1)
          DC(i,1)=cff*DC(i,1)
        END DO
        DO k=2,N(ng)-1
          DO i=IstrU,Iend
            cff=1.0_r8/(BC(i,k)-FC(i,k-1)*CF(i,k-1))
            CF(i,k)=cff*FC(i,k)
            DC(i,k)=cff*(DC(i,k)-FC(i,k-1)*DC(i,k-1))
          END DO
        END DO

        # Compute new solution by back substitution.

        DO i=IstrU,Iend

          DC(i,N(ng))=(DC(i,N(ng))-FC(i,N(ng)-1)*DC(i,N(ng)-1))/        &
     &                (BC(i,N(ng))-FC(i,N(ng)-1)*CF(i,N(ng)-1))
          u(i,j,N(ng),nnew)=DC(i,N(ng))

        END DO
        DO k=N(ng)-1,1,-1
          DO i=IstrU,Iend

            DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
            u(i,j,k,nnew)=DC(i,k)

          END DO
        END DO
# endif

# if defined OMEGA_IMPLICIT
!
!  Implicit vertical advection part for u.
!
        DO i=IstrU,Iend
          WK(i,0)=0.5_r8*(Wi(i-1,j,0)+                                  &
     &                    Wi(i  ,j,0))
          DO k=1,N(ng)
            WK(i,k)=0.5_r8*(Wi(i-1,j,k)+                                &
     &                      Wi(i  ,j,k))
            Hzk(i,k)=0.5_r8*(Hz(i-1,j,k)+                               &
     &                       Hz(i  ,j,k))
          END DO
        END DO
!
!  Compute off-diagonal coefficients [dt*Wi*pm*pn] for the
!  implicit vertical viscosity term at horizontal U-points and
!  vertical W-points.
!
        cff=dt(ng)
        DO k=1,N(ng)-1
          DO i=IstrU,Iend
            cff1=cff/(on_u(i,j)*om_u(i,j))
            FCmax(i,k)=MAX(WK(i,k),0.0_r8)*cff1
            FCmin(i,k)=MIN(WK(i,k),0.0_r8)*cff1
          END DO
        END DO
        DO i=IstrU,Iend
          FCmax(i,0)=0.0_r8
          FCmin(i,0)=0.0_r8
          FCmax(i,N(ng))=0.0_r8
          FCmin(i,N(ng))=0.0_r8
        END DO
!
!  Solve the tridiagonal system.
!
        DO k=1,N(ng)
          DO i=IstrU,Iend
            BC(i,k)=Hzk(i,k)+FCmax(i,k)-FCmin(i,k-1)
            DC(i,k)=u(i,j,k,nnew)*Hzk(i,k)
          END DO
        END DO
        DO i=IstrU,Iend
          cff=1.0_r8/BC(i,1)
          CF(i,1)=cff*FCmin(i,1)
          DC(i,1)=cff*DC(i,1)
        END DO
        DO k=2,N(ng)-1
          DO i=IstrU,Iend
            cff=1.0_r8/(BC(i,k)+FCmax(i,k-1)*CF(i,k-1))
            CF(i,k)=cff*FCmin(i,k)
            DC(i,k)=cff*(DC(i,k)+FCmax(i,k-1)*DC(i,k-1))
          END DO
        END DO
!
!  Compute new solution by back substitution.
!
        DO i=IstrU,Iend
#  ifdef DIAGNOSTICS_UV
          cff1=u(i,j,N(ng),nnew)
#  endif
          cff=1.0_r8/(BC(i,N(ng))+FCmax(i,N(ng)-1)*CF(i,N(ng)-1))
          DC(i,N(ng))=cff*(DC(i,N(ng))+                                 &
     &                     FCmax(i,N(ng)-1)*DC(i,N(ng)-1))
          u(i,j,N(ng),nnew)=DC(i,N(ng))
#  ifdef DIAGNOSTICS_UV
          DiaRU(i,j,N(ng),nrhs,M3vadv)=DiaRU(i,j,N(ng),nrhs,M3vadv)+    &
     &                                 u(i,j,N(ng),nnew)-cff1
#  endif
        END DO
!
        DO k=N(ng)-1,1,-1
          DO i=IstrU,Iend
#  ifdef DIAGNOSTICS_UV
            cff1=u(i,j,k,nnew)
#  endif
            DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
            u(i,j,k,nnew)=DC(i,k)
#  ifdef DIAGNOSTICS_UV
            DiaRU(i,j,k,nrhs,M3vadv)=DiaRU(i,j,k,nrhs,M3vadv)+          &
     &                               u(i,j,k,nnew)-cff1
#  endif
          END DO
        END DO
# endif
!
!  Replace INTERIOR POINTS incorrect vertical mean with more accurate
!  barotropic component, ubar=DU_avg1/(D*on_u). Recall that, D=CF(:,0).
!
        DO i=IstrU,Iend
          CF(i,0)=Hzk(i,1)
          DC(i,0)=u(i,j,1,nnew)*Hzk(i,1)

        END DO
        DO k=2,N(ng)
          DO i=IstrU,Iend
            CF(i,0)=CF(i,0)+Hzk(i,k)
            DC(i,0)=DC(i,0)+u(i,j,k,nnew)*Hzk(i,k)

          END DO
        END DO
        DO i=IstrU,Iend
          cff1=1.0_r8/(CF(i,0)*on_u(i,j))
          DC(i,0)=(DC(i,0)*on_u(i,j)-DU_avg1(i,j))*cff1      ! recursive

        END DO
!
!  Couple and update new solution.
!
        DO k=1,N(ng)
          DO i=IstrU,Iend
            u(i,j,k,nnew)=u(i,j,k,nnew)-DC(i,0)
# ifdef MASKING
            u(i,j,k,nnew)=u(i,j,k,nnew)*umask(i,j)
# endif

          END DO
        END DO

# if defined DIAGNOSTICS_UV && defined MASKING
        DO k=1,N(ng)
          DO i=IstrU,Iend
            DO idiag=1,NDM3d
              DiaU3wrk(i,j,k,idiag)=DiaU3wrk(i,j,k,idiag)*umask(i,j)
            END DO
          END DO
        END DO
# endif
!
!-----------------------------------------------------------------------
!  Time step momentum equation in the ETA-direction.
!-----------------------------------------------------------------------
!
        IF (j.ge.JstrV) THEN
          DO i=Istr,Iend
            AK(i,0)=0.5_r8*(Akv(i,j-1,0)+                               &
     &                      Akv(i,j  ,0))
            DO k=1,N(ng)
              AK(i,k)=0.5_r8*(Akv(i,j-1,k)+                             &
     &                        Akv(i,j  ,k))
              Hzk(i,k)=0.5_r8*(Hz(i,j-1,k)+                             &
     &                         Hz(i,j  ,k))
# if defined SPLINES_VVISC || defined DIAGNOSTICS_UV
              oHz(i,k)=1.0_r8/Hzk(i,k)
# endif
            END DO
          END DO
!
!  Time step right-hand-side terms.
!
          IF (iic(ng).eq.ntfirst(ng)) THEN
            cff=0.25_r8*dt(ng)
          ELSE IF (iic(ng).eq.(ntfirst(ng)+1)) THEN
            cff=0.25_r8*dt(ng)*3.0_r8/2.0_r8
          ELSE
            cff=0.25_r8*dt(ng)*23.0_r8/12.0_r8
          END IF
          DO i=Istr,Iend
            DC(i,0)=cff*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
          END DO
          DO k=1,N(ng)
            DO i=Istr,Iend
              v(i,j,k,nnew)=v(i,j,k,nnew)+DC(i,0)*rv(i,j,k,nrhs)
# ifdef SPLINES_VVISC
              v(i,j,k,nnew)=v(i,j,k,nnew)*oHz(i,k)
# endif

            END DO
          END DO

# ifdef SPLINES_VVISC
!
!  Use conservative, parabolic spline reconstruction of vertical
!  viscosity derivatives.  Then, time step vertical viscosity term
!  implicitly by solving a tridiagonal system.
!
          cff1=1.0_r8/6.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FC(i,k)=cff1*Hzk(i,k  )-dt(ng)*AK(i,k-1)*oHz(i,k  )
              CF(i,k)=cff1*Hzk(i,k+1)-dt(ng)*AK(i,k+1)*oHz(i,k+1)
            END DO
          END DO
          DO i=Istr,Iend
            CF(i,0)=0.0_r8
            DC(i,0)=0.0_r8
          END DO
!
!  LU decomposition and forward substitution.
!
          cff1=1.0_r8/3.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              BC(i,k)=cff1*(Hzk(i,k)+Hzk(i,k+1))+                       &
     &                dt(ng)*AK(i,k)*(oHz(i,k)+oHz(i,k+1))
              cff=1.0_r8/(BC(i,k)-FC(i,k)*CF(i,k-1))
              CF(i,k)=cff*CF(i,k)
              DC(i,k)=cff*(v(i,j,k+1,nnew)-v(i,j,k,nnew)-               &
     &                     FC(i,k)*DC(i,k-1))
            END DO
          END DO
!
!  Backward substitution.
!
          DO i=Istr,Iend
            DC(i,N(ng))=0.0_r8
          END DO
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
              DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
            END DO
          END DO
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              DC(i,k)=DC(i,k)*AK(i,k)
              cff=dt(ng)*oHz(i,k)*(DC(i,k)-DC(i,k-1))
              v(i,j,k,nnew)=v(i,j,k,nnew)+cff

            END DO
          END DO
# else
!
!  Compute off-diagonal coefficients [lambda*dt*Akv/Hz] for the
!  implicit vertical viscosity term at horizontal V-points and
!  vertical W-points.
!
          cff=-lambda*dt(ng)/0.5_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              cff1=1.0_r8/(z_r(i,j,k+1)+z_r(i,j-1,k+1)-                 &
     &                     z_r(i,j,k  )-z_r(i,j-1,k  ))
              FC(i,k)=cff*cff1*AK(i,k)
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,0)=0.0_r8
            FC(i,N(ng))=0.0_r8
          END DO
!
!  Solve the tridiagonal system.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              DC(i,k)=v(i,j,k,nnew)
              BC(i,k)=Hzk(i,k)-FC(i,k)-FC(i,k-1)
            END DO
          END DO
          DO i=Istr,Iend
            cff=1.0_r8/BC(i,1)
            CF(i,1)=cff*FC(i,1)
            DC(i,1)=cff*DC(i,1)
          END DO
          DO k=2,N(ng)-1
            DO i=Istr,Iend
              cff=1.0_r8/(BC(i,k)-FC(i,k-1)*CF(i,k-1))
              CF(i,k)=cff*FC(i,k)
              DC(i,k)=cff*(DC(i,k)-FC(i,k-1)*DC(i,k-1))
            END DO
          END DO
!
!  Compute new solution by back substitution.
!
          DO i=Istr,Iend
            DC(i,N(ng))=(DC(i,N(ng))-FC(i,N(ng)-1)*DC(i,N(ng)-1))/      &
     &                  (BC(i,N(ng))-FC(i,N(ng)-1)*CF(i,N(ng)-1))
            v(i,j,N(ng),nnew)=DC(i,N(ng))

          END DO
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend

              DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
              v(i,j,k,nnew)=DC(i,k)

            END DO
          END DO
# endif

# if defined OMEGA_IMPLICIT
!
!  Implicit vertical advection part for v.
!
          DO i=Istr,Iend
            WK(i,0)=0.5_r8*(Wi(i,j-1,0)+                                &
     &                      Wi(i,j  ,0))
            DO k=1,N(ng)
              WK(i,k)=0.5_r8*(Wi(i,j-1,k)+                              &
     &                        Wi(i,j  ,k))
              Hzk(i,k)=0.5_r8*(Hz(i,j-1,k)+                             &
     &                         Hz(i,j  ,k))
            END DO
          END DO
!
!  Compute off-diagonal coefficients [dt*Wi*pm*pn] for the
!  implicit vertical viscosity term at horizontal V-points and
!  vertical W-points.
!
          cff=dt(ng)
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              cff1=cff/(on_v(i,j)*om_v(i,j))
              FCmax(i,k)=MAX(WK(i,k),0.0_r8)*cff1
              FCmin(i,k)=MIN(WK(i,k),0.0_r8)*cff1
            END DO
          END DO
          DO i=Istr,Iend
            FCmax(i,0)=0.0_r8
            FCmin(i,0)=0.0_r8
            FCmax(i,N(ng))=0.0_r8
            FCmin(i,N(ng))=0.0_r8
          END DO
!
!  Solve the tridiagonal system.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              BC(i,k)=Hzk(i,k)+FCmax(i,k)-FCmin(i,k-1)
              DC(i,k)=v(i,j,k,nnew)*Hzk(i,k)
            END DO
          END DO
          DO i=Istr,Iend
            cff=1.0_r8/BC(i,1)
            CF(i,1)=cff*FCmin(i,1)
            DC(i,1)=cff*DC(i,1)
          END DO
          DO k=2,N(ng)-1
            DO i=Istr,Iend
              cff=1.0_r8/(BC(i,k)+FCmax(i,k-1)*CF(i,k-1))
              CF(i,k)=cff*FCmin(i,k)
              DC(i,k)=cff*(DC(i,k)+FCmax(i,k-1)*DC(i,k-1))
            END DO
          END DO
!
!  Compute new solution by back substitution.
!
          DO i=Istr,Iend
            cff=1.0_r8/(BC(i,N(ng))+FCmax(i,N(ng)-1)*CF(i,N(ng)-1))
            DC(i,N(ng))=cff*(DC(i,N(ng))+                               &
     &                       FCmax(i,N(ng)-1)*DC(i,N(ng)-1))
            v(i,j,N(ng),nnew)=DC(i,N(ng))
          END DO
!
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
              DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
              v(i,j,k,nnew)=DC(i,k)
            END DO
          END DO
# endif
!
!  Replace INTERIOR POINTS incorrect vertical mean with more accurate
!  barotropic component, vbar=DV_avg1/(D*om_v). Recall that, D=CF(:,0).
!
          DO i=Istr,Iend
            CF(i,0)=Hzk(i,1)
            DC(i,0)=v(i,j,1,nnew)*Hzk(i,1)

          END DO
          DO k=2,N(ng)
            DO i=Istr,Iend
              CF(i,0)=CF(i,0)+Hzk(i,k)
              DC(i,0)=DC(i,0)+v(i,j,k,nnew)*Hzk(i,k)

            END DO
          END DO
          DO i=Istr,Iend
            cff1=1.0_r8/(CF(i,0)*om_v(i,j))
            DC(i,0)=(DC(i,0)*om_v(i,j)-DV_avg1(i,j))*cff1    ! recursive

          END DO
!
!  Couple and update new solution.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              v(i,j,k,nnew)=v(i,j,k,nnew)-DC(i,0)
# ifdef MASKING
              v(i,j,k,nnew)=v(i,j,k,nnew)*vmask(i,j)
# endif

            END DO
          END DO


        END IF
      END DO

!
!-----------------------------------------------------------------------
!  Set lateral boundary conditions.
!-----------------------------------------------------------------------
!
      CALL u3dbc_tile (ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj, N(ng),                       &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
     &                 nstp, nnew,                                      &
     &                 u)
      CALL v3dbc_tile (ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj, N(ng),                       &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
     &                 nstp, nnew,                                      &
     &                 v)
!
!-----------------------------------------------------------------------
!  Apply momentum transport point sources (like river runoff), if any.
!-----------------------------------------------------------------------
!
      IF (LuvSrc(ng)) THEN
        DO is=1,Nsrc(ng)
          i=SOURCES(ng)%Isrc(is)
          j=SOURCES(ng)%Jsrc(is)
          IF (((IstrR.le.i).and.(i.le.IendR)).and.                      &
     &        ((JstrR.le.j).and.(j.le.JendR))) THEN
            IF (INT(SOURCES(ng)%Dsrc(is)).eq.0) THEN
              DO k=1,N(ng)
                cff1=1.0_r8/(on_u(i,j)*                                 &
     &                       0.5_r8*(z_w(i-1,j,k)-z_w(i-1,j,k-1)+       &
     &                               z_w(i  ,j,k)-z_w(i  ,j,k-1)))
                u(i,j,k,nnew)=SOURCES(ng)%Qsrc(is,k)*cff1
              END DO
            ELSE IF (INT(SOURCES(ng)%Dsrc(is)).eq.1) THEN
              DO k=1,N(ng)
                cff1=1.0_r8/(om_v(i,j)*                                 &
     &                       0.5_r8*(z_w(i,j-1,k)-z_w(i,j-1,k-1)+       &
     &                               z_w(i,j  ,k)-z_w(i,j  ,k-1)))
                v(i,j,k,nnew)=SOURCES(ng)%Qsrc(is,k)*cff1
              END DO
            END IF
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Couple 2D and 3D momentum equations.
!-----------------------------------------------------------------------
!
!  Couple velocity component in the XI-direction.
!
      DO j=JstrT,JendT
        DO i=IstrP,IendT
          DC(i,0)=0.0_r8
          CF(i,0)=0.0_r8

          FC(i,0)=0.0_r8
        END DO
!
!  Compute thicknesses of U-boxes DC(i,1:N), total depth of the water
!  column DC(i,0), and incorrect vertical mean CF(i,0).  Notice that
!  barotropic component is replaced with its fast-time averaged
!  values.
!
        DO k=1,N(ng)
          DO i=IstrP,IendT
            cff=0.5_r8*on_u(i,j)
            DC(i,k)=cff*(Hz(i,j,k)+Hz(i-1,j,k))
            DC(i,0)=DC(i,0)+DC(i,k)
            CF(i,0)=CF(i,0)+                                            &
     &              DC(i,k)*u(i,j,k,nnew)

          END DO
        END DO
        DO i=IstrP,IendT
          cff1=DC(i,0)                                  ! intermediate

          DC(i,0)=1.0_r8/DC(i,0)                        ! recursive
          CF(i,0)=DC(i,0)*(CF(i,0)-DU_avg1(i,j))        ! recursive

          ubar(i,j,1)=DC(i,0)*DU_avg1(i,j)

          ubar(i,j,2)=ubar(i,j,1)

        END DO

!
!  Replace only BOUNDARY POINTS incorrect vertical mean with more
!  accurate barotropic component, ubar=DU_avg1/(D*on_u). Recall that,
!  D=CF(:,0).
!
!  NOTE:  Only the BOUNDARY POINTS need to be replaced. Avoid redundant
!         update in the interior again for computational purposes which
!         will not affect the nonlinear code.  However, the adjoint
!         code is wrong because the interior solution is corrected
!         twice. The replacement is avoided if the boundary edge is
!         periodic. The J-loop is pipelined, so we need to do a special
!         test for the southern and northern domain edges.
!
        IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            DO k=1,N(ng)
              u(Istr,j,k,nnew)=u(Istr,j,k,nnew)-CF(Istr,0)
# ifdef MASKING
              u(Istr,j,k,nnew)=u(Istr,j,k,nnew)*                        &
     &                         umask(Istr,j)
# endif

            END DO
          END IF
        END IF
!
        IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO k=1,N(ng)
              u(Iend+1,j,k,nnew)=u(Iend+1,j,k,nnew)-CF(Iend+1,0)
# ifdef MASKING
              u(Iend+1,j,k,nnew)=u(Iend+1,j,k,nnew)*                    &
     &                           umask(Iend+1,j)
# endif

            END DO
          END IF
        END IF
!
        IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
          IF (j.eq.0) THEN                           ! southern boundary
            DO k=1,N(ng)                             ! J-loop pipelined
              DO i=IstrU,Iend
                u(i,j,k,nnew)=u(i,j,k,nnew)-CF(i,0)
# ifdef MASKING
                u(i,j,k,nnew)=u(i,j,k,nnew)*                            &
     &                        umask(i,j)
# endif

              END DO
            END DO
          END IF
        END IF
!
        IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
          IF (j.eq.Mm(ng)+1) THEN                    ! northern boundary
            DO k=1,N(ng)                             ! J-loop pipelined
              DO i=IstrU,Iend
                u(i,j,k,nnew)=u(i,j,k,nnew)-CF(i,0)
# ifdef MASKING
                u(i,j,k,nnew)=u(i,j,k,nnew)*                            &
     &                        umask(i,j)
# endif

              END DO
            END DO
          END IF
        END IF
!
!  Compute correct mass flux, Hz*u/n.
!
        DO k=N(ng),1,-1
          DO i=IstrP,IendT
            Huon(i,j,k)=0.5_r8*(Huon(i,j,k)+u(i,j,k,nnew)*DC(i,k))

          END DO
        END DO
        DO i=IstrP,IendT
          FC(i,0)=DC(i,0)*(FC(i,0)-DU_avg2(i,j))        ! recursive
        END DO
        DO k=1,N(ng)
          DO i=IstrP,IendT
            Huon(i,j,k)=Huon(i,j,k)-DC(i,k)*FC(i,0)
          END DO
        END DO
!
!  Couple velocity component in the ETA-direction.
!
        IF (j.ge.Jstr) THEN
          DO i=IstrT,IendT
            DC(i,0)=0.0_r8
            CF(i,0)=0.0_r8

            FC(i,0)=0.0_r8
          END DO
!
!  Compute thicknesses of V-boxes DC(i,1:N), total depth of the water
!  column DC(i,0), and incorrect vertical mean CF(i,0).  Notice that
!  barotropic component is replaced with its fast-time averaged
!  values.
!
          DO k=1,N(ng)
            DO i=IstrT,IendT
              cff=0.5_r8*om_v(i,j)
              DC(i,k)=cff*(Hz(i,j,k)+Hz(i,j-1,k))
              DC(i,0)=DC(i,0)+DC(i,k)
              CF(i,0)=CF(i,0)+                                          &
     &                DC(i,k)*v(i,j,k,nnew)

            END DO
          END DO
          DO i=IstrT,IendT
            cff1=DC(i,0)                                 ! Intermediate

            DC(i,0)=1.0_r8/DC(i,0)                       ! recursive
            CF(i,0)=DC(i,0)*(CF(i,0)-DV_avg1(i,j))       ! recursive

            vbar(i,j,1)=DC(i,0)*DV_avg1(i,j)

            vbar(i,j,2)=vbar(i,j,1)

          END DO

!
!  Replace only BOUNDARY POINTS incorrect vertical mean with more
!  accurate barotropic component, vbar=DV_avg1/(D*om_v).  Recall that,
!  D=CF(:,0).
!
!  NOTE:  Only the BOUNDARY POINTS need to be replaced. Avoid redundant
!         update in the interior again for computational purposes which
!         will not affect the nonlinear code.  However, the adjoint
!         code is wrong because the interior solution is corrected
!         twice. The replacement is avoided if the boundary edge is
!         periodic. The J-loop is pipelined, so we need to do a special
!         test for the southern and northern domain edges.
!
          IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
            IF (DOMAIN(ng)%Western_Edge(tile)) THEN
              DO k=1,N(ng)
                v(Istr-1,j,k,nnew)=v(Istr-1,j,k,nnew)-CF(Istr-1,0)
# ifdef MASKING
                v(Istr-1,j,k,nnew)=v(Istr-1,j,k,nnew)*                  &
     &                             vmask(Istr-1,j)
# endif

              END DO
            END IF
          END IF
!
          IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
            IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
              DO k=1,N(ng)
                v(Iend+1,j,k,nnew)=v(Iend+1,j,k,nnew)-CF(Iend+1,0)
# ifdef MASKING
                v(Iend+1,j,k,nnew)=v(Iend+1,j,k,nnew)*                  &
     &                             vmask(Iend+1,j)
# endif

              END DO
            END IF
          END IF
!
          IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
            IF (j.eq.1) THEN                         ! southern boundary
              DO k=1,N(ng)                           ! J-loop pipelined
                DO i=Istr,Iend
                  v(i,j,k,nnew)=v(i,j,k,nnew)-CF(i,0)
# ifdef MASKING
                  v(i,j,k,nnew)=v(i,j,k,nnew)*                          &
     &                          vmask(i,j)
# endif

                END DO
              END DO
            END IF
          END IF
!
          IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
            IF (j.eq.Mm(ng)+1) THEN                  ! northern boundary
              DO k=1,N(ng)                           ! J-loop pipelined
                DO i=Istr,Iend
                  v(i,j,k,nnew)=v(i,j,k,nnew)-CF(i,0)
# ifdef MASKING
                  v(i,j,k,nnew)=v(i,j,k,nnew)*                          &
     &                          vmask(i,j)
# endif

                END DO
              END DO
            END IF
          END IF
!
!  Compute correct mass flux, Hz*v/m.
!
          DO k=N(ng),1,-1
            DO i=IstrT,IendT
              Hvom(i,j,k)=0.5_r8*(Hvom(i,j,k)+v(i,j,k,nnew)*DC(i,k))

              FC(i,0)=FC(i,0)+Hvom(i,j,k)

            END DO
          END DO
          DO i=IstrT,IendT
            FC(i,0)=DC(i,0)*(FC(i,0)-DV_avg2(i,j))      ! recursive
          END DO
          DO k=1,N(ng)
            DO i=IstrT,IendT
              Hvom(i,j,k)=Hvom(i,j,k)-DC(i,k)*FC(i,0)
            END DO
          END DO
        END IF
      END DO
