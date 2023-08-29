                                                       !
# This subroutine evaluates right-hand-side terms for 3D momentum and tracers equations.


def rhs3d():

    # Initialize computations for new time step of the 3D primitive
    # variables.
    pre_step3d()

    # Compute baroclinic pressure gradient.
    prsgrd()



    # Compute right-hand-side terms for the 3D momentum equations.
    # -----------------------------------------------------------------------




    # ***** ifdef BODYFORCE

    # Apply surface stress as a bodyforce: determine the thickness (m) of the surface layer; then add in surface stress as a bodyfoce.

    # ***** FOR JstrV - 1, Jend,  IstrU - 1, Iend
    wrk = Hz[levsfrc:N].sum(dim = 0)



      DO j=Jstr,Jend
        DO i=IstrU,Iend
          cff=0.25_r8*(pm(i-1,j)+pm(i,j))*                              &
     &                (pn(i-1,j)+pn(i,j))
          cff1=1.0_r8/(cff*(wrk(i-1,j)+wrk(i,j)))
          Uwrk(i,j)=sustr(i,j)*cff1
        END DO
      END DO

      DO j=JstrV,Jend
        DO i=Istr,Iend
          cff=0.25*(pm(i,j-1)+pm(i,j))*                                 &
     &             (pn(i,j-1)+pn(i,j))
          cff1=1.0_r8/(cff*(wrk(i,j-1)+wrk(i,j)))
          Vwrk(i,j)=svstr(i,j)*cff1
        END DO
      END DO


      DO k=levsfrc(ng),N(ng)
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            cff=Uwrk(i,j)*(Hz(i  ,j,k)+                                 &
     &                     Hz(i-1,j,k))
            ru(i,j,k,nrhs)=ru(i,j,k,nrhs)+cff

          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend
            cff=Vwrk(i,j)*(Hz(i,j  ,k)+                                 &
     &                     Hz(i,j-1,k))
            rv(i,j,k,nrhs)=rv(i,j,k,nrhs)+cff

          END DO
        END DO
      END DO

    # Apply bottom stress as a bodyforce: determine the thickness (m) of the bottom layer; then add in bottom stress as a bodyfoce.
!
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
          wrk(i,j)=0.0_r8
        END DO
      END DO
      DO k=1,levbfrc(ng)
        DO j=JstrV-1,Jend
          DO i=IstrU-1,Iend
            wrk(i,j)=wrk(i,j)+Hz(i,j,k)
          END DO
        END DO
      END DO
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          cff=0.25_r8*(pm (i-1,j)+pm (i,j))*                            &
     &                (pn (i-1,j)+pn (i,j))
          cff1=1.0_r8/(cff*(wrk(i-1,j)+wrk(i,j)))
          Uwrk(i,j)=bustr(i,j)*cff1
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
          cff=0.25_r8*(pm (i,j-1)+pm (i,j))*                            &
     &                (pn (i,j-1)+pn (i,j))
          cff1=1.0_r8/(cff*(wrk(i,j-1)+wrk(i,j)))
          Vwrk(i,j)=bvstr(i,j)*cff1
        END DO
      END DO
      DO k=1,levbfrc(ng)
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            cff=Uwrk(i,j)*(Hz(i  ,j,k)+                                 &
     &                     Hz(i-1,j,k))
            ru(i,j,k,nrhs)=ru(i,j,k,nrhs)-cff

          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend
            cff=Vwrk(i,j)*(Hz(i,j  ,k)+                                 &
     &                     Hz(i,j-1,k))
            rv(i,j,k,nrhs)=rv(i,j,k,nrhs)-cff


          END DO
        END DO
      END DO
# endif



      K_LOOP : DO k=1,N(ng)

# ifdef UV_COR


    # Add in Coriolis terms.
        DO j=JstrV-1,Jend
          DO i=IstrU-1,Iend
            cff=0.5_r8*Hz(i,j,k)*fomn(i,j)
            UFx(i,j)=cff*(v(i,j  ,k,nrhs)+                              &
     &                    v(i,j+1,k,nrhs))
            VFe(i,j)=cff*(u(i  ,j,k,nrhs)+                              &
     &                    u(i+1,j,k,nrhs))
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            cff1=0.5_r8*(UFx(i,j)+UFx(i-1,j))
            ru(i,j,k,nrhs)=ru(i,j,k,nrhs)+cff1

          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend
            cff1=0.5_r8*(VFe(i,j)+VFe(i,j-1))
            rv(i,j,k,nrhs)=rv(i,j,k,nrhs)-cff1

          END DO
        END DO
!

# endif



# if defined CURVGRID && defined UV_ADV
!
!-----------------------------------------------------------------------
!  Add in curvilinear transformation terms.
!-----------------------------------------------------------------------
!
        DO j=JstrV-1,Jend
          DO i=IstrU-1,Iend
            cff1=0.5_r8*(v(i,j  ,k,nrhs)+                               &

     &                   v(i,j+1,k,nrhs))
            cff2=0.5_r8*(u(i  ,j,k,nrhs)+                               &

     &                   u(i+1,j,k,nrhs))
            cff3=cff1*dndx(i,j)
            cff4=cff2*dmde(i,j)

            cff=Hz(i,j,k)*(cff3-cff4)
            UFx(i,j)=cff*cff1
            VFe(i,j)=cff*cff2

          END DO
        END DO
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            cff1=0.5_r8*(UFx(i,j)+UFx(i-1,j))
            ru(i,j,k,nrhs)=ru(i,j,k,nrhs)+cff1

          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend
            cff1=0.5_r8*(VFe(i,j)+VFe(i,j-1))
            rv(i,j,k,nrhs)=rv(i,j,k,nrhs)-cff1

          END DO
        END DO
# endif



!  Add in nudging of 3D momentum climatology.
!-----------------------------------------------------------------------
!
        IF (LnudgeM3CLM(ng)) THEN
          DO j=Jstr,Jend
            DO i=IstrU,Iend
              cff=0.25_r8*(CLIMA(ng)%M3nudgcof(i-1,j,k)+                &
     &                     CLIMA(ng)%M3nudgcof(i  ,j,k))*               &
     &            om_u(i,j)*on_u(i,j)
              ru(i,j,k,nrhs)=ru(i,j,k,nrhs)+                            &
     &                       cff*(Hz(i-1,j,k)+Hz(i,j,k))*               &
     &                       (CLIMA(ng)%uclm(i,j,k)-                    &
     &                        u(i,j,k,nrhs))
            END DO
          END DO
          DO j=JstrV,Jend
            DO i=Istr,Iend
              cff=0.25_r8*(CLIMA(ng)%M3nudgcof(i,j-1,k)+                &
     &                     CLIMA(ng)%M3nudgcof(i,j  ,k))*               &
     &            om_v(i,j)*on_v(i,j)
              rv(i,j,k,nrhs)=rv(i,j,k,nrhs)+                            &
     &                       cff*(Hz(i,j-1,k)+Hz(i,j,k))*               &
     &                       (CLIMA(ng)%vclm(i,j,k)-                    &
     &                        v(i,j,k,nrhs))
            END DO
          END DO
        END IF

# ifdef UV_ADV
!
!-----------------------------------------------------------------------
!  Add in horizontal advection of momentum.
!-----------------------------------------------------------------------
!
!  Compute diagonal [UFx,VFe] and off-diagonal [UFe,VFx] components
!  of tensor of momentum flux due to horizontal advection.
!
#  ifdef UV_C2ADVECTION
!
!  Second-order, centered differences advection.
!
        DO j=Jstr,Jend
          DO i=IstrU-1,Iend
            UFx(i,j)=0.25_r8*(u(i  ,j,k,nrhs)+                          &


     &                        u(i+1,j,k,nrhs))*                         &
     &                       (Huon(i  ,j,k)+                            &
     &                        Huon(i+1,j,k))
          END DO
        END DO
        DO j=Jstr,Jend+1
          DO i=IstrU,Iend
            UFe(i,j)=0.25_r8*(u(i,j-1,k,nrhs)+                          &


     &                        u(i,j  ,k,nrhs))*                         &
     &                       (Hvom(i-1,j,k)+                            &
     &                        Hvom(i  ,j,k))
          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend+1
            VFx(i,j)=0.25_r8*(v(i-1,j,k,nrhs)+                          &


     &                        v(i  ,j,k,nrhs))*                         &
     &                       (Huon(i,j-1,k)+                            &
     &                        Huon(i,j  ,k))
          END DO
        END DO
        DO j=JstrV-1,Jend
          DO i=Istr,Iend
            VFe(i,j)=0.25_r8*(v(i,j  ,k,nrhs)+                          &


     &                        v(i,j+1,k,nrhs))*                         &
     &                       (Hvom(i,j  ,k)+                            &
     &                        Hvom(i,j+1,k))
          END DO
        END DO
#  else
        DO j=Jstr,Jend
          DO i=IstrUm1,Iendp1
            uxx(i,j)=u(i-1,j,k,nrhs)-2.0_r8*u(i,j,k,nrhs)+              &


     &               u(i+1,j,k,nrhs)
            Huxx(i,j)=Huon(i-1,j,k)-2.0_r8*Huon(i,j,k)+Huon(i+1,j,k)
          END DO
        END DO
        IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            DO j=Jstr,Jend
              uxx (Istr,j)=uxx (Istr+1,j)
              Huxx(Istr,j)=Huxx(Istr+1,j)
            END DO
          END IF
        END IF
        IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO j=Jstr,Jend
              uxx (Iend+1,j)=uxx (Iend,j)
              Huxx(Iend+1,j)=Huxx(Iend,j)
            END DO
          END IF
        END IF
#   ifdef UV_C4ADVECTION
!
!  Fourth-order, centered differences u-momentum horizontal advection.
!
        cff=1.0_r8/6.0_r8
        DO j=Jstr,Jend
          DO i=IstrU-1,Iend
            UFx(i,j)=0.25_r8*(u(i  ,j,k,nrhs)+                          &


     &                        u(i+1,j,k,nrhs)-                          &
     &                        cff*(uxx (i  ,j)+                         &
     &                             uxx (i+1,j)))*                       &
     &                       (Huon(i  ,j,k)+                            &
     &                        Huon(i+1,j,k)-                            &
     &                        cff*(Huxx(i  ,j)+                         &
     &                             Huxx(i+1,j)))
          END DO
        END DO
#   else
!
!  Third-order, upstream bias u-momentum advection with velocity
!  dependent hyperdiffusion.
!
        DO j=Jstr,Jend
          DO i=IstrU-1,Iend
            cff1=u(i  ,j,k,nrhs)+                                       &


     &           u(i+1,j,k,nrhs)
            IF (cff1.gt.0.0_r8) THEN
              cff=uxx(i,j)
            ELSE
              cff=uxx(i+1,j)
            END IF
            UFx(i,j)=0.25_r8*(cff1+Gadv*cff)*                           &
     &               (Huon(i  ,j,k)+                                    &
     &                Huon(i+1,j,k)+                                    &
     &                Gadv*0.5_r8*(Huxx(i  ,j)+                         &
     &                             Huxx(i+1,j)))
          END DO
        END DO
#   endif
        DO j=Jstrm1,Jendp1
          DO i=IstrU,Iend
            uee(i,j)=u(i,j-1,k,nrhs)-2.0_r8*u(i,j,k,nrhs)+              &


     &               u(i,j+1,k,nrhs)
          END DO
        END DO
        IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO i=IstrU,Iend
              uee(i,Jstr-1)=uee(i,Jstr)
            END DO
          END IF
        END IF
        IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO i=IstrU,Iend
              uee(i,Jend+1)=uee(i,Jend)
            END DO
          END IF
        END IF
        DO j=Jstr,Jend+1
          DO i=IstrU-1,Iend
           Hvxx(i,j)=Hvom(i-1,j,k)-2.0_r8*Hvom(i,j,k)+Hvom(i+1,j,k)
          END DO
        END DO
#   ifdef UV_C4ADVECTION
        cff=1.0_r8/6.0_r8
        DO j=Jstr,Jend+1
          DO i=IstrU,Iend
            UFe(i,j)=0.25_r8*(u(i,j  ,k,nrhs)+                          &

     &                        u(i,j-1,k,nrhs)-                          &
     &                        cff*(uee (i,j  )+                         &
     &                             uee (i,j-1)))*                       &
     &                       (Hvom(i  ,j,k)+                            &
     &                        Hvom(i-1,j,k)-                            &
     &                        cff*(Hvxx(i  ,j)+                         &
     &                             Hvxx(i-1,j)))
          END DO
        END DO
#   else
        DO j=Jstr,Jend+1
          DO i=IstrU,Iend
            cff1=u(i,j  ,k,nrhs)+                                       &

     &           u(i,j-1,k,nrhs)
            cff2=Hvom(i,j,k)+Hvom(i-1,j,k)
            IF (cff2.gt.0.0_r8) THEN
              cff=uee(i,j-1)
            ELSE
              cff=uee(i,j)
            END IF
            UFe(i,j)=0.25_r8*(cff1+Gadv*cff)*                           &
     &               (cff2+Gadv*0.5_r8*(Hvxx(i  ,j)+                    &
     &                                  Hvxx(i-1,j)))
          END DO
        END DO
#   endif
        DO j=JstrV,Jend
          DO i=Istrm1,Iendp1
            vxx(i,j)=v(i-1,j,k,nrhs)-2.0_r8*v(i,j,k,nrhs)+              &

     &               v(i+1,j,k,nrhs)
          END DO
        END DO
        IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            DO j=JstrV,Jend
              vxx(Istr-1,j)=vxx(Istr,j)
            END DO
          END IF
        END IF
        IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO j=JstrV,Jend
              vxx(Iend+1,j)=vxx(Iend,j)
            END DO
          END IF
        END IF
        DO j=JstrV-1,Jend
          DO i=Istr,Iend+1
           Huee(i,j)=Huon(i,j-1,k)-2.0_r8*Huon(i,j,k)+Huon(i,j+1,k)
          END DO
        END DO
#   ifdef UV_C4ADVECTION

    # Fourth-order, centered differences v-momentum horizontal advection.

        cff=1.0_r8/6.0_r8
        DO j=JstrV,Jend
          DO i=Istr,Iend+1
            VFx(i,j)=0.25_r8*(v(i  ,j,k,nrhs)+                          &

     &                        v(i-1,j,k,nrhs)-                          &
     &                        cff*(vxx (i  ,j)+                         &
     &                             vxx (i-1,j)))*                       &
     &                       (Huon(i,j  ,k)+                            &
     &                        Huon(i,j-1,k)-                            &
     &                        cff*(Huee(i,j  )+                         &
     &                             Huee(i,j-1)))
          END DO
        END DO
#   else
!
!  Third-order, upstream bias v-momentum advection with velocity
!  dependent hyperdiffusion.
!
        DO j=JstrV,Jend
          DO i=Istr,Iend+1
            cff1=v(i  ,j,k,nrhs)+                                       &

     &           v(i-1,j,k,nrhs)
            cff2=Huon(i,j,k)+Huon(i,j-1,k)
            IF (cff2.gt.0.0_r8) THEN
              cff=vxx(i-1,j)
            ELSE
              cff=vxx(i,j)
            END IF
            VFx(i,j)=0.25_r8*(cff1+Gadv*cff)*                           &
     &               (cff2+Gadv*0.5_r8*(Huee(i,j  )+                    &
     &                                  Huee(i,j-1)))
          END DO
        END DO
#   endif
        DO j=JstrVm1,Jendp1
          DO i=Istr,Iend
            vee(i,j)=v(i,j-1,k,nrhs)-2.0_r8*v(i,j,k,nrhs)+              &

     &               v(i,j+1,k,nrhs)
            Hvee(i,j)=Hvom(i,j-1,k)-2.0_r8*Hvom(i,j,k)+Hvom(i,j+1,k)
          END DO
        END DO
        IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO i=Istr,Iend
              vee (i,Jstr)=vee (i,Jstr+1)
              Hvee(i,Jstr)=Hvee(i,Jstr+1)
            END DO
          END IF
        END IF
        IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO i=Istr,Iend
              vee (i,Jend+1)=vee (i,Jend)
              Hvee(i,Jend+1)=Hvee(i,Jend)
            END DO
          END IF
        END IF
#   ifdef UV_C4ADVECTION
        cff=1.0_r8/6.0_r8
        DO j=JstrV-1,Jend
          DO i=Istr,Iend
            VFe(i,j)=0.25_r8*(v(i,j  ,k,nrhs)+                          &

     &                        v(i,j+1,k,nrhs)-                          &
     &                        cff*(vee (i,j  )+                         &
     &                             vee (i,j+1)))*                       &
     &                       (Hvom(i,j  ,k)+                            &
     &                        Hvom(i,j+1,k)-                            &
     &                        cff*(Hvee(i,j  )+                         &
     &                             Hvee(i,j+1)))
          END DO
        END DO
#   else
        DO j=JstrV-1,Jend
          DO i=Istr,Iend
            cff1=v(i,j  ,k,nrhs)+                                       &

     &           v(i,j+1,k,nrhs)
            IF (cff1.gt.0.0_r8) THEN
              cff=vee(i,j)
            ELSE
              cff=vee(i,j+1)
            END IF
            VFe(i,j)=0.25_r8*(cff1+Gadv*cff)*                           &
     &               (Hvom(i,j  ,k)+                                    &
     &                Hvom(i,j+1,k)+                                    &
     &                Gadv*0.5_r8*(Hvee(i,j  )+                         &
     &                             Hvee(i,j+1)))
          END DO
        END DO
#   endif
#  endif
!
!  Add in horizontal advection.
!
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            cff1=UFx(i,j)-UFx(i-1,j)
            cff2=UFe(i,j+1)-UFe(i,j)
            cff=cff1+cff2
            ru(i,j,k,nrhs)=ru(i,j,k,nrhs)-cff

          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend
            cff1=VFx(i+1,j)-VFx(i,j)
            cff2=VFe(i,j)-VFe(i,j-1)
            cff=cff1+cff2
            rv(i,j,k,nrhs)=rv(i,j,k,nrhs)-cff

          END DO
        END DO
# endif


      END DO K_LOOP
!
      J_LOOP : DO j=Jstr,Jend
# ifdef UV_ADV
!
!-----------------------------------------------------------------------
!  Add in vertical advection.
!-----------------------------------------------------------------------
!
#  ifdef UV_SADVECTION
!
!  Construct conservative parabolic splines for the vertical
!  derivatives "CF" of u-momentum.
!
        cff1=9.0_r8/16.0_r8
        cff2=1.0_r8/16.0_r8
        DO k=1,N(ng)
          DO i=IstrU,Iend
            DC(i,k)=cff1*(Hz(i  ,j,k)+                                  &
     &                    Hz(i-1,j,k))-                                 &
     &              cff2*(Hz(i+1,j,k)+                                  &
     &                    Hz(i-2,j,k))
          END DO
        END DO
        DO i=IstrU,Iend
          FC(i,0)=0.0_r8
          CF(i,0)=0.0_r8
        END DO
        DO k=1,N(ng)-1
          DO i=IstrU,Iend
            cff=1.0_r8/(2.0_r8*DC(i,k+1)+DC(i,k)*(2.0_r8-FC(i,k-1)))
            FC(i,k)=cff*DC(i,k+1)
            CF(i,k)=cff*(6.0_r8*(u(i,j,k+1,nrhs)-                       &

     &                           u(i,j,k  ,nrhs))-                      &
     &                   DC(i,k)*CF(i,k-1))
          END DO
        END DO
        DO i=IstrU,Iend
          CF(i,N(ng))=0.0_r8
        END DO
        DO k=N(ng)-1,1,-1
          DO i=IstrU,Iend
            CF(i,k)=CF(i,k)-FC(i,k)*CF(i,k+1)
          END DO
        END DO
!
! Compute spline-interpolated, vertical advective u-momentum flux.
!
        cff3=1.0_r8/3.0_r8
        cff4=1.0_r8/6.0_r8
        DO k=1,N(ng)-1
          DO i=IstrU,Iend
            FC(i,k)=(cff1*(W(i  ,j,k)+                                  &
     &                     W(i-1,j,k))-                                 &
     &               cff2*(W(i+1,j,k)+                                  &
     &                     W(i-2,j,k)))*                                &
     &              (u(i,j,k,nrhs)+                                     &

     &               DC(i,k)*(cff3*CF(i,k  )+                           &
     &                        cff4*CF(i,k-1)))
          END DO
        END DO
        DO i=IstrU,Iend
          FC(i,N(ng))=0.0_r8
          FC(i,0)=0.0_r8
        END DO
#  elif defined UV_C2ADVECTION
        DO k=1,N(ng)-1
          DO i=IstrU,Iend
            FC(i,k)=0.25_r8*(u(i,j,k  ,nrhs)+                           &

     &                       u(i,j,k+1,nrhs))*                          &
     &              (W(i  ,j,k)+                                        &
     &               W(i-1,j,k))
          END DO
        END DO
        DO i=IstrU,Iend
          FC(i,0)=0.0_r8
          FC(i,N(ng))=0.0_r8
        END DO
#  elif defined UV_C4ADVECTION
        cff1=9.0_r8/32.0_r8
        cff2=1.0_r8/32.0_r8
        DO k=2,N(ng)-2
          DO i=IstrU,Iend
            FC(i,k)=(cff1*(u(i,j,k  ,nrhs)+                             &

     &                     u(i,j,k+1,nrhs))-                            &
     &               cff2*(u(i,j,k-1,nrhs)+                             &

     &                     u(i,j,k+2,nrhs)))*                           &
     &              (W(i  ,j,k)+                                        &
     &               W(i-1,j,k))
          END DO
        END DO
        DO i=IstrU,Iend
          FC(i,N(ng))=0.0_r8
          FC(i,N(ng)-1)=(cff1*(u(i,j,N(ng)-1,nrhs)+                     &

     &                         u(i,j,N(ng)  ,nrhs))-                    &
     &                   cff2*(u(i,j,N(ng)-2,nrhs)+                     &

     &                         u(i,j,N(ng)  ,nrhs)))*                   &
     &                  (W(i  ,j,N(ng)-1)+                              &
     &                   W(i-1,j,N(ng)-1))
          FC(i,1)=(cff1*(u(i,j,1,nrhs)+                                 &

     &                   u(i,j,2,nrhs))-                                &
     &             cff2*(u(i,j,1,nrhs)+                                 &

     &                   u(i,j,3,nrhs)))*                               &
     &            (W(i  ,j,1)+                                          &
     &             W(i-1,j,1))
          FC(i,0)=0.0_r8
        END DO
#  else
        cff1=9.0_r8/16.0_r8
        cff2=1.0_r8/16.0_r8
        DO k=2,N(ng)-2
          DO i=IstrU,Iend
            FC(i,k)=(cff1*(u(i,j,k  ,nrhs)+                             &

     &                     u(i,j,k+1,nrhs))-                            &
     &               cff2*(u(i,j,k-1,nrhs)+                             &

     &                     u(i,j,k+2,nrhs)))*                           &
     &              (cff1*(W(i  ,j,k)+                                  &
     &                     W(i-1,j,k))-                                 &
     &               cff2*(W(i+1,j,k)+                                  &
     &                     W(i-2,j,k)))
          END DO
        END DO
        DO i=IstrU,Iend
          FC(i,N(ng))=0.0_r8
          FC(i,N(ng)-1)=(cff1*(u(i,j,N(ng)-1,nrhs)+                     &

     &                         u(i,j,N(ng)  ,nrhs))-                    &
     &                   cff2*(u(i,j,N(ng)-2,nrhs)+                     &

     &                         u(i,j,N(ng)  ,nrhs)))*                   &
     &                  (cff1*(W(i  ,j,N(ng)-1)+                        &
     &                         W(i-1,j,N(ng)-1))-                       &
     &                   cff2*(W(i+1,j,N(ng)-1)+                        &
     &                         W(i-2,j,N(ng)-1)))
          FC(i,1)=(cff1*(u(i,j,1,nrhs)+                                 &

     &                   u(i,j,2,nrhs))-                                &
     &             cff2*(u(i,j,1,nrhs)+                                 &

     &                   u(i,j,3,nrhs)))*                               &
     &            (cff1*(W(i  ,j,1)+                                    &
     &                   W(i-1,j,1))-                                   &
     &             cff2*(W(i+1,j,1)+                                    &
     &                   W(i-2,j,1)))
          FC(i,0)=0.0_r8
        END DO
#  endif
        DO k=1,N(ng)
          DO i=IstrU,Iend
            cff=FC(i,k)-FC(i,k-1)
            ru(i,j,k,nrhs)=ru(i,j,k,nrhs)-cff

          END DO
        END DO
        IF (j.ge.JstrV) THEN
#  ifdef UV_SADVECTION
!
!  Construct conservative parabolic splines for the vertical
!  derivatives "CF" of v-momentum.
!
          cff1=9.0_r8/16.0_r8
          cff2=1.0_r8/16.0_r8
          DO k=1,N(ng)
            DO i=Istr,Iend
              DC(i,k)=(cff1*(Hz(i,j  ,k)+                               &
     &                       Hz(i,j-1,k))-                              &
     &                 cff2*(Hz(i,j+1,k)+                               &
     &                       Hz(i,j-2,k)))
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,0)=0.0_r8
            CF(i,0)=0.0_r8
          END DO
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              cff=1.0_r8/(2.0_r8*DC(i,k+1)+DC(i,k)*(2.0_r8-FC(i,k-1)))
              FC(i,k)=cff*DC(i,k+1)
              CF(i,k)=cff*(6.0_r8*(v(i,j,k+1,nrhs)-                     &

     &                             v(i,j,k  ,nrhs))-                    &
     &                     DC(i,k)*CF(i,k-1))
            END DO
          END DO
          DO i=Istr,Iend
            CF(i,N(ng))=0.0_r8
          END DO
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
              CF(i,k)=CF(i,k)-FC(i,k)*CF(i,k+1)
            END DO
          END DO
!
! Compute spline-interpolated, vertical advective v-momentum flux.
!
          cff3=1.0_r8/3.0_r8
          cff4=1.0_r8/6.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FC(i,k)=(cff1*(W(i,j  ,k)+                                &
     &                       W(i,j-1,k))-                               &
     &                 cff2*(W(i,j+1,k)+                                &
     &                       W(i,j-2,k)))*                              &
     &                (v(i,j,k,nrhs)+                                   &

     &                 DC(i,k)*(cff3*CF(i,k  )+                         &
     &                          cff4*CF(i,k-1)))
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,N(ng))=0.0_r8
            FC(i,0)=0.0_r8
          END DO
#  elif defined UV_C2ADVECTION
!
!  Second-order, centered differences vertical advection.
!
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FC(i,k)=0.25_r8*(v(i,j,k  ,nrhs)+                         &

     &                         v(i,j,k+1,nrhs))*                        &
     &                (W(i,j  ,k)+                                      &
     &                 W(i,j-1,k))
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,0)=0.0_r8
            FC(i,N(ng))=0.0_r8
          END DO
#  elif defined UV_C4ADVECTION
!
!  Forth-order, centered differences vertical advection.
!
          cff1=9.0_r8/32.0_r8
          cff2=1.0_r8/32.0_r8
          DO k=2,N(ng)-2
            DO i=Istr,Iend
              FC(i,k)=(cff1*(v(i,j,k  ,nrhs)+                           &

     &                       v(i,j,k+1,nrhs))-                          &
     &                 cff2*(v(i,j,k-1,nrhs)+                           &

     &                       v(i,j,k+2,nrhs)))*                         &
     &                (W(i,j  ,k)+                                      &
     &                 W(i,j-1,k))
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,N(ng))=0.0_r8
            FC(i,N(ng)-1)=(cff1*(v(i,j,N(ng)-1,nrhs)+                   &

     &                           v(i,j,N(ng)  ,nrhs))-                  &
     &                     cff2*(v(i,j,N(ng)-2,nrhs)+                   &

     &                           v(i,j,N(ng)  ,nrhs)))*                 &
     &                    (W(i,j  ,N(ng)-1)+                            &
     &                     W(i,j-1,N(ng)-1))
            FC(i,1)=(cff1*(v(i,j,1,nrhs)+                               &

     &                     v(i,j,2,nrhs))-                              &
     &               cff2*(v(i,j,1,nrhs)+                               &

     &                     v(i,j,3,nrhs)))*                             &
     &              (W(i,j  ,1)+                                        &
     &               W(i,j-1,1))
            FC(i,0)=0.0_r8
          END DO
#  else
          cff1=9.0_r8/16.0_r8
          cff2=1.0_r8/16.0_r8
          DO k=2,N(ng)-2
            DO i=Istr,Iend
              FC(i,k)=(cff1*(v(i,j,k  ,nrhs)+                           &

     &                       v(i,j,k+1,nrhs))-                          &
     &                 cff2*(v(i,j,k-1,nrhs)+                           &

     &                       v(i,j,k+2,nrhs)))*                         &
     &                (cff1*(W(i,j  ,k)+                                &
     &                       W(i,j-1,k))-                               &
     &                 cff2*(W(i,j+1,k)+                                &
     &                       W(i,j-2,k)))
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,N(ng))=0.0_r8
            FC(i,N(ng)-1)=(cff1*(v(i,j,N(ng)-1,nrhs)+                   &

     &                           v(i,j,N(ng)  ,nrhs))-                  &
     &                     cff2*(v(i,j,N(ng)-2,nrhs)+                   &

     &                           v(i,j,N(ng)  ,nrhs)))*                 &
     &                    (cff1*(W(i,j  ,N(ng)-1)+                      &
     &                           W(i,j-1,N(ng)-1))-                     &
     &                     cff2*(W(i,j+1,N(ng)-1)+                      &
     &                           W(i,j-2,N(ng)-1)))
            FC(i,1)=(cff1*(v(i,j,1,nrhs)+                               &

     &                     v(i,j,2,nrhs))-                              &
     &               cff2*(v(i,j,1,nrhs)+                               &

     &                     v(i,j,3,nrhs)))*                             &
     &              (cff1*(W(i,j  ,1)+                                  &
     &                     W(i,j-1,1))-                                 &
     &               cff2*(W(i,j+1,1)+                                  &
     &                     W(i,j-2,1)))
            FC(i,0)=0.0_r8
          END DO
#  endif
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff=FC(i,k)-FC(i,k-1)
              rv(i,j,k,nrhs)=rv(i,j,k,nrhs)-cff

            END DO
          END DO
        END IF

!
!-----------------------------------------------------------------------
!  Compute forcing term for the 2D momentum equations.
!-----------------------------------------------------------------------
!
!  Vertically integrate baroclinic right-hand-side terms. If not
!  body force stresses, add in the difference between surface and
!  bottom stresses.
!
        DO i=IstrU,Iend

          rufrc(i,j)=ru(i,j,1,nrhs)

        END DO
        DO k=2,N(ng)
          DO i=IstrU,Iend

            rufrc(i,j) = rufrc(i,j) + ru(i,j,k,nrhs)

          END DO
        END DO

# ifndef BODYFORCE
        DO i=IstrU,Iend
          cff=om_u(i,j)*on_u(i,j)
          cff1= sustr(i,j)*cff
          cff2=-bustr(i,j)*cff
          rufrc(i,j)=rufrc(i,j)+cff1+cff2

        IF (j.ge.JstrV) THEN
          DO i=Istr,Iend

            rvfrc(i,j)=rv(i,j,1,nrhs)

          END DO
          DO k=2,N(ng)
            DO i=Istr,Iend

            END DO
          END DO

# ifndef BODYFORCE
          DO i=Istr,Iend
            cff=om_v(i,j)*on_v(i,j)
            cff1= svstr(i,j)*cff
            cff2=-bvstr(i,j)*cff
            rvfrc(i,j)=rvfrc(i,j)+cff1+cff2

        END IF
      END DO J_LOOP


#endif
