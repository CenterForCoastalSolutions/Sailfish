                                                       !
# This subroutine evaluates right-hand-side terms for 3D momentum and tracers equations.


def rhs3d():

    # Initialize computations for new time step of the 3D primitive variables.
    pre_step3d()

    # Compute baroclinic pressure gradient.
    prsgrd()


    # Compute right-hand-side terms for the 3D momentum equations.
    # -----------------------------------------------------------------------


    K_LOOP : DO k=1,N(ng)




    # Add in Coriolis terms.

    addCoriolis()
    # ifdef UV_COR
        # for j=JstrV-1,Jend:
        #   for i=IstrU-1,Iend:
        #     cff=Hz(i,j,k)*fomn(i,j)
        #     UFξ(i,j)=cff*0.5*(v(i,j,k,nrhs) + v(i,j+1,k,nrhs))
        #     VFη(i,j)=cff*0.5*(u(i,j,k,nrhs) + u(i+1,j,k,nrhs))

        # for j=Jstr,Jend:
        #   for i=IstrU,Iend:
        #     cff1=0.5*(UFx(i,j) + UFx(i-1,j))
        #     ru(i,j,k,nrhs) = ru(i,j,k,nrhs) + cff1
        #
        #
        # for j=JstrV,Jend:
        #   for i=Istr,Iend:
        #     cff1=0.5*(VFe(i,j) + VFe(i,j-1))
        #     rv(i,j,k,nrhs) = rv(i,j,k,nrhs) - cff1

    # endif



    pass
    # if defined CURVGRID && defined UV_ADV
    # These are the terms (u*d_xi(1/n) - v*d_eta(1/m))*Hz*u  and  (v*d_xi(1/n) - u*d_eta(1/m))*Hz*v
    pass
    # Add in curvilinear transformation terms.
    # -----------------------------------------------------------------------



        # DO j=JstrV-1,Jend
        #   DO i=IstrU-1,Iend
        #     cff1=0.5*(v(i,j  ,k,nrhs) + v(i,j+1,k,nrhs))
        #     cff2=0.5*(u(i  ,j,k,nrhs) + u(i+1,j,k,nrhs))
        #     cff3=cff1*dndx(i,j)
        #     cff4=cff2*dmde(i,j)
        #
        #     cff=Hz(i,j,k)*(cff3-cff4)
        #     UFx(i,j)=cff*cff1
        #     VFe(i,j)=cff*cff2
        #   END DO
        # END DO
        #
        # DO j=Jstr,Jend
        #   DO i=IstrU,Iend
        #     cff1=0.5*(UFx(i,j)+UFx(i-1,j))
        #     ru(i,j,k,nrhs) = ru(i,j,k,nrhs) + cff1
        #   END DO
        # END DO
        #
        # DO j=JstrV,Jend
        #   DO i=Istr,Iend
        #     cff1=0.5*(VFe(i,j)+VFe(i,j-1))
        #     rv(i,j,k,nrhs) = rv(i,j,k,nrhs) - cff1
        #   END DO
        # END DO
# endif

    pass
    # Add in nudging towards 3D momentum climatology.
    # -----------------------------------------------------------------------
        IF (LnudgeM3CLM(ng)) THEN
          DO j=Jstr,Jend
            DO i=IstrU,Iend
              cff=0.25*(CLIMA(ng)%M3nudgcof(i-1,j,k)  CLIMA(ng)%M3nudgcof(i  ,j,k))*om_u(i,j)*on_u(i,j)
              ru(i,j,k,nrhs)=ru(i,j,k,nrhs) + cff*(Hz(i-1,j,k)+Hz(i,j,k))*(CLIMA(ng)%uclm(i,j,k) - u(i,j,k,nrhs))
            END DO
          END DO

          DO j=JstrV,Jend
            DO i=Istr,Iend
              cff=0.25*(CLIMA(ng)%M3nudgcof(i,j-1,k) + CLIMA(ng)%M3nudgcof(i,j  ,k))*om_v(i,j)*on_v(i,j)
              rv(i,j,k,nrhs)=rv(i,j,k,nrhs) + cff*(Hz(i,j-1,k)+Hz(i,j,k))*(CLIMA(ng)%vclm(i,j,k) - v(i,j,k,nrhs))
            END DO
          END DO
        END IF



# ifdef UV_ADV


# Add in horizontal advection of momentum.
# -----------------------------------------------------------------------

# Compute diagonal [UFx,VFe] and off-diagonal [UFe,VFx] components
# of tensor of momentum flux due to horizontal advection.

#  ifdef UV_C2ADVECTION
#
#   Second-order, centered differences advection.
#
    UFx = UtoR(u[nrhs,k,:,:])*UtoR(Huon[k,:,:])
    UFe = UtoR(u[nrhs,k,:,:])*UtoR(Huon[k,:,:])
        DO j=Jstr,Jend
          DO i=IstrU-1,Iend
            UFx(i,j) = 0.25*(u(i  ,j,k,nrhs) + u(i+1,j,k,nrhs))*(Huon(i  ,j,k) + Huon(i+1,j,k))
          END DO
        END DO

        DO j=Jstr,Jend+1
          DO i=IstrU,Iend
            UFe(i,j) = 0.25*(u(i,j-1,k,nrhs) + u(i,j  ,k,nrhs))*(Hvom(i-1,j,k) + Hvom(i  ,j,k))
          END DO
        END DO

        DO j=JstrV,Jend
          DO i=Istr,Iend+1
            VFx(i,j) = 0.25*(v(i-1,j,k,nrhs) + v(i  ,j,k,nrhs))*(Huon(i,j-1,k) + Huon(i,j  ,k))
          END DO
        END DO

        DO j=JstrV-1,Jend
          DO i=Istr,Iend
            VFe(i,j) = 0.25*(v(i,j  ,k,nrhs) + v(i,j+1,k,nrhs))*(Hvom(i,j  ,k) + Hvom(i,j+1,k))
          END DO
        END DO

#  else
        DO j=Jstr,Jend
          DO i=IstrUm1,Iendp1
            uxx(i,j)  = u(i-1,j,k,nrhs) - 2.0*u(i,j,k,nrhs) + u(i+1,j,k,nrhs)
            Huxx(i,j) = Huon(i-1,j,k) - 2.0*Huon(i,j,k) + Huon(i+1,j,k)
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
        cff=1.0/6.0
        DO j=Jstr,Jend
          DO i=IstrU-1,Iend
            UFx(i,j)=0.25*(u(i  ,j,k,nrhs)+                          &


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
            IF (cff1.gt.0.0) THEN
              cff=uxx(i,j)
            ELSE
              cff=uxx(i+1,j)
            END IF
            UFx(i,j)=0.25*(cff1+Gadv*cff)*                           &
     &               (Huon(i  ,j,k)+                                    &
     &                Huon(i+1,j,k)+                                    &
     &                Gadv*0.5*(Huxx(i  ,j)+                         &
     &                             Huxx(i+1,j)))
          END DO
        END DO
#   endif
        DO j=Jstrm1,Jendp1
          DO i=IstrU,Iend
            uee(i,j)=u(i,j-1,k,nrhs)-2.0*u(i,j,k,nrhs)+              &


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
           Hvxx(i,j)=Hvom(i-1,j,k)-2.0*Hvom(i,j,k)+Hvom(i+1,j,k)
          END DO
        END DO




#   ifdef UV_C4ADVECTION
        cff=1.0/6.0
        DO j=Jstr,Jend+1
          DO i=IstrU,Iend
            UFe(i,j)=0.25*(u(i,j  ,k,nrhs)+                          &

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
            IF (cff2.gt.0.0) THEN
              cff=uee(i,j-1)
            ELSE
              cff=uee(i,j)
            END IF
            UFe(i,j)=0.25*(cff1+Gadv*cff)*                           &
     &               (cff2+Gadv*0.5*(Hvxx(i  ,j)+                    &
     &                                  Hvxx(i-1,j)))
          END DO
        END DO
#   endif
        DO j=JstrV,Jend
          DO i=Istrm1,Iendp1
            vxx(i,j)=v(i-1,j,k,nrhs)-2.0*v(i,j,k,nrhs)+              &

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
           Huee(i,j)=Huon(i,j-1,k)-2.0*Huon(i,j,k)+Huon(i,j+1,k)
          END DO
        END DO
#   ifdef UV_C4ADVECTION

    # Fourth-order, centered differences v-momentum horizontal advection.

        cff=1.0/6.0
        DO j=JstrV,Jend
          DO i=Istr,Iend+1
            VFx(i,j)=0.25*(v(i  ,j,k,nrhs)+                          &

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
            IF (cff2.gt.0.0) THEN
              cff=vxx(i-1,j)
            ELSE
              cff=vxx(i,j)
            END IF
            VFx(i,j)=0.25*(cff1+Gadv*cff)*                           &
     &               (cff2+Gadv*0.5*(Huee(i,j  )+                    &
     &                                  Huee(i,j-1)))
          END DO
        END DO
#   endif
        DO j=JstrVm1,Jendp1
          DO i=Istr,Iend
            vee(i,j)=v(i,j-1,k,nrhs)-2.0*v(i,j,k,nrhs)+              &

     &               v(i,j+1,k,nrhs)
            Hvee(i,j)=Hvom(i,j-1,k)-2.0*Hvom(i,j,k)+Hvom(i,j+1,k)
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
        cff=1.0/6.0
        DO j=JstrV-1,Jend
          DO i=Istr,Iend
            VFe(i,j)=0.25*(v(i,j  ,k,nrhs) + v(i,j+1,k,nrhs) - cff*(vee (i,j  ) + vee (i,j+1)))*                       &
     &                    (Hvom(i,j  ,k) + Hvom(i,j+1,k) - cff*(Hvee(i,j  ) + Hvee(i,j+1)))
          END DO
        END DO
#   else
        DO j=JstrV-1,Jend
          DO i=Istr,Iend
            cff1 = v(i,j  ,k,nrhs) + v(i,j+1,k,nrhs)
            IF (cff1.gt.0.0) THEN
              cff=vee(i,j)
            ELSE
              cff=vee(i,j+1)
            END IF

            VFe(i,j)=0.25*(cff1+Gadv*cff)*                           &
     &               (Hvom(i,j  ,k)+                                    &
     &                Hvom(i,j+1,k)+                                    &
     &                Gadv*0.5*(Hvee(i,j  )+                         &
     &                             Hvee(i,j+1)))
          END DO
        END DO
#   endif
#  endif



        # Add in horizontal advection.

        DO j=Jstr,Jend
          DO i=IstrU,Iend
            cff1 = UFx(i,j  ) - UFx(i-1,j)
            cff2 = UFe(i,j+1) - UFe(i,j)
            ru(i,j,k,nrhs) = ru(i,j,k,nrhs) - (cff1+cff2)
          END DO
        END DO

        DO j=JstrV,Jend
          DO i=Istr,Iend
            cff1 = VFx(i+1,j) - VFx(i,j  )
            cff2 = VFe(i  ,j) - VFe(i,j-1)
            rv(i,j,k,nrhs) = rv(i,j,k,nrhs) - (cff1+cff2)

          END DO
        END DO
# endif


      END DO K_LOOP

      J_LOOP : DO j=Jstr,Jend
# ifdef UV_ADV


    # Add in vertical advection.
    # -----------------------------------------------------------------------


        cff1=9.0/16.0
        cff2=1.0/16.0
        DO k=2,N(ng)-2
          DO i=IstrU,Iend
            FC(i,k)=(cff1*(u(i,j,k  ,nrhs) + u(i,j,k+1,nrhs)) - cff2*(u(i,j,k-1,nrhs) + u(i,j,k+2,nrhs)))*                           &
                    (cff1*(W(i  ,j,k) + W(i-1,j,k)) - cff2*(W(i+1,j,k) + W(i-2,j,k)))
          END DO
        END DO

        DO i=IstrU,Iend
          FC(i,N(ng))=0.0
          FC(i,N(ng)-1)=(cff1*(u(i,j,N(ng)-1,nrhs) + u(i,j,N(ng)  ,nrhs)) - cff2*(u(i,j,N(ng)-2,nrhs) + u(i,j,N(ng)  ,nrhs)))*                   &
     &                  (cff1*(W(i  ,j,N(ng)-1) + W(i-1,j,N(ng)-1)) - cff2*(W(i+1,j,N(ng)-1) +W(i-2,j,N(ng)-1)))
          FC(i,1)=(cff1*(u(i,j,1,nrhs) + u(i,j,2,nrhs)) - cff2*(u(i,j,1,nrhs) +  u(i,j,3,nrhs)))*                               &
     &            (cff1*(W(i  ,j,1) + W(i-1,j,1)) - cff2*(W(i+1,j,1) + W(i-2,j,1)))
          FC(i,0)=0.0
        END DO



        DO k=1,N(ng)
          DO i=IstrU,Iend
            cff = FC(i,k) - FC(i,k-1)
            ru(i,j,k,nrhs) = ru(i,j,k,nrhs) - cff
          END DO
        END DO



        IF (j.ge.JstrV) THEN

          cff1=9.0/16.0
          cff2=1.0/16.0
          DO k=2,N(ng)-2
            DO i=Istr,Iend
              FC(i,k)=(cff1*(v(i,j,k  ,nrhs) + v(i,j,k+1,nrhs)) - cff2*(v(i,j,k-1,nrhs) + v(i,j,k+2,nrhs)))*                         &
                      (cff1*(W(i,j  ,k) + W(i,j-1,k)) - cff2*(W(i,j+1,k) + W(i,j-2,k)))
            END DO
          END DO

          DO i=Istr,Iend
            FC(i,N(ng))=0.0
            FC(i,N(ng)-1)=(cff1*(v(i,j,N(ng)-1,nrhs) + v(i,j,N(ng)  ,nrhs)) - cff2*(v(i,j,N(ng)-2,nrhs) + v(i,j,N(ng)  ,nrhs)))*                 &
                          (cff1*(W(i,j  ,N(ng)-1) + W(i,j-1,N(ng)-1)) - cff2*(W(i,j+1,N(ng)-1) W(i,j-2,N(ng)-1)))
            FC(i,1)=(cff1*(v(i,j,1,nrhs) + v(i,j,2,nrhs)) - cff2*(v(i,j,1,nrhs) + v(i,j,3,nrhs)))*                             &
     &              (cff1*(W(i,j  ,1) + W(i,j-1,1)) - cff2*(W(i,j+1,1) + W(i,j-2,1)))
            FC(i,0)=0.0
          END DO

          DO k=1,N(ng)
            DO i=Istr,Iend
              cff=FC(i,k)-FC(i,k-1)
              rv(i,j,k,nrhs)=rv(i,j,k,nrhs)-cff

            END DO
          END DO
        END IF


        # Compute forcing term for the 2D momentum equations.
        # -----------------------------------------------------------------------

        # Vertically integrate baroclinic right-hand-side terms. If not body force stresses, add in the difference between surface and
        # bottom stresses.

        DO i=IstrU,Iend
          rufrc(i,j)=ru(i,j,1,nrhs)
        END DO

        DO k=2,N(ng)
          DO i=IstrU,Iend
            rufrc(i,j) = rufrc(i,j) + ru(i,j,k,nrhs)
          END DO
        END DO


          DO i=Istr,Iend
            cff=om_v(i,j)*on_v(i,j)
            cff1= svstr(i,j)*cff
            cff2=-bvstr(i,j)*cff
            rvfrc(i,j) = rvfrc(i,j)+cff1+cff2

        END IF
      END DO J_LOOP


#endif
