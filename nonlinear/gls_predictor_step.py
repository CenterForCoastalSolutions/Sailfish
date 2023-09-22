
#if defined NONLINEAR && defined GLS_MIXING && defined SOLVE3D
# !=======================================================================
# !  Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
# !    Licensed under a MIT/X style license           Hernan G. Arango   !
# !    See License_ROMS.txt                   Alexander F. Shchepetkin   !
# !==================================================== John C. Warner ===
# !                                                                      !
# !  This routine perfoms the predictor step for turbulent kinetic       !
# !  energy prognostic variables, tke and gls. A NON-conservative,       !
# !  but constancy preserving, auxiliary advective substep for tke       !
# !  gls equations is carried out. The result of this substep will       !
# !  be used to compute advective terms in the corrector substep.        !
# !  No dissipation terms are included here.                             !
# !                                                                      !
# !=======================================================================

def gls_prestep:

# !-----------------------------------------------------------------------
# !  Predictor step for advection of turbulent kinetic energy variables.
# !-----------------------------------------------------------------------
# !
# ! Start computation of auxiliary time step fields tke(:,:,:,n+1/2) and
# ! gls(:,:,:,n+1/2) with computation of horizontal advection terms and
# ! auxiliary grid-box height field Hz_new()=Hz(:,:,k+1/2,n+1/2);
# ! This is effectivey an LF step with subsequent interpolation of the
# ! result half step back, using AM3 weights. The LF step and
# ! interpolation are perfomed as a single operation, which results in
# ! weights cff1,cff2,cff3 below.
# !
# ! Either centered fourth-order accurate or standard second order
# ! accurate versions are supported.
# !
# ! At the same time prepare for corrector step for tke,gls: set tke,
# ! gls(:,:,:,nnew) to  tke,gls(:,:,:,nstp) multiplied by the
# ! corresponding grid-box height. This needs done at this time because
# ! array Hz(:,:,:) will overwritten after 2D time stepping with the
# ! values computed from zeta(:,:,n+1) rather than zeta(:,:,n), so that
# ! the old-time-step Hz will be no longer available.
# !


    if K_C2ADVECTION:

# !
# !  Second-order, centered differences advection.
# !

            XF  = WtoR(Huon)
            FX  = XF*RtoU(tke[:,:,:,nspt])
            FXL = XF*RtoU(gls[:,:,:,nspt])

            EF  = WtoR(Hvom)
            FE  = EF*RtoV(tke[:,:,:,nspt])
            FEL = EF*RtoV(gls[:,:,:,nspt])

    else
# !
# !  Fourth-order, centered differences advection.
# !
#         DO j=Jstr,Jend
#           DO i=Istrm1,Iendp2
            grad (i,j) = tke(i,j,k,nstp) - tke(i-1,j,k,nstp)
            gradL(i,j) = gls(i,j,k,nstp) - gls(i-1,j,k,nstp)

#           END DO
#         END DO


        # BOUNDARY THINGS
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            DO j=Jstr,Jend
              grad (Istr-1,j)=grad (Istr,j)
              gradL(Istr-1,j)=gradL(Istr,j)
            END DO
          END IF

          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO j=Jstr,Jend
              grad (Iend+2,j)=grad (Iend+1,j)
              gradL(Iend+2,j)=gradL(Iend+1,j)
            END DO
          END IF



        cff=1.0/6.0
        DO j=Jstr,Jend
          DO i=Istr,Iend+1
            XF(i,j) = 0.5*(Huon(i,j,k) + Huon(i,j,k+1))
            FX (i,j)=XF(i,j)*0.5*(tke(i-1,j,k,nstp) + tke(i,j,k,nstp) - cff*( grad(i+1,j) -  grad(i-1,j)))
            FXL(i,j)=XF(i,j)*0.5*(gls(i-1,j,k,nstp)  +gls(i,j,k,nstp) - cff*(gradL(i+1,j) - gradL(i-1,j)))
          END DO
        END DO

        DO j=Jstrm1,Jendp2
          DO i=Istr,Iend
            grad (i,j) = (tke(i,j,k,nstp) - tke(i,j-1,k,nstp))
            gradL(i,j) = (gls(i,j,k,nstp) - gls(i,j-1,k,nstp))

          END DO
        END DO


          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO i=Istr,Iend
              grad (i,Jstr-1)=grad (i,Jstr)
              gradL(i,Jstr-1)=gradL(i,Jstr)
            END DO
          END IF
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO i=Istr,Iend
              grad (i,Jend+2)=grad (i,Jend+1)
              gradL(i,Jend+2)=gradL(i,Jend+1)
            END DO
          END IF


        cff=1.0/6.0
        DO j=Jstr,Jend+1
          DO i=Istr,Iend
            EF(i,j)=0.5*(Hvom(i,j,k)+Hvom(i,j,k+1))
            FE (i,j)=EF(i,j)*0.5*(tke(i,j-1,k,nstp) + tke(i,j,k,nstp) - cff*(grad (i,j+1) - grad (i,j-1)))
            FEL(i,j)=EF(i,j)*0.5*(gls(i,j-1,k,nstp) + gls(i,j,k,nstp) - cff*(gradL(i,j+1) - gradL(i,j-1)))
          END DO
        END DO
# endif

# !
# !  Time-step horizontal advection.
# !
        IF (iic(ng).eq.ntfirst(ng)) THEN
          cff1=1.0
          cff2=0.0
          cff3=0.5*dt(ng)
          indx=nstp
        ELSE
          cff1=0.5 + Gamma
          cff2=0.5 - Gamma
          cff3=(1.0 - Gamma)*dt
          indx=3 - nstp
        END IF

        DO j=Jstr,Jend
          DO i=Istr,Iend
            cff = 0.5*(Hz(i,j,k) + Hz(i,j,k+1))
            cff4 = cff3*pm(i,j)*pn(i,j)
            Hz_half(i,j,k) = cff - cff4*(XF(i+1,j) - XF(i,j) + EF(i,j+1) - EF(i,j))
            tke(i,j,k,3) = cff*(cff1*tke(i,j,k,nstp) + cff2*tke(i,j,k,indx)) -
                           cff4*(FX (i+1,j) - FX (i,j) + FE (i,j+1) - FE (i,j))
            gls(i,j,k,3) = cff*(cff1*gls(i,j,k,nstp) + cff2*gls(i,j,k,indx)) -
                           cff4*(FXL(i+1,j) - FXL(i,j) + FEL(i,j+1) - FEL(i,j))
            tke(i,j,k,nnew) = cff*tke(i,j,k,nstp)
            gls(i,j,k,nnew) = cff*gls(i,j,k,nstp)
          END DO
        END DO
      END DO

    # Compute vertical advection term.

      DO j=Jstr,Jend

        cff1=7.0/12.0
        cff2=1.0/12.0
        DO k=2,N(ng)-1
          DO i=Istr,Iend
            CF(i,k)=0.5*(W(i,j,k)+W(i,j,k-1))

            FC (i,k) = CF(i,k)*(cff1*(tke(i,j,k-1,nstp) + tke(i,j,k  ,nstp)) -
                                cff2*(tke(i,j,k-2,nstp) + tke(i,j,k+1,nstp)))
            FCL(i,k) = CF(i,k)*(cff1*(gls(i,j,k-1,nstp) + gls(i,j,k  ,nstp)) -
                                cff2*(gls(i,j,k-2,nstp) + gls(i,j,k+1,nstp)))
          END DO
        END DO

        # Edge cases (top and down)
        cff1=1.0/3.0
        cff2=5.0/6.0
        cff3=1.0/6.0
        DO i=Istr,Iend
          CF(i,1)=0.5*(W(i,j,0)+W(i,j,1))

          FC (i,1) = CF(i,1)*(cff1*tke(i,j,0,nstp) + cff2*tke(i,j,1,nstp) - cff3*tke(i,j,2,nstp))
          FCL(i,1) = CF(i,1)*(cff1*gls(i,j,0,nstp) + cff2*gls(i,j,1,nstp) - cff3*gls(i,j,2,nstp))
          CF(i,N(ng)) = 0.5*(W(i,j,N(ng))+W(i,j,N(ng)-1))

          FC (i,N(ng)) = CF(i,N(ng))*(cff1*tke(i,j,N(ng)  ,nstp) + cff2*tke(i,j,N(ng)-1,nstp) - cff3*tke(i,j,N(ng)-2,nstp))
          FCL(i,N(ng)) = CF(i,N(ng))*(cff1*gls(i,j,N(ng)  ,nstp) + cff2*gls(i,j,N(ng)-1,nstp) - cff3*gls(i,j,N(ng)-2,nstp))
        END DO
# endif

    # Time-step vertical advection term.
        IF (iic(ng).eq.ntfirst(ng)) THEN
          cff3 = 0.5*dt
        ELSE
          cff3 = (1.0 - Gamma)*dt
        END IF

        DO k=1,N(ng)-1
          DO i=Istr,Iend
            cff4=cff3*pm(i,j)*pn(i,j)
            Hz_half(i,j,k)=Hz_half(i,j,k)-cff4*(CF(i,k+1)-CF(i,k))
            cff1=1.0_r8/Hz_half(i,j,k)
            tke(i,j,k,3) = cff1*(tke(i,j,k,3) - cff4*(FC (i,k+1) - FC (i,k)))
            gls(i,j,k,3) = cff1*(gls(i,j,k,3) - cff4*(FCL(i,k+1) - FCL(i,k)))
          END DO
        END DO
      END DO

    # Apply lateral boundary conditions.
      CALL tkebc_tile (ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj, N(ng),                       &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
     &                 3, nstp,                                         &
     &                 gls, tke)



