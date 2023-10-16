

    // The Eq. are:
    // $$\frac{\partial\boldsymbol{u}}{\partial t}+\frac{\partial}{\partial z}K_{v}\frac{\partial\boldsymbol{u}}{\partial z}=\boldsymbol{ru}$$
    // Which can be discretized as (for each component)
    // $$\frac{u_{k}^{n+1}-u_{k}^{n}}{\Delta t}+\frac{1}{\Delta z_{k}}\left(K_{v}^{^{k-\nicefrac{1}{2}}}\frac{u_{k}^{n+1}-u_{k-1}^{n+1}}{\Delta z_{k-\nicefrac{1}{2}}}-K_{v}^{^{k+\nicefrac{1}{2}}}\frac{u_{k+1}^{n+1}-u_{k}^{n+1}}{\Delta z_{k+\nicefrac{1}{2}}}\right)=ru$$
    // After some algebra:
    // $$\underbrace{\left(\Delta z_{k}+\underbrace{\frac{\Delta tK_{v}^{^{k-\nicefrac{1}{2}}}}{\Delta z_{k-\nicefrac{1}{2}}}}_{-FC[k-1]}+\underbrace{\frac{\Delta tK_{v}^{^{k+\nicefrac{1}{2}}}}{\Delta z_{k+\nicefrac{1}{2}}}}_{-FC[k]}\right)}_{BC[k]}u_{k}^{n+1}\underbrace{-\Delta t\frac{K_{v}^{^{k-\nicefrac{1}{2}}}}{\Delta z_{k-\nicefrac{1}{2}}}}_{FC[k-1]}u_{k-1}^{n+1}\underbrace{-\Delta t\frac{K_{v}^{^{k+\nicefrac{1}{2}}}}{\Delta z_{k+\nicefrac{1}{2}}}}_{FC[k]}u_{k+1}^{n+1}	=\underbrace{\Delta z_{k}\left(\Delta t\,ru+u^{n}\right)}_{DC[k]}$$
    // The matrix equation is $M\,U=DC$, where:
    // $$M=\left(\begin{array}{ccccccc}BC[1] & FC[1]\\FC[1] & BC[2] & FC[2]\\ & FC[2] & BC[3] & FC[3]\\ &  & \ddots & \ddots & \ddots\\ &  &  & \ddots & \ddots & \ddots\\ &  &  &  & FC[N-2] & BC[N-1] & FC[N-1]\\ &  &  &  &  & FC[N-1] & BC[N]\end{array}\right)$$
    // In reality, the terms with η or ξ derivatives are divided bi m and n respectively, while the terms with no derivatives are divided by m*n.


//    u  = u(nnew,:,:,:)
//    ru = ru(nrhs,:,:,:)

void XXXXXXX(const double Δt, const double lambda)
{
    const double coef = ;

    double *FC = bufFC[i*N];

    // If the index is not a U-node leaves.
    if (!isUnode(i)) return;

//      STENCIL3D(h);

    // Builds matrix M and vector DC.
    //-------------------------------
    const auto AKvU = RtoU(Akv);
    const auto HzU  = RtoU(Hz);

    const auto pmnU = RtoU(pm)*RtoU(pn)

    const auto cΔt_mn = coef*Δt*RtoU(pm)*RtoU(pn);

    FC[0] = 0.0;
    FC[N] = 0.0;
    for (int k=1; k<N; k++)
    {
        // Δz is the vertical distance between two U nodes.
        const double Δz = RtoU(z_r(k+1,0,0) - z_r(k,0,0));

        // RHS
        DC[k] = u(nnew,k,0,0) + cΔt_mn*ru(k,0,0);

        // Off-diagonal elements
        FC[k] = -lambda*Δt*AKvU(k,0,0)/Δz;

        //Diagonal elements.
        BC[k] = HzU(0,0,k) - FC[k] - FC[k-1];
    }
    BC[N] = HzU(0,0,N) - FC[N] - FC[N-1];


    // Up to here, DC contains the RHS, BC is the main diagonal and FC[1:-1], FC[:,-2] the other two diagonals.


    // Solve tridiagonal system.
    // -------------------------

    // See https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm.

    // Forward substitution.
    CF[1] = FC[1]/BC[1];
    DC[1] = DC[1]/BC[1];

    for (int k=2; k<=N; k++)
    {
        double denom = 1.0/(BC[k] - FC[k-1]*CF[k-1]);

        CF(i,k) = denom*FC[k];
        DC(i,k) = denom*(DC[k] - FC[k-1]*DC[k-1]);
    }


    // Back substitution.
    DC[N] -= FC[N-1]*DC[N-1]/(BC[N] - FC[N-1]*CF[N-1]);
    u(0,0,N,nnew) = DC[N];

    for (int k=N-1; k>1; k--)
    {
        u(0,0,k,nnew) = DC[k] - CF[k]*u(0,0,k+1,nnew);
    }





# This subroutine time-steps the nonlinear  horizontal  momentum
# equations. The vertical viscosity terms are time-stepped using
# an implicit algorithm.

def step3d_UV():

    # Time step momentum equations
    # ----------------------------

    # ξ-direction.
    AK  = RtoU(Akv)
    Hzk = RtoU(Hz)
    DO j = Jstr, Jend
        # DO i=IstrU,Iend
        #
        #   AK(i,0) = 0.5*(Akv(i-1,j,0)+ Akv(i  ,j,0))
        #   DO k=1,N(ng)
        #     AK (i,k) = 0.5*(Akv(i-1,j,k) + Akv(i,j,k))
        #     Hzk(i,k) = 0.5*(Hz (i-1,j,k) + Hz (i,j,k))
        #
        #
        #   END DO
        # END DO


    # Time step right-hand-side terms.



        # Create vertical viscous term operator matrix.
        #         IF (iic(ng).eq.ntfirst(ng)) THEN
        #   cff = 0.25*dt
        # ELSE IF (iic(ng).eq.(ntfirst(ng)+1)) THEN
        #   cff = 0.25*dt*3.0/2.0
        # ELSE
        #   cff = 0.25*dt*23.0/12.0
        # END IF


    #     DO i=IstrU,Iend
    #       DC(i,0)=cff*(pm(i,j) + pm(i-1,j))*(pn(i,j) + pn(i-1,j))
    #     END DO
    #
    #     DO k=1,N
    #       DO i=IstrU,Iend
    #         u(i,j,k, nnew) = u(i,j,k, nnew) + DC(i,0)*ru(i,j,k, nrhs)
    #       END DO
    #     END DO
    #
    #
    #
    # # Compute off-diagonal coefficients [lambda*dt*Akv/Hz] for the implicit vertical viscosity term at horizontal
    # # U-points and vertical W-points.
    #
    #     cff = -lambda*dt/0.5
    #     DO k=1,N(ng)-1
    #       DO i=IstrU,Iend
    #         cff1 = 1.0/(z_r(i,j,k+1) + z_r(i-1,j,k+1) - z_r(i,j,k  ) - z_r(i-1,j,k  ))
    #         FC(i,k)=cff*cff1*AK(i,k)
    #       END DO
    #     END DO
    #
    #     DO i=IstrU,Iend
    #       FC(i,0)=0.0
    #       FC(i,N(ng))=0.0
    #     END DO

    MvU = createVertViscousOpMatrix()


    # Solve the tridiagonal system.
    u = solveTri(MvV, AK, z_r, u)
    # Solve the tridiagonal system.

        # DO k=1,N(ng)
        #   DO i=IstrU,Iend
        #     DC(i,k)=u(i,j,k,nnew)
        #     BC(i,k)= Hzk(i,k) - FC(i,k) - FC(i,k-1)
        #   END DO
        # END DO
        #
        # DO i=IstrU,Iend
        #   cff=1.0/BC(i,1)
        #   CF(i,1)=cff*FC(i,1)
        #   DC(i,1)=cff*DC(i,1)
        # END DO
        #
        # DO k=2,N(ng)-1
        #   DO i=IstrU,Iend
        #     cff=1.0_r8/(BC(i,k)-FC(i,k-1)*CF(i,k-1))
        #     CF(i,k)=cff*FC(i,k)
        #     DC(i,k)=cff*(DC(i,k)-FC(i,k-1)*DC(i,k-1))
        #   END DO
        # END DO
        #
        # # Compute new solution by back substitution.
        #
        # DO i=IstrU,Iend
        #
        #   DC(i,N(ng))=(DC(i,N(ng))-FC(i,N(ng)-1)*DC(i,N(ng)-1))/(BC(i,N(ng))-FC(i,N(ng)-1)*CF(i,N(ng)-1))
        #   u(i,j,N(ng),nnew) = DC(i,N(ng))
        #
        # END DO
        #
        #
        # DO k=N(ng)-1,1,-1
        #   DO i=IstrU,Iend
        #
        #     DC(i,k) = DC(i,k) - CF(i,k)*DC(i,k+1)
        #     u(i,j,k,nnew) = DC(i,k)
        #
        #   END DO
        # END DO


    # CORRECTS the baroclinic velocity average
    pass
    # Replace INTERIOR POINTS incorrect vertical mean with more accurate barotropic component, ubar=DU_avg1/(D*on_u).
    # Recall that, D=CF(:,0).
    # INTERIOR POINTS here means not the point at the bottom. In other words, ensures that the average baroclinic velocity = barotropic velocity.



    #     DO i=IstrU,Iend
    #       CF(i,0)=Hzk(i,1)
    #       DC(i,0)=u(i,j,1,nnew)*Hzk(i,1)
    #
    #     END DO
    #     DO k=2,N(ng)
    #       DO i=IstrU,Iend
    #         CF(i,0)=CF(i,0)+Hzk(i,k)
    #         DC(i,0)=DC(i,0)+u(i,j,k,nnew)*Hzk(i,k)
    #
    #       END DO
    #     END DO
    #     DO i=IstrU,Iend
    #       cff1=1.0/(CF(i,0)*on_u(i,j))
    #       DC(i,0)=(DC(i,0)*on_u(i,j) - DU_avg1(i,j))*cff1     # recursive
    #
    #     END DO
    #
    # # Couple and update new solution.
    #
    #     DO k=1,N(ng)
    #       DO i=IstrU,Iend
    #         u(i,j,k, nnew) = u(i,j,k,nnew) - DC(i,0)
    #       END DO
    #     END DO


    pass


    MvV = createVertViscousOpMatrix()
    #     IF (j.ge.JstrV) THEN
    #       DO i=Istr,Iend
    #         AK(i,0) = 0.5*(Akv(i,j-1,0) + Akv(i,j  ,0))
    #         DO k=1,N(ng)
    #           AK(i,k) = 0.5*(Akv(i,j-1,k) + Akv(i,j  ,k))
    #           Hzk(i,k) = 0.5*(Hz(i,j-1,k) + Hz(i,j  ,k))
    #       END DO
    #
    # # Time step right-hand-side terms.
    #
    #       IF (iic(ng).eq.ntfirst(ng)) THEN
    #         cff=0.25*dt(ng)
    #       ELSE IF (iic(ng).eq.(ntfirst(ng)+1)) THEN
    #         cff=0.25*dt(ng)*3.0/2.0
    #       ELSE
    #         cff=0.25*dt(ng)*23.0/12.0
    #       END IF


    # Create vertical viscous term operator matrix.

    #       DO i=Istr,Iend
    #         DC(i,0) = cff*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
    #       END DO
    #
    #       DO k=1,N(ng)
    #         DO i=Istr,Iend
    #           v(i,j,k,nnew)=v(i,j,k,nnew)+DC(i,0)*rv(i,j,k,nrhs)
    #         END DO
    #       END DO
    #
    #
    #
    # # Compute off-diagonal coefficients [lambda*dt*Akv/Hz] for the
    # # implicit vertical viscosity term at horizontal V-points and
    # # vertical W-points.
    #
    #       cff=-lambda*dt/0.5
    #       DO k=1,N(ng)-1
    #         DO i=Istr,Iend
    #           cff1=1.0_r8/(z_r(i,j,k+1) + z_r(i,j-1,k+1) - z_r(i,j,k  ) - z_r(i,j-1,k  ))
    #           FC(i,k)=cff*cff1*AK(i,k)
    #         END DO
    #       END DO
    #
    #       DO i=Istr,Iend
    #         FC(i,0)=0.0_r8
    #         FC(i,N(ng))=0.0_r8
    #       END DO

    # Solve the tridiagonal system.
    v = solveTri(MvV, AK, z_r, v)

    #       DO k=1,N(ng)
    #         DO i=Istr,Iend
    #           DC(i,k)=v(i,j,k,nnew)
    #           BC(i,k)=Hzk(i,k)-FC(i,k)-FC(i,k-1)
    #         END DO
    #       END DO
    #
    #       DO i=Istr,Iend
    #         cff=1.0/BC(i,1)
    #         CF(i,1)=cff*FC(i,1)
    #         DC(i,1)=cff*DC(i,1)
    #       END DO
    #
    #       DO k=2,N(ng)-1
    #         DO i=Istr,Iend
    #           cff=1.0/(BC(i,k)-FC(i,k-1)*CF(i,k-1))
    #           CF(i,k)=cff*FC(i,k)
    #           DC(i,k)=cff*(DC(i,k)-FC(i,k-1)*DC(i,k-1))
    #         END DO
    #       END DO
    #
    # # Compute new solution by back substitution.
    #       DO i=Istr,Iend
    #         DC(i,N(ng))=(DC(i,N(ng)) - FC(i,N(ng)-1)*DC(i,N(ng)-1))/ (BC(i,N(ng)) - FC(i,N(ng)-1)*CF(i,N(ng)-1))
    #         v(i,j,N(ng),nnew) = DC(i,N(ng))
    #
    #       END DO
    #       DO k=N(ng)-1,1,-1
    #         DO i=Istr,Iend
    #
    #           DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
    #           v(i,j,k,nnew)=DC(i,k)
    #
    #         END DO
    #       END DO


    # CORRECTS the baroclinic velocity average

    correctBaroclinicMeanVel()
    # Replace INTERIOR POINTS incorrect vertical mean with more accurate barotropic component, vbar=DV_avg1/(D*om_v).
    # Recall that, D=CF(:,0).

    #       DO i=Istr,Iend
    #         CF(i,0)=Hzk(i,1)
    #         DC(i,0)=v(i,j,1,nnew)*Hzk(i,1)
    #
    #       END DO
    #       DO k=2,N(ng)
    #         DO i=Istr,Iend
    #           CF(i,0)=CF(i,0) + Hzk(i,k)
    #           DC(i,0)=DC(i,0) + v(i,j,k,nnew)*Hzk(i,k)
    #
    #         END DO
    #       END DO

    #       DO i=Istr,Iend
    #         DC(i,0) = (DC(i,0)*om_v(i,j) - DV_avg1(i,j))/(CF(i,0)*om_v(i,j))    ! recursive
    #
    #       END DO
    #
    # # Couple and update new solution.
    #
    #       DO k=1,N(ng)
    #         DO i=Istr,Iend
    #           v(i,j,k,nnew) = v(i,j,k,nnew) - DC(i,0)
    #
    #
    #         END DO
    #       END DO
    #
    #
    #     END IF
    #   END DO


    setBCs______()
    # # Set lateral boundary conditions.
    # # -----------------------------------------------------------------------
    #
    #   CALL u3dbc_tile (LBi, UBi, LBj, UBj, N(ng),                       &
    #  &                 IminS, ImaxS, JminS, JmaxS,                      &
    #  &                 nstp, nnew,                                      &
    #  &                 u)
    #   CALL v3dbc_tile (LBi, UBi, LBj, UBj, N(ng),                       &
    #  &                 IminS, ImaxS, JminS, JmaxS,                      &
    #  &                 nstp, nnew,                                      &
    #  &                 v)
    #
    #
    #

    applyPointSources()
    # # Apply momentum transport point sources (like river runoff), if any.
    # # -----------------------------------------------------------------------
    #   IF (LuvSrc(ng)) THEN
    #     DO is=1,Nsrc(ng)
    #       i=SOURCES(ng)%Isrc(is)
    #       j=SOURCES(ng)%Jsrc(is)
    #       IF (((IstrR.le.i).and.(i.le.IendR)).and.                      &
    #  &        ((JstrR.le.j).and.(j.le.JendR))) THEN
    #         IF (INT(SOURCES(ng)%Dsrc(is)).eq.0) THEN
    #           DO k=1,N(ng)
    #             cff1=1.0_r8/(on_u(i,j)*                                 &
    #  &                       0.5_r8*(z_w(i-1,j,k)-z_w(i-1,j,k-1)+       &
    #  &                               z_w(i  ,j,k)-z_w(i  ,j,k-1)))
    #             u(i,j,k,nnew)=SOURCES(ng)%Qsrc(is,k)*cff1
    #           END DO
    #         ELSE IF (INT(SOURCES(ng)%Dsrc(is)).eq.1) THEN
    #           DO k=1,N(ng)
    #             cff1=1.0_r8/(om_v(i,j)*                                 &
    #  &                       0.5_r8*(z_w(i,j-1,k)-z_w(i,j-1,k-1)+       &
    #  &                               z_w(i,j  ,k)-z_w(i,j  ,k-1)))
    #             v(i,j,k,nnew)=SOURCES(ng)%Qsrc(is,k)*cff1
    #           END DO
    #         END IF
    #       END IF
    #     END DO
    #   END IF


    pass
    # Couple 2D and 3D momentum equations.
    # -----------------------------------------------------------------------

    correctBarotropicVel()
#     # Couple velocity component in the XI-direction.
#       DO j=JstrT,JendT
#         DO i=IstrP,IendT
#           DC(i,0)=0.0
#           CF(i,0)=0.0
#           FC(i,0)=0.0
#         END DO
#
#     # Compute thicknesses of U-boxes DC(i,1:N), total depth of the water column DC(i,0), and incorrect vertical mean CF(i,0).
#     # Notice that barotropic component is replaced with its fast-time averaged values.
#     #  DC = integrated height, CF integrated flux.
#         DO k=1,N(ng)
#           DO i=IstrP,IendT
#             DC(i,k)=on_u(i,j)*RtoU(Hz(i,j,k))
#             DC(i,0)=DC(i,0) + DC(i,k)
#             CF(i,0)=CF(i,0) + DC(i,k)*u(i,j,k,nnew)
#           END DO
#         END DO
#
#         DO i=IstrP,IendT
#           DC(i,0)=1.0/DC(i,0)                           # recursive
#           CF(i,0)=DC(i,0)*(CF(i,0)-DU_avg1(i,j))        # recursive
#
#           ubar(i,j,1) = DC(i,0)*DU_avg1(i,j)
#           ubar(i,j,2) = ubar(i,j,1)
#
#         END DO
#
#
#   # Replace only BOUNDARY POINTS incorrect vertical mean with more accurate barotropic component, ubar = DU_avg1/(D*on_u). Recall that, D=CF(:,0).
#
# # NOTE:  Only the BOUNDARY POINTS need to be replaced. Avoid redundant update in the interior again for computational purposes which
# #        will not affect the nonlinear code.  However, the adjoint code is wrong because the interior solution is corrected
# #        twice. The replacement is avoided if the boundary edge is periodic. The J-loop is pipelined, so we need to do a special
# #        test for the southern and northern domain edges.
# #        IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
#           IF (DOMAIN(ng)%Western_Edge(tile)) THEN
#             DO k=1,N(ng)
#               u(Istr,j,k,nnew)=u(Istr,j,k,nnew)-CF(Istr,0)
#             END DO
#           END IF
#         END IF
#
#           IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
#             DO k=1,N(ng)
#               u(Iend+1,j,k,nnew)=u(Iend+1,j,k,nnew)-CF(Iend+1,0)
#             END DO
#           END IF
# !
#           IF (j.eq.0) THEN                           ! southern boundary
#             DO k=1,N(ng)                             ! J-loop pipelined
#               DO i=IstrU,Iend
#                 u(i,j,k,nnew)=u(i,j,k,nnew)-CF(i,0)
#               END DO
#             END DO
#           END IF
# !
#           IF (j.eq.Mm(ng)+1) THEN                    # northern boundary
#             DO k=1,N(ng)                             # J-loop pipelined
#               DO i=IstrU,Iend
#                 u(i,j,k,nnew) = u(i,j,k,nnew) - CF(i,0)
#               END DO
#             END DO
#           END IF
#
#
#         # Compute correct mass flux, Hz*u/n.
#
#         DO k=N(ng),1,-1
#           DO i=IstrP,IendT
#             Huon(i,j,k) = 0.5*(Huon(i,j,k) + u(i,j,k,nnew)*DC(i,k))
#           END DO
#         END DO
#
#         DO i=IstrP,IendT
#           FC(i,0) = DC(i,0)*(FC(i,0) - DU_avg2(i,j))        # recursive
#         END DO
#
#         DO k=1,N(ng)
#           DO i=IstrP,IendT
#             Huon(i,j,k) = Huon(i,j,k) - DC(i,k)*FC(i,0)
#           END DO
#         END DO
#
#
#
#
#
#     # Couple velocity component in the ETA-direction.
#
#         IF (j.ge.Jstr) THEN
#           DO i=IstrT,IendT
#             DC(i,0)=0.0
#             CF(i,0)=0.0
#
#             FC(i,0)=0.0
#           END DO
#
#
#     # Compute thicknesses of V-boxes DC(i,1:N), total depth of the water column DC(i,0), and incorrect vertical mean CF(i,0).
#     # Notice that barotropic component is replaced with its fast-time averaged values.
#     #  DC = integrated height, CF integrated flux.
#           DO k=1,N(ng)
#             DO i=IstrT,IendT
#               DC(i,k)=on_v(i,j)*RtoV(Hz(i,j,k))
#               DC(i,0) = DC(i,0) + DC(i,k)
#               CF(i,0) = CF(i,0) + DC(i,k)*v(i,j,k,nnew)
#             END DO
#           END DO
#
#           DO i=IstrT,IendT
#             DC(i,0) = 1.0/DC(i,0)                          # recursive
#             CF(i,0) = DC(i,0)*(CF(i,0)-DV_avg1(i,j))       # recursive
#
#             vbar(i,j,1) = DC(i,0)*DV_avg1(i,j)
#
#             vbar(i,j,2) = vbar(i,j,1)
#
#           END DO
#
#     # Replace only BOUNDARY POINTS incorrect vertical mean with more accurate barotropic component, vbar=DV_avg1/(D*om_v).  Recall that,
#     # D=CF(:,0).
#
#     # NOTE:  Only the BOUNDARY POINTS need to be replaced. Avoid redundant update in the interior again for computational purposes which
#     #        will not affect the nonlinear code.  However, the adjoint code is wrong because the interior solution is corrected
#     #        twice. The replacement is avoided if the boundary edge is periodic. The J-loop is pipelined, so we need to do a special
#     #        test for the southern and northern domain edges.
#             IF (DOMAIN(ng)%Western_Edge(tile)) THEN
#               DO k=1,N(ng)
#                 v(Istr-1,j,k,nnew)=v(Istr-1,j,k,nnew)-CF(Istr-1,0)
#               END DO
#            END IF
#
#             IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
#               DO k=1,N(ng)
#                 v(Iend+1,j,k,nnew)=v(Iend+1,j,k,nnew)-CF(Iend+1,0)
#               END DO
#             END IF
#           END IF
#
#             IF (j.eq.1) THEN                         ! southern boundary
#               DO k=1,N(ng)                           ! J-loop pipelined
#                 DO i=Istr,Iend
#                   v(i,j,k,nnew)=v(i,j,k,nnew)-CF(i,0)
#                 END DO
#               END DO
#           END IF
#
#             IF (j.eq.Mm(ng)+1) THEN                  ! northern boundary
#               DO k=1,N(ng)                           ! J-loop pipelined
#                 DO i=Istr,Iend
#                   v(i,j,k,nnew)=v(i,j,k,nnew)-CF(i,0)
#                 END DO
#               END DO
#           END IF
#
#     # Compute correct mass flux, Hz*v/m.
#
#           DO k=N(ng),1,-1
#             DO i=IstrT,IendT
#               Hvom(i,j,k)=0.5*(Hvom(i,j,k)+v(i,j,k,nnew)*DC(i,k))
#
#               FC(i,0)=FC(i,0)+Hvom(i,j,k)
#
#             END DO
#           END DO
#           DO i=IstrT,IendT
#             FC(i,0)=DC(i,0)*(FC(i,0)-DV_avg2(i,j))      ! recursive
#           END DO
#           DO k=1,N(ng)
#             DO i=IstrT,IendT
#               Hvom(i,j,k)=Hvom(i,j,k)-DC(i,k)*FC(i,0)
#             END DO
#           END DO
#         END IF
#       END DO
