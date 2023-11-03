

// The Eq. are:
// $$\frac{\partial\boldsymbol{u}}{\partial t}+\frac{\partial}{\partial z}K_{v}\frac{\partial\boldsymbol{u}}{\partial z}=\boldsymbol{ru}$$
// Which can be discretized as (for each component)
// $$\frac{u_{k}^{n+1}-u_{k}^{n}}{\Delta t}+\frac{1}{\Delta z_{k}}\left(K_{v}^{^{k-\nicefrac{1}{2}}}\frac{u_{k}^{n+1}-u_{k-1}^{n+1}}{\Delta z_{k-\nicefrac{1}{2}}}-K_{v}^{^{k+\nicefrac{1}{2}}}\frac{u_{k+1}^{n+1}-u_{k}^{n+1}}{\Delta z_{k+\nicefrac{1}{2}}}\right)=ru$$
// After some algebra:
// $$\underbrace{\left(\Delta z_{k}+\underbrace{\frac{\Delta tK_{v}^{^{k-\nicefrac{1}{2}}}}{\Delta z_{k-\nicefrac{1}{2}}}}_{-FC[k-1]}+\underbrace{\frac{\Delta tK_{v}^{^{k+\nicefrac{1}{2}}}}{\Delta z_{k+\nicefrac{1}{2}}}}_{-FC[k]}\right)}_{BC[k]}u_{k}^{n+1}\underbrace{-\Delta t\frac{K_{v}^{^{k-\nicefrac{1}{2}}}}{\Delta z_{k-\nicefrac{1}{2}}}}_{FC[k-1]}u_{k-1}^{n+1}\underbrace{-\Delta t\frac{K_{v}^{^{k+\nicefrac{1}{2}}}}{\Delta z_{k+\nicefrac{1}{2}}}}_{FC[k]}u_{k+1}^{n+1}	=\underbrace{\Delta z_{k}\left(\Delta t\,ru+u^{n}\right)}_{RHS[k]}$$
// The matrix equation is $M\,U=RHS$, where:
// $$M=\left(\begin{array}{ccccccc}BC[1] & FC[1]\\FC[1] & BC[2] & FC[2]\\ & FC[2] & BC[3] & FC[3]\\ &  & \ddots & \ddots & \ddots\\ &  &  & \ddots & \ddots & \ddots\\ &  &  &  & FC[N-2] & BC[N-1] & FC[N-1]\\ &  &  &  &  & FC[N-1] & BC[N]\end{array}\right)$$
// In reality, the terms with η or ξ derivatives are divided bi m and n respectively, while the terms with no derivatives are divided by m*n.

void applyPointSources()
{
//# # Apply momentum transport point sources (like river runoff), if any.
//    # # -----------------------------------------------------------------------
//    #   IF (LuvSrc(ng)) THEN
//    #     DO is=1,Nsrc(ng)
//    #       i=SOURCES(ng)%Isrc(is)
//    #       j=SOURCES(ng)%Jsrc(is)
//    #       IF (((IstrR.le.i).and.(i.le.IendR)).and.                      &
//    #  &        ((JstrR.le.j).and.(j.le.JendR))) THEN
//    #         IF (INT(SOURCES(ng)%Dsrc(is)).eq.0) THEN
//    #           DO k=1,N(ng)
//    #             cff1=1.0_r8/(on_u(i,j)*                                 &
//    #  &                       0.5_r8*(z_w(i-1,j,k)-z_w(i-1,j,k-1)+       &
//    #  &                               z_w(i  ,j,k)-z_w(i  ,j,k-1)))
//    #             u(i,j,k,nnew)=SOURCES(ng)%Qsrc(is,k)*cff1
//    #           END DO
//    #         ELSE IF (INT(SOURCES(ng)%Dsrc(is)).eq.1) THEN
//    #           DO k=1,N(ng)
//    #             cff1=1.0_r8/(om_v(i,j)*                                 &
//    #  &                       0.5_r8*(z_w(i,j-1,k)-z_w(i,j-1,k-1)+       &
//    #  &                               z_w(i,j  ,k)-z_w(i,j  ,k-1)))
//    #             v(i,j,k,nnew)=SOURCES(ng)%Qsrc(is,k)*cff1
//    #           END DO
//    #         END IF
//    #       END IF
//    #     END DO
//    #   END IF
}

template<typename RtoU>
void correctBarotropicVel()
{
    RtoU
    on_u
    u(i,j,k,nnew)


    // Compute thicknesses of U-boxes DC(i,1:N), total depth of the water column DC(i,0), and incorrect vertical mean CF(i,0).
    // Notice that barotropic component is replaced with its fast-time averaged values.
    // DC = integrated height, CF integrated flux.
    CF(0)=0.0;
    FC(0)=0.0;
    DCsum = 0.0;
    for (K=0; K<N; K++)
    {
         DC = on_u(0,0)*RtoU(Hz(0,0,0));
         DCsum += DC(0,0,0);
         CF(0) = CF(i, + DC(0,0,0)*u(0,0,0);
    }


    invDCsum =1.0/DCsum;                            // recursive
    CF(0) = invDCsum*(CF(0) - DU_avg1(0,0));        // recursive

    ubar(i,j,1) = invDCsum*DU_avg1(0,0);
    ubar(i,j,2) = ubar(i,j,1);

    // Replace only BOUNDARY POINTS incorrect vertical mean with more accurate barotropic component, ubar = DU_avg1/(D*on_u).
    // Recall that, D=CF(:,0).
    //



//        IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
//         IF (DOMAIN(ng)%Western_Edge(tile)) THEN
//           DO k=1,N(ng)
//             u(Istr,j,k,nnew)=u(Istr,j,k,nnew)-CF(Istr,0)
//           END DO
//         END IF
//       END IF
///         IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
//           DO k=1,N(ng)
//             u(Iend+1,j,k,nnew)=u(Iend+1,j,k,nnew)-CF(Iend+1,0)
//           END DO
//         END IF
//
//         IF (j.eq.0) THEN                           ! southern boundary
//           DO k=1,N(ng)                             ! J-loop pipelined
//             DO i=IstrU,Iend
//               u(i,j,k,nnew)=u(i,j,k,nnew)-CF(i,0)
//             END DO
//           END DO
//         END IF
//
//         IF (j.eq.Mm(ng)+1) THEN                    # northern boundary
//           DO k=1,N(ng)                             # J-loop pipelined
//             DO i=IstrU,Iend
//               u(i,j,k,nnew) = u(i,j,k,nnew) - CF(i,0)
//             END DO
//           END DO
//         END IF
//
//#         # Compute correct mass flux, Hz*u/n.
//#
//#         DO k=N(ng),1,-1
//#           DO i=IstrP,IendT
//#             Huon(i,j,k) = 0.5*(Huon(i,j,k) + u(i,j,k,nnew)*DC(i,k))
//#           END DO
//#         END DO
//#
//#         DO i=IstrP,IendT
//#           FC(i,0) = DC(i,0)*(FC(i,0) - DU_avg2(i,j))        # recursive
//#         END DO
//#
//#         DO k=1,N(ng)
//#           DO i=IstrP,IendT
//#             Huon(i,j,k) = Huon(i,j,k) - DC(i,k)*FC(i,0)
//#           END DO
//#         END DO
//#
//#
//#
//#
//#
//#     # Couple velocity component in the ETA-direction.
//#
//#         IF (j.ge.Jstr) THEN
//#           DO i=IstrT,IendT
//#             DC(i,0)=0.0
//#             CF(i,0)=0.0
//#
//#             FC(i,0)=0.0
//#           END DO
//

}

template<typename RtoU, typename isUnode>
void correctBaroclinicMeanVel(const double *_u, const double *_Hz)
{
    // Replace INTERIOR POINTS incorrect vertical mean with more accurate barotropic component, vbar=DV_avg1/(D*om_v).
    // Recall that, D=CF(:,0).
    //
    // NOTE:  Only the BOUNDARY POINTS need to be replaced. Avoid redundant update in the interior again for computational purposes
    //        which will not affect the nonlinear code.
v(i,j,k,nnew)
om_v
    if !isUnode(i) return;

    STENCIL3D(u,   K);
    STENCIL3D(Hz,  K);

    Hzk = RtoU(Hz)

    double H = 0.0;
    double Hu = 0.0;

    for (K=0; K<N; K++)
    {
        H += Hz(0,0,0)
        Hu += u(0,0,0)*RtoU(Hz)
    }

    const double udiff = (Hu*on_u(0,0) - DU_avg1(0,0))/(H*on_u(0,0));    // recursive

    // Couple and update new solution.
    for (K=0; K<N; K++)
    {
        u -= udiff;
    }

}


void solveTri()
// Solve tridiagonal system.
// -------------------------
{
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
//# Solve the tridiagonal system.
//
//        # DO k=1,N(ng)
//        #   DO i=IstrU,Iend
//        #     DC(i,k)=u(i,j,k,nnew)
//        #     BC(i,k)= Hzk(i,k) - FC(i,k) - FC(i,k-1)
//        #   END DO
//        # END DO
//        #
//        # DO i=IstrU,Iend
//        #   cff=1.0/BC(i,1)
//        #   CF(i,1)=cff*FC(i,1)
//        #   DC(i,1)=cff*DC(i,1)
//        # END DO
//        #
//        # DO k=2,N(ng)-1
//        #   DO i=IstrU,Iend
//        #     cff=1.0_r8/(BC(i,k)-FC(i,k-1)*CF(i,k-1))
//        #     CF(i,k)=cff*FC(i,k)
//        #     DC(i,k)=cff*(DC(i,k)-FC(i,k-1)*DC(i,k-1))
//        #   END DO
//        # END DO
//        #
//        # # Compute new solution by back substitution.
//        #
//        # DO i=IstrU,Iend
//        #
//        #   DC(i,N(ng))=(DC(i,N(ng))-FC(i,N(ng)-1)*DC(i,N(ng)-1))/(BC(i,N(ng))-FC(i,N(ng)-1)*CF(i,N(ng)-1))
//        #   u(i,j,N(ng),nnew) = DC(i,N(ng))
//        #
//        # END DO
//        #
//        #
//        # DO k=N(ng)-1,1,-1
//        #   DO i=IstrU,Iend
//        #
//        #     DC(i,k) = DC(i,k) - CF(i,k)*DC(i,k+1)
//        #     u(i,j,k,nnew) = DC(i,k)
//        #
//        #   END DO
//        # END DO
}



template<typename RtoU, typename isUnode>
void createVertViscousOpMatrix(int &K, const double cff, const double Δt, const double lambda,
                               const double *_z_r, const double *_u, const double *_ru)
{
    //    u  = u(nnew,:,:,:)
    //    ru = ru(nrhs,:,:,:)
    STENCIL3D(u,   K);
    STENCIL3D(ru,  K);
    STENCIL3D(Akv, K);
    STENCIL3D(z_r, K);
    STENCIL3D(Hz,  K);

    double *FC = bufFC[i*N];

    // If the index is not a U-node leaves.
    if (!isUnode(i)) return;



    // Builds matrix M and vector RHS.
    //-------------------------------
    const auto AKvU = RtoU(Akv);
    const auto HzU  = RtoU(Hz);

    const auto pmnU = RtoU(pm)*RtoU(pn)
    const auto cΔt_mn = coef*Δt*pmnU;

    FC[0] = 0.0;
    FC[N] = 0.0;
    for (K=1; K < N; K++)
    {
        // Δz is the vertical distance between two U nodes.
        auto Δz = RtoU(z_r(1,0,0) - z_r(0,0,0));

        RHS = u(0,0,0) + cΔt_mn*ru(0,0,0);

        // Off-diagonal elements
        FC = -lambda*Δt*AKvU(0,0,0)/Δz;

        //Diagonal elements.
        BC = HzU(0,0,0) - FC(0,0,0) - FC(-1,0,0);
    }
    BC[N] = HzU(0,0,0) - FC(0,0,0) - FC(0,0,-1);


    // Now RHS contains the RHS, BC is the main diagonal and FC[1:-1], FC[:,-2] the other two diagonals.


}




// This subroutine time-steps the nonlinear  horizontal  momentum equations.
// The vertical viscosity terms are time-stepped using an implicit algorithm.

void step3d_UV()
{

    u(i,j,k, nnew);
    ru(i,j,k, nrhs);



    // Time step momentum equations
    // ----------------------------



    # IF (iic(ng).eq.ntfirst(ng)) THEN
    #   cff = 1
    # ELSE IF (iic(ng).eq.(ntfirst(ng)+1)) THEN
    #   cff = 3.0/2.0
    # ELSE
    #   cff = 23.0/12.0
    # END IF


    // ξ-direction.
    createVertViscousOpMatrix<RtoU, isUnode>(cff, lambda, MuU);

    u = solveTri(MuU, AK, z_r, u);

    correctBaroclinicMeanVel(u);



    // η-direction.
    createVertViscousOpMatrix<RtoU, isUnode>(cff, lambda, MvV)

    v = solveTri(MvV, AK, z_r, v)

    correctBaroclinicMeanVel(v)



    setUVBCs(u,v)
//    # # Set lateral boundary conditions.
//    # # -----------------------------------------------------------------------
//    #
//    #   CALL u3dbc_tile (LBi, UBi, LBj, UBj, N(ng),                       &
//    #  &                 IminS, ImaxS, JminS, JmaxS,                      &
//    #  &                 nstp, nnew,                                      &
//    #  &                 u)
//    #   CALL v3dbc_tile (LBi, UBi, LBj, UBj, N(ng),                       &
//    #  &                 IminS, ImaxS, JminS, JmaxS,                      &
//    #  &                 nstp, nnew,                                      &
//    #  &                 v)
//    #
//    #
//    #

    applyPointSources()



    // Couple 2D and 3D momentum equations.
    // -----------------------------------------------------------------------
    correctBarotropicVel<UtoR>(ubar);
    correctBarotropicVel<VtoR>(vbar);


}