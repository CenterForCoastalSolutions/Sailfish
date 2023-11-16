

// The Eq. are:
// $$\frac{\partial\boldsymbol{u}}{\partial t}+\frac{\partial}{\partial z}K_{v}\frac{\partial\boldsymbol{u}}{\partial z}=\boldsymbol{ru}$$
// Which can be discretized as (for each component)
// $$\frac{u_{k}^{n+1}-u_{k}^{n}}{\Delta t}+\frac{1}{\Delta z_{k}}\left(K_{v}^{^{k-\nicefrac{1}{2}}}\frac{u_{k}^{n+1}-u_{k-1}^{n+1}}{\Delta z_{k-\nicefrac{1}{2}}}-K_{v}^{^{k+\nicefrac{1}{2}}}\frac{u_{k+1}^{n+1}-u_{k}^{n+1}}{\Delta z_{k+\nicefrac{1}{2}}}\right)=ru$$
// After some algebra:
// $$\underbrace{\left(\Delta z_{k}+\underbrace{\frac{\Delta tK_{v}^{^{k-\nicefrac{1}{2}}}}{\Delta z_{k-\nicefrac{1}{2}}}}_{-FC[k-1]}+\underbrace{\frac{\Delta tK_{v}^{^{k+\nicefrac{1}{2}}}}{\Delta z_{k+\nicefrac{1}{2}}}}_{-FC[k]}\right)}_{BC[k]}u_{k}^{n+1}\underbrace{-\Delta t\frac{K_{v}^{^{k-\nicefrac{1}{2}}}}{\Delta z_{k-\nicefrac{1}{2}}}}_{FC[k-1]}u_{k-1}^{n+1}\underbrace{-\Delta t\frac{K_{v}^{^{k+\nicefrac{1}{2}}}}{\Delta z_{k+\nicefrac{1}{2}}}}_{FC[k]}u_{k+1}^{n+1}	=\underbrace{\Delta z_{k}\left(\Delta t\,ru+u^{n}\right)}_{RHS[k]}$$
// The matrix equation is $M\,U=RHS$, where:
// $$M=\left(\begin{array}{ccccccc}BC[1] & FC[1]\\FC[1] & BC[2] & FC[2]\\ & FC[2] & BC[3] & FC[3]\\ &  & \ddots & \ddots & \ddots\\ &  &  & \ddots & \ddots & \ddots\\ &  &  &  & FC[N-2] & BC[N-1] & FC[N-1]\\ &  &  &  &  & FC[N-1] & BC[N]\end{array}\right)$$
// In reality, the terms with η or ξ derivatives are divided bi m and n respectively, while the terms with no derivatives are divided by m*n.



#include "mod_cppkernels.h"



VerticalVelEq *velEq;


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



extern "C"  __global__
void omega(const double *_W, const double *_u, const double *_v, const double *_Huon, const double *_Hvom, const double *_z_w, const int *idxFieldBC)
// This routine computes S-coordinate vertical velocity (m^3/s),
//
//                 W = [Hz/(m*n)]*omega,
//
// diagnostically at horizontal RHO-points and vertical W-points.
//
// Added implicit vertical advection from an adaptive, Courant-number-dependent implicit scheme for vertical advection in oceanic modeling, Alexander F. Shchepetkin, pp 38-69.
{
    const unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;

//    if (((i % szI) == 0 || ((i % szI) >= (szI - 2) || (i/szI) == 0) || (i/szI) >= (szJ - 2)) || i >= sz2D) return;

    if (i >= sz2D) return;

    const int N = szK;
    int K = 0;

    STENCIL3D(W,  K);

    if (idxFieldBC[i] >= 0)
    // This is a homogenous gradient BC
    {
        W = _W[idxFieldBC[i]]; // Makes W at the BC equal to W at cell BC[i], which is the closest one normal to the BC.
        return;
    }


    STENCIL3D(u,    K);
    STENCIL3D(v,    K);
    STENCIL3D(Huon, K);
    STENCIL3D(Hvom, K);
    STENCIL3D(z_w,  K);



    //Vertically integrate horizontal mass flux divergence.
    // Starting with zero vertical velocity at the bottom, integrate from the bottom (k=0) to the free-surface (k=N).  The w(:,:,N)
    // contains the vertical velocity at the free-surface, d(zeta)/d(t).
    // Notice that barotropic mass flux divergence is not used directly.
    W[0] = 0.0;

    for (K = 1; K < N; K++)
    {
        W = W(-1,0,0) - (DξUtoR(Huon) + DηVtoR(Hvom));
    }

// TODO : Finish this.
//        #  Apply mass point sources (volume vertical influx), if any.
//        #
//        #  Overwrite W(Isrc,Jsrc,k) with the same divergence of Huon,Hvom as
//        #  above but add in point source Qsrc(k) and reaccumulate the vertical
//        #  sum to obtain the correct net Qbar given in user input - J. Levin
//        #  (Jupiter Intelligence Inc.) and J. Wilkin
//        #
//        #    Dsrc(is) = 2,  flow across grid cell w-face (positive or negative)
//
//        IF (LwSrc(ng)) THEN
//          DO is=1,Nsrc(ng)
//            IF (INT(SOURCES(ng)%Dsrc(is)).eq.2) THEN
//              ii=SOURCES(ng)%Isrc(is)
//              jj=SOURCES(ng)%Jsrc(is)
//              IF (((IstrR.le.ii).and.(ii.le.IendR)).and.                &
//     &            ((JstrR.le.jj).and.(jj.le.JendR)).and.                &
//     &            (j.eq.jj)) THEN
//
//                DO k=1,N(ng)
//                  W(ii,jj,k)=W(ii,jj,k-1)-                              &
//
//
//     &                       (Huon(ii+1,jj,k)-Huon(ii,jj,k)+            &
//     &                        Hvom(ii,jj+1,k)-Hvom(ii,jj,k))+           &
//     &                       SOURCES(ng)%Qsrc(is,k)
//                END DO
//              END IF
//            END IF
//          END DO
//        END IF
//
//        DO i=Istr,Iend
//          wrk(i)=W(i,j,N(ng))/(z_w(i,j,N(ng))-z_w(i,j,0))
//
//        END DO


    // In order to ensure zero vertical velocity at the free-surface, subtract the vertical velocities of the moving S-coordinates
    // isosurfaces. These isosurfaces are proportional to d(zeta)/d(t).
    // The proportionally coefficients are a linear function of the S-coordinate with zero value at the bottom (K=0) and unity at
    // the free-surface (K=N-1).

    K = 0;
    auto z_w_bed = z_w.Eval(0,0,0);
    auto h = z_w.Eval(0,0,N-1) - z_w_bed;
    auto Ws = W.Eval(0,0,N-1);  // Vertical velocity at the surface
    for (K = N-1; K > 0; K--)
    {
        W -= Ws*(z_w(0,0,0) - z_w_bed)/h;
    }
    W[N - 1] = 0.0;

}


//extern "C"  __global__
//void wvelocity(const double *_u, const double *_v, const double *_Huon, const double *_Hvom)
//// This routine computes horizontal mass fluxes, Hz*u/n and Hz*v/m.
//{
////  Compute "true" vertical velocity (m/s).
////
////      In ROMS, the terrain-following vertical velocity, omega, is given by:
////
////             Hz * omega = w - dz/dt - div(z)
////
////      where w is the "true" vertical velocity and
////
////             div(z) = pm*u*(dz/dξ) + pn*v*(dz/dη)
////
////      The vertical coordinate is a function of several parameters but only
////      the free-surface is time dependent. However, in sediment applications
////      with stratigraphy, the bathymetry (h) also evolves in time.
//
//    const int N = szK;
//    int K = 0;
//
//    // Compute contribution due to quasi-horizontal motions along S-coordinate surfaces:  U·GRADs(z).
////    auto vert = UtoR(u[Ninp,:,:,:]*dξRtoU(z_r)) + VtoR(v[Ninp,:,:,:]*DξRtoU(z_r))
//
//
//    // Compute contribution due to time tendency of the free-surface, d(zeta)/d(t), which is the vertical velocity at the free-surface
//    // and it is expressed in terms of barotropic mass flux divergence.
//    // Notice that it is divided by the total depth of the water column.
//    // This is needed because this contribution is linearly distributed throughout the water column by multiplying it by the distance from
//    // the bottom to the depth at which the vertical velocity is computed.
//
//    double invWcDepth = 1.0/(z_w(0,0,N-1) - z_w(0,0,0));   // Water Column (wc) depth (2D), remember that K=0;
//    double z_w_bed = z_w.Eval(0,0,0);                      // Location of the bed.
//    auto   lcDepth = z_w(0,0,0) - z_w_bed;                 // Local (lc) depth (3D). *NOTE* z_w(0,0,0) means the local value of s_w at (i,j,K) calculated lazily later
//
//    auto vertW = RtoW(UtoR(u*dξRtoU(z_r)) + VtoR(v*dηRtoV(z_r)));
//
////    auto wvel = pm*pn*(W - ηξdivR(DU_avg1)*lcDepth/wcDepth) + RtoW(vert);
////
////    krnl(u, v, z_r, pm, pn);
//
//
//
//    // Compute contribution due to time tendency of the free-surface, d(zeta)/d(t), which is the vertical velocity at the free-surface
//    // and it is expressed in terms of barotropic mass flux divergence.
//    // Notice that it is divided by the total depth of the water column.
//    // This is needed because this contribution is linearly distributed throughout the water column by multiplying it by the distance from
//    // the bottom to the depth at which the vertical velocity is computed.
//
//    for (K=0; K<N; K++)
//    {
//        wvel = pm*pn*(W - UVdivR(DU_avg1)*lcDepth*invWcDepth) + vertW;
//    }
//
//}


extern "C"  __global__
void set_depth(const int Vtransform, const double *_Zt_avg1, const double *_z_w, const double *_z_r, const double *_h, const double hc, const double *_Hz,
               const double *sc_r,  const double *sc_w, const double *Cs_r,  const double *Cs_w)
//This routine computes the time evolving depths of the model grid and its associated vertical transformation metric (thickness).      !
//Currently, two vertical coordinate transformations are available with various possible vertical stretching, C(s), functions, (see
//routine "set_scoord.F" for details).
{
    const unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
    const int N = szK;
    int K = 0;

    if (!isRNode(i)) return;

    STENCIL(Zt_avg1);
    STENCIL(z_w);
    STENCIL(z_r);
    STENCIL(h);
    STENCIL3D(Hz, K);

    if (Vtransform == 1)
    //Original formulation: Compute vertical depths (meters, negative) at RHO- and W-points, and vertical grid
    //thicknesses. Various stretching functions are possible.
    //        z_w(x,y,s,t)  =  Zo_w  +  zeta(x,y,t) * [1.0 + Zo_w / h(x,y)]
    //                Zo_w  =  hc * [s(k) - C(k)]  +  C(k) * h(x,y)
    {
        z_w = -h;

        for (K=1; K<=N; K++) // TODO: Is it K=1..N or 1..N-1?
        {
            double cff_r  = hc*(sc_r[K] - Cs_r[K]);
            double cff_w  = hc*(sc_w[K] - Cs_w[K]);
            double cff1_r = Cs_r[K];
            double cff1_w = Cs_w[K];

            auto z_w0 = cff_w + cff1_w*h(0,0);
            auto z_r0 = cff_r + cff1_r*h(0,0);

            z_w = z_w0 + Zt_avg1(0,0)*(1.0 + z_w0/h(0,0));
            z_r = z_r0 + Zt_avg1(0,0)*(1.0 + z_r0/h(0,0));

            Hz = DσWtoR(z_w);
        }
    }
    else if (Vtransform == 2)
    // New formulation: Compute vertical depths (meters, negative) at RHO -  and W - points, and vertical grid thicknesses.
    // Various stretching functions are possible.
    //       z_w(x,y,s,t)  =  zeta(x,y,t)  +  [zeta(x,y,t) +  h(x,y)] * Zo_w
    {
        z_w = -h;

        for (K=1; K<=N; K++)
        {
            double cff_r = hc*sc_r[K];
            double cff_w = hc*sc_w[K];

            double cff1_r = Cs_r[K];
            double cff1_w = Cs_w[K];

            auto cff2_r = (cff_r + cff1_r*h)/(hc + h);
            auto cff2_w = (cff_w + cff1_w*h)/(hc + h);

            z_w = Zt_avg1 + (Zt_avg1 + h)*cff2_w;
            z_r = Zt_avg1 + (Zt_avg1 + h)*cff2_r;


            Hz = DσWtoR(z_w);
        }

    }

}


extern "C"  __global__
void set_zeta(const double *_zeta1, const double *_zeta2, const double *_Zt_avg1)
// This routine sets free-surface to its fast-time averaged value.
{
    const unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (!isRNode(i)) return;

    STENCILR(zeta2);
    STENCILR(zeta1);
    STENCILR(Zt_avg1);

    // Prepare to time-step 2D equations:  set initial free-surface to its fast-time averaged values (which corresponds
    // to the time step "n").
    zeta2 = Zt_avg1;
    zeta1 = Zt_avg1;
}



extern "C"  __global__
void setLateralUVBCs(double *u, double *v, const int *bcUIdxFieldIdx2, const int *bcVIdxFieldIdx2, const int *bcUIdxFieldType, const int *bcVIdxFieldType)
{
    const unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
    const int N = szK;

    if (i >= sz2D) return;

    const int bcUType = bcUIdxFieldType[i];
    const int bcVType = bcVIdxFieldType[i];

    if ((bcUType == bcNone) && (bcVType == bcNone)) return;


    switch (bcUType)
    {

        case bcClosed:
            for (int k = 0; k < N; k++) u[i + k*sz2D] = 0.0;

        break;

        case bcGradient:
            for (int k = 0; k < N; k++) u[i + k*sz2D] = u[bcUIdxFieldIdx2[i + k*sz2D]];

        break;
    }

    switch (bcVType)
    {

        case bcClosed:
            for (int k = 0; k < N; k++) v[i + k*sz2D] = 0.0;

        break;

        case bcGradient:
            for (int k = 0; k < N; k++) v[i + k*sz2D] = v[bcVIdxFieldIdx2[i + k*sz2D]];

        break;
    }


}



extern "C"  __global__
void set_maxflux(const double *_u, const double *_v, const double *_Huon, const double *_Hvom, const double *_Hz)
// This routine computes horizontal mass fluxes, Hz*u/n and Hz*v/m.
{
    const unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;

    int K = 0;  // All stencils share this variable.
    const int N = szK;

    STENCIL3D(u,    K);
    STENCIL3D(v,    K);
    STENCIL3D(Huon, K);
    STENCIL3D(Hvom, K);
    STENCIL3D(on_u, K);
    STENCIL3D(om_v, K);
    STENCIL3D(Hz,   K);

    if (isUNode(i))
    {
        for (K = 0; K < N; K++)
        {
            Huon = RtoU(Hz)*u*on_u;
        }
    }

    if (isVNode(i))
    {
        for (K = 0; K < N; K++)
        {
            Hvom = RtoV(Hz)*v*om_v;
        }
    }

}

extern "C"  __global__
void adjustBarotropicVelocity(const double *_Hz, const double *_u, const double *_v,
                              const double *_DU_avg1, const double *_DV_avg1)
// Replace INTERIOR POINTS? incorrect vertical mean velocity with more accurate barotropic component, ubar=DU_avg1/(D*om_u)
// and vbar=DV_avg1/(D*om_v).
//
//      _HzU, _HzV: Hz interpolated at U and V points.
//      _u, _v    : velocity components.
{
    const unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
    const int N = szK;
    int K = 0;   // This variable is used by all stencils to define the level they act upon.

    STENCIL3D(Hz,      K);
    STENCIL3D(u,       K);
    STENCIL3D(v,       K);
    STENCIL(on_u);
    STENCIL(om_v);
    STENCIL(DU_avg1);
    STENCIL(DV_avg1);

    auto HzU = RtoU(Hz);
    auto HzV = RtoV(Hz);

    double HU  = 0.0, HV  = 0.0;
    double HuU = 0.0, HvV = 0.0;

    for (K = 0; K < N; K++)
    {
        HU  += (HzU  ).Eval(0,0,0);
        HV  += (HzV  ).Eval(0,0,0);

        HuU += (HzU*u).Eval(0,0,0);
        HvV += (HzV*v).Eval(0,0,0);
    }

    auto uUpdate = (HuU - DU_avg1/HU)*on_u;
    auto vUpdate = (HvV - DV_avg1/HV)*om_v;


    // Update new solution.
    for (K=0; K<N; K++)
    {
        u -= uUpdate;
        v -= vUpdate;
    }
}


//template<typename RtoU, typename isUNode>

template<decltype(isUNode) isUNode, template<typename> class RtoU, typename T>
void correctBaroclinicMeanVel(const double *_u, const double *_Hz, const double *_DU_avg1, double *_on_u)
{
    // Replace INTERIOR POINTS incorrect vertical mean with more accurate barotropic component, vbar=DV_avg1/(D*om_v).
    // Recall that, D=CF(:,0).
    //
    // NOTE:  Only the BOUNDARY POINTS need to be replaced. Avoid redundant update in the interior again for computational purposes
    //        which will not affect the nonlinear code.
    const unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
    const int N = szK;
    int K = 0;

    if (!isUNode(i)) return;

    STENCIL3D(u,   K);
    STENCIL3D(Hz,  K);
    STENCIL(on_u);
    STENCIL(DU_avg1);

    auto HzU = RtoU<T>(Hz);

    double H = 0.0;
    double Hu = 0.0;

    for (K=0; K<N; K++)
    {
        H += HzU.Eval(0,0,0);
        Hu += (u*HzU).Eval(0,0,0);
    }

    const double udiff = ((Hu*on_u - DU_avg1)/(H*on_u)).Eval(0,0);    // recursive

    // Couple and update new solution.
    for (K=0; K<N; K++)
    {
        u -= udiff;
    }

}


double solveTri(const void *velEq, double AK, const double *z_r, const double *u)
// Solve tridiagonal system.
// -------------------------
{
    // See https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm.


//    // Forward substitution.
//    CF[1] = FC[1]/BC[1];
//    DC[1] = DC[1]/BC[1];
//
//    for (int k=2; k<=N; k++)
//    {
//        double denom = 1.0/(BC[k] - FC[k-1]*CF[k-1]);
//
//        CF(i,k) = denom*FC[k];
//        DC(i,k) = denom*(DC[k] - FC[k-1]*DC[k-1]);
//    }
//
//
//    // Back substitution.
//    DC[N] -= FC[N-1]*DC[N-1]/(BC[N] - FC[N-1]*CF[N-1]);
//    u(0,0,N,nnew) = DC[N];
//
//    for (int k=N-1; k>1; k--)
//    {
//        u(0,0,k,nnew) = DC[k] - CF[k]*u(0,0,k+1,nnew);
//    }
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
    return 0.0;
}


template<NodeType nt, decltype(isUNode) isNode>
void createVertViscousOpMatrix(int &K, double cff, double Δt, double lambda, VerticalVelEq *velEq,
                               const double *_Hz, const double *_Akv, const double *_z_r, const double *_u, const double *_ru)
{
    const double *_FC  = velEq->SD;
    const double *_BC  = velEq->D;
    const double *_RHS = velEq->RHS;

    const unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
    const int N = szK;

    STENCIL3D(u,   K);
    STENCIL3D(ru,  K);
    STENCILR3D(Akv, K);
    STENCILR3D(z_r, K);
    STENCILR3D(Hz,  K);
    STENCILR3D(FC,  K);
    STENCILR3D(BC,  K);
    STENCILR3D(RHS, K);
    STENCILR(pm);
    STENCILR(pn);



    // If the index is not a U or V-node, leaves.
    if (!isNode(i)) return;

    // Builds matrix M and vector RHS.
    //-------------------------------
    const auto AKvU = Akv.to<nt>();  // to<nt>() interpolates Akv from R node to the node defined by nt (either U or V).
    const auto HzU  = Hz.to<nt>();

    const auto pmU  = pm.to<nt>();
    const auto pnU  = pn.to<nt>();

    const auto pmnU = pmU*pnU;
    const auto cΔt_mn = cff*Δt*pmnU;



    FC[0] = 0.0;
    FC[N] = 0.0;
    for (K=1; K < N; K++)
    {
        // Δz is the vertical distance between two U nodes.
        auto Δz = DσR(z_r).to<nt>();

        RHS = u(0,0,0) + cΔt_mn*ru(0,0,0);

        // Off-diagonal elements
        FC = -lambda*Δt*AKvU(0,0,0)/Δz;

        //Diagonal elements.
        BC = HzU(0,0,0) - FC(0,0,0) - FC(-1,0,0);
    }
    BC[N] = HzU(0,0,0) - FC(0,0,0) - FC(0,0,-1);

    // At this point, RHS contains the RHS, BC is the main diagonal and FC[1:-1], FC[:,-2] the other two diagonals.

}




// This subroutine time-steps the nonlinear  horizontal  momentum equations.
// The vertical viscosity terms are time-stepped using an implicit algorithm.
extern "C"  __global__
void step3d_UV(const double *_u, const double *_v, const double *ru, const double *rv, const double *ubar_t2, const double *vbar_t2,
               const double *Hz, const double *Akv, const double *z_r, const double *DU_avg1, const double *DV_avg1,
               int iic, int ntfirst, double lambda, double AK, double Dt)
{
    const unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
//    const int N = szK;

    // Time step momentum equations
    // ----------------------------

    double cff;
    if (iic == ntfirst)           cff = 1;
    else if (iic == ntfirst + 1)  cff = 3.0/2.0;
    else                          cff = 23.0/12.0;

    int K = 0;

    STENCIL3D(u, K);
    STENCIL3D(v, K);

//    for (K=0 ; K<N; K++)
//    {

    // ξ-direction.
    createVertViscousOpMatrix<ntU, isUNode>(K, cff, Δt, lambda, velEq, Hz, Akv, z_r, _u, ru);

    u = solveTri(velEq, AK, z_r, _u);

//    correctBaroclinicMeanVel<decltype(u), isUNode, RtoU<decltype(u)> >(_u, Hz, DU_avg1, _on_u);



    // η-direction.
//    createVertViscousOpMatrix<RtoV<double>, isVNode>(K, cff, Δt, lambda, velEq, Hz, Akv, z_r, _v, rv);

    v = solveTri(velEq, AK, z_r, _v);

//    correctBaroclinicMeanVel<isVNode, RtoV>(_v, Hz, DV_avg1, _on_u);

//    }

    applyPointSources();


    // Couple 2D and 3D momentum equations.
    // -----------------------------------------------------------------------
    adjustBarotropicVelocity(Hz, _u, _v, DU_avg1, DV_avg1);


}