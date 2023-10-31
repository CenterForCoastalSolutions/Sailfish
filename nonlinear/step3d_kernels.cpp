#include "mod_cppkernels.h"


extern "C"  __global__
void omega(const double *_W, const double *_u, const double *_v, const double *_Huon, const double *_Hvom, const double *_z_w)
// This routine computes S-coordinate vertical velocity (m^3/s),
//
//                 W = [Hz/(m*n)]*omega,
//
// diagnostically at horizontal RHO-points and vertical W-points.
//
// Added implicit vertical advection from an adaptive, Courant-number-dependent implicit scheme for vertical advection in oceanic modeling, Alexander F. Shchepetkin, pp 38-69.

    STENCIL3D(W);
    STENCIL3D(u);
    STENCIL3D(v);
    STENCIL3D(Huon);
    STENCIL3D(Hvom);
    STENCIL3D(z_w);

    //Vertically integrate horizontal mass flux divergence.
    // Starting with zero vertical velocity at the bottom, integrate from the bottom (k=0) to the free-surface (k=N).  The w(:,:,N)
    // contains the vertical velocity at the free-surface, d(zeta)/d(t).
    // Notice that barotropic mass flux divergence is not used directly.
    int K = 0;
    W = 0.0;

    for (K = 1; K < N; K++)
    {
        W = W(-1,0,0) - (Dξ(Huon) + Dη(Hvom));
    }

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
    for (K = N-1; K >0; K--)
    {
        W -= Ws*(z_w(0,0,0) - z_w_bed)/h;
    }
    K = N - 1;
    W = 0.0;

}



extern "C"  __global__
void set_maxflux(const double *_u, const double *_v, const double *_Huon, const double *_Hvom)
// This routine computes horizontal mass fluxes, Hz*u/n and Hz*v/m.
{
//  Compute "true" vertical velocity (m/s).
//
//      In ROMS, the terrain-following vertical velocity, omega, is given by:
//
//             Hz * omega = w - dz/dt - div(z)
//
//      where w is the "true" vertical velocity and
//
//             div(z) = pm*u*(dz/dξ) + pn*v*(dz/dη)
//
//      The vertical coordinate is a function of several parameters but only
//      the free-surface is time dependent. However, in sediment applications
//      with stratigraphy, the bathymetry (h) also evolves in time.

    int K = 0;

    // Compute contribution due to quasi-horizontal motions along S-coordinate surfaces:  U·GRADs(z).
    vert = UtoR(u[Ninp,:,:,:]*dξRtoU(z_r)) + VtoR(v[Ninp,:,:,:]*DξRtoU(z_r))


    // Compute contribution due to time tendency of the free-surface, d(zeta)/d(t), which is the vertical velocity at the free-surface
    // and it is expressed in terms of barotropic mass flux divergence.
    // Notice that it is divided by the total depth of the water column.
    // This is needed because this contribution is linearly distributed throughout the water column by multiplying it by the distance from
    // the bottom to the depth at which the vertical velocity is computed.

    double invWcDepth = 1.0/(z_w(0,0,N-1) - z_w(0,0,0));   // Water Column (wc) depth (2D), remember that K=0;
    double z_w_bed = z_w.Eval(0,0,0);                      // Location of the bed.
    auto   lcDepth = z_w(0,0,0) - z_w_bed;                 // Local (lc) depth (3D). *NOTE* z_w(0,0,0) means the local value of s_w at (i,j,K) calculated lazily later

    auto vertW = RtoW(UtoR(u*dξRtoU(z_r)) + VtoR(v*dηRtoV(z_r)));

//    auto wvel = pm*pn*(W - ηξdivR(DU_avg1)*lcDepth/wcDepth) + RtoW(vert);
//
//    krnl(u, v, z_r, pm, pn);



    // Compute contribution due to time tendency of the free-surface, d(zeta)/d(t), which is the vertical velocity at the free-surface
    // and it is expressed in terms of barotropic mass flux divergence.
    // Notice that it is divided by the total depth of the water column.
    // This is needed because this contribution is linearly distributed throughout the water column by multiplying it by the distance from
    // the bottom to the depth at which the vertical velocity is computed.

    for (K=0; K<N; K++)
    {
        wvel = pm*pn*(W - UVdivR(DU_avg1)*lcDepth*invWcDepth) + vertW;
    }

}


extern "C"  __global__
void set_depth(const int Vtransform, const double *_z)
//This routine computes the time evolving depths of the model grid and its associated vertical transformation metric (thickness).      !
//Currently, two vertical coordinate transformations are available with various possible vertical stretching, C(s), functions, (see
//routine "set_scoord.F" for details).
{


    if !isTnode(i) return;

    STENCIL(Zt_avg1);

    if (Vtransform == 1)
    //Original formulation: Compute vertical depths (meters, negative) at RHO- and W-points, and vertical grid
    //thicknesses. Various stretching functions are possible.
    //        z_w(x,y,s,t)  =  Zo_w  +  zeta(x,y,t) * [1.0 + Zo_w / h(x,y)]
    //                Zo_w  =  hc * [s(k) - C(k)]  +  C(k) * h(x,y)
    {
        K = 0;
        z_w(0,0,0) = -h;

        for (K=1; K<=N; K++) // TODO: Is it K=1..N or 0..N-1
        {
            double cff_r  = hc*(SCALARS.sc_r[k] - SCALARS.Cs_r[k]);
            double cff_w  = hc*(SCALARS.sc_w[k] - SCALARS.Cs_w[k]);
            double cff1_r = SCALARS.Cs_r[k];
            double cff1_w = SCALARS.Cs_w[k];

            auto z_w0 = cff_w + cff1_w*h(0,0);
            auto z_r0 = cff_r + cff1_r*h(0,0);

            z_w = z_w0 + Zt_avg1(0,0)*(1.0 + z_w0/h(0,0));
            z_r = z_r0 + Zt_avg1(0,0)*(1.0 + z_r0/h(0,0));

            Hz = z_w(0,0,0) - z_w(0,0, -1);
        }
    }
    else if (Vtransform == 2)
    // New formulation: Compute vertical depths (meters, negative) at RHO -  and W - points, and vertical grid thicknesses.
    // Various stretching functions are possible.
    //       z_w(x,y,s,t)  =  zeta(x,y,t)  +  [zeta(x,y,t) +  h(x,y)] * Zo_w
    {
        K = 0;
        z_w = - h(0,0);

        for (K=1; K<=N; K++)
        {
            double cff_r = hc*SCALARS.sc_r[k];
            double cff_w = hc*SCALARS.sc_w[k];

            double cff1_r = SCALARS.Cs_r[k];
            double cff1_w = SCALARS.Cs_w[k];

            auto cff2_r = (cff_r + cff1_r*h)/(hc + h);
            auto cff2_w = (cff_w + cff1_w*h)/(hc + h);

            z_w = Zt_avg1(0,0) + (Zt_avg1(0,0) + h(0,0))*cff2_w(0,0,0)
            z_r = Zt_avg1(0,0) + (Zt_avg1(0,0) + h(0,0))*cff2_r(0,0,0)


            Hz = z_w(0,0,0) - z_w(0,0,-1)
        }

    }

}




extern "C"  __global__
void set_maxflux(const double *_u, const double *_v, const double *_Huon, const double *_Hvom)
// This routine computes horizontal mass fluxes, Hz*u/n and Hz*v/m.
{
    const unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;

    int K = 0;  // All stencils share this variable.

    STENCIL3D(u,    K);
    STENCIL3D(v,    K);
    STENCIL3D(Huon, K);
    STENCIL3D(Hvom, K);
    STENCIL3D(on_u, K);
    STENCIL3D(om_v, K);

    if (isUnode(i))
    {
        for (K = 0; K < N; K++)
        {
            Huon = RtoU(Hz)*u*on_u;
        }
    }

    if (isVnode(i))
    {
        for (K = 0; K < N; K++)
        {
            Hvom = RtoV(Hz)*v*om_v;
        }
    }

}

extern "C"  __global__
void createVertViscousOpMatrixU(const double *_DC,  const double *_BC, const double *_FC,
                                const double *_Akv, const double *_Hz, const double *_z_r,
                                const double *_u,   const double *_ru, const int N)
{
    //The Eq. is: Un+1 = Un - Δt*∂z(Akv*∂zUn) = [1 - ∂z(Δt*Akv*∂z)]Un = Ru
    //Integrating vertically between nodes k and k+1 we have (using Gauss central quadrature for the velocity):
    //
    //Hz*Un+1 = Hz*Un + Δt*Akv*(∂zU) = Hz*Un + Akv*(Un,k+1 - Un,k)/(Zk+1 - Zk)
    u  = u(nnew,:,:,:)
    ru = ru(nrhs,:,:,:)

    int K;
    STENCIL3D(DC, K);
    STENCIL3D(BC, K);
    STENCIL3D(FC, K);

    if (isUnode(i))
    {
        auto AK  = RtoU(Akv);
        auto HzU = RtoU(Hz)
        auto zU  = RtoU(z_r);

        DC = cff*RtoU(pm)*RtoU(pn);
        for (int k=1; k<N; k++)
        {
            u += DC*ru(k,0,0);

            double Δz = zU(k+1,0,0) - zU(k,0,0);

            FC(k) = (-lambda*Δt)*AK(i,k)/Δz;
        }

        FC(0,0,0) = 0.0;
        FC(N,0,0) = 0.0;

        for (int k=1; k<=N; k++)
        {
            // RHS
            DC[k] = u(i,j,k,nnew);

            // Diagonal includes Hz and two of the components of (Un,k+1 - Un,k)/(Zk+1 - Zk) coming from the elements up and down
            BC[k] = HzU(i,j,k) - FC(k,0,0) - FC(k-1,0,0);
        }
    }

}



extern "C"  __global__
void adjustBarotropicVelocity(const double *_Hz, const double *_u, const double *_v,
                              const double *_DU_avg1, const double *_DV_avg1, const int N)
// Replace INTERIOR POINTS? incorrect vertical mean velocity with more accurate barotropic component, ubar=DU_avg1/(D*om_u)
// and vbar=DV_avg1/(D*om_v).
//
//      _HzU, _HzV: Hz interpolated at U and V points.
//      _u, _v    : velocity components.
{
    v := v(i,j,1,nnew);

    int K = 0;   // This variable is used by all stencils to define the level they act upon.

    STENCIL3D(Hz,      K);
    STENCIL3D(u,       K);
    STENCIL3D(v,       K);
    STENCIL3D(om_u,    K);
    STENCIL3D(on_v,    K);
    STENCIL3D(DU_avg1, K);
    STENCIL3D(DV_avg1, K);

    auto HzU = RtoU(Hz);
    auto HzV = RtoV(Hz);

    double HU  = 0.0, HV  = 0.0;
    double HuU = 0.0, HvV = 0.0;

    for (K = 0; k < N; k++)
    {
        HU  += (HzU  ).Eval(0,0,0);
        HV  += (HzV  ).Eval(0,0,0);

        HuU += (HzU*u).Eval(0,0,0);
        HvV += (HzV*v).Eval(0,0,0);
    }

    double uUpdate = (HuU - DU_avg1/HU)*on_u;
    double vUpdate = (HvV - DV_avg1/HV)*om_v;


    // Update new solution.
    for (K=0; K<N; K++)
    {
        u -= uUpdate;
        v -= vUpdate;
    }
}