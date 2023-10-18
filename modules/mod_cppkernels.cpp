#include "./mod_cppkernels.h"










extern "C"  __global__
void computeMomentumRHSPred(const double *_h,
                            const double *_rhs_ubar, const double *_rhs_vbar,
                            const double *_zeta_t1, const double *_zeta_t2, const double g, const double weight)
{
    const unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;

//    if (((i % szJ) == 0 || (i/szJ) == 0) || i >= sz2D) return;
    if (((i % szI) == 0 || ((i % szI) == (szI - 1) || (i/szI) == 0) || (i/szI) == (szJ - 1)) || i >= sz2D)
    {
//        rhs_ubar = 0.0;
//        rhs_vbar = 0.0;
        return;
    }

    STENCIL(rhs_ubar);
    STENCIL(rhs_vbar);
    STENCIL(h);
    STENCIL(on_u);
    STENCIL(om_v);
    STENCIL(zeta_t1);
    STENCIL(zeta_t2);

    auto gzeta  = (1 - weight)*zeta_t2 + weight*zeta_t1;
    auto gzeta2 = gzeta*gzeta;   // TODO : sqr expression.

    rhs_ubar = 0.5*g*(RtoU(h)*DERtoU(gzeta,on_u) + DERtoU(gzeta2,on_u));
    rhs_vbar = 0.5*g*(RtoV(h)*DXRtoV(gzeta,om_v) + DXRtoV(gzeta2,om_v));

}


extern "C"  __global__
void computeMomentumRHSCorr(const double *_h,
                            const double *_rhs_ubar, const double *_rhs_vbar,
                            const double *_zeta_t0, const double *_zeta_t1, const double *_zeta_t2, const double g)
{
    const unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;



//    if (((i % szJ) == 0 || (i/szJ) == 0) || i >= sz2D) return;
    if (((i % szI) == 0 || ((i % szI) == (szI - 1) || (i/szI) == 0) || (i/szI) == (szJ - 1)) || i >= sz2D)
    {
//        rhs_ubar = 0.0;
//        rhs_vbar = 0.0;
        return;
    }

    STENCIL(rhs_ubar);
    STENCIL(rhs_vbar);
    STENCIL(h);
    STENCIL(on_u);
    STENCIL(om_v);
    STENCIL(zeta_t0);
    STENCIL(zeta_t1);
    STENCIL(zeta_t2);

    constexpr double weight = 4.0/25.0;
    auto gzeta = (1 - weight)*zeta_t1 + weight*0.5*(zeta_t2 + zeta_t0);

    auto gzeta2 = gzeta*gzeta;
    rhs_ubar = 0.5*g*(RtoU(h)*DERtoU(gzeta,on_u) + DERtoU(gzeta2,on_u));
    rhs_vbar = 0.5*g*(RtoV(h)*DXRtoV(gzeta,om_v) + DXRtoV(gzeta2,om_v));
}


extern "C"  __global__
void computeZetaRHS(const double *_zeta, const double *_h, double *_ubar, const double *_vbar, double *_res)
{
    const unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
//TODO: also (i % szJ) == N || (i/szJ) == M)
//    if (((i % szJ) == 0 || (i/szJ) == 0) || i >= sz2D) return;

    if (((i % szI) == 0 || ((i % szI) == (szI - 1) || (i/szI) == 0) || (i/szI) == (szJ - 1)) || i >= sz2D)
    {

        return;
    }

    STENCIL(res);
    STENCIL(zeta);
    STENCIL(h);
    STENCIL(ubar);
    STENCIL(vbar);
    STENCIL(on_u);
    STENCIL(om_v);
    STENCIL(pn);
    STENCIL(pm);


    // Water column depth
    auto D = zeta + h;

    // Fluxes.
    auto DU = ubar*RtoU(D);
    auto DV = vbar*RtoV(D);

    res = divUVtoR(DU, DV, on_u, om_v, pn, pm);

}




extern "C"  __global__
void computeZetaPred(const double Dt, const double *_zeta_t0, const double *_zeta_t1, const double *_zeta_t2,
                     const double *_rzeta_t1)
{
    const unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (((i % szI) == 0 || (i/szI) == 0) || i >= sz2D) return;

    STENCIL(zeta_t0);
    STENCIL(zeta_t1);
    STENCIL(zeta_t2);
    STENCIL(rzeta_t1);

    zeta_t2 = zeta_t0(0,0) + 2.0*Dt*rzeta_t1(0,0);


}



extern "C"  __global__
void Pred(const double Dt, const double *_v_t1, const double *_v_t2, const double *_rhs, const double *_D_t1, const double *_D_t2)
{
    const unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (((i % szJ) == 0 || (i/szJ) == 0) || i >= sz2D) return;

    STENCIL(v_t1);
    STENCIL(v_t2);
    STENCIL(D_t1);
    STENCIL(D_t2);
    STENCIL(rhs);

    v_t2 = (v_t1(0,0)*D_t1(0,0) + Dt*rhs(0,0))/D_t2(0,0);

}

extern "C"  __global__
void computeMomentumPred(const double Dt, const double *_u_t1, const double *_u_t2, const double *_v_t1, const double *_v_t2,
                         const double *_rhsu, const double *_rhsv, const double *_h,
                         const double *_zeta_t1, const double *_zeta_t2)
{
    const unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (((i % szJ) == 0 || (i/szJ) == 0) || i >= sz2D) return;

    STENCIL(u_t1);
    STENCIL(u_t2);
    STENCIL(v_t1);
    STENCIL(v_t2);
    STENCIL(zeta_t1);
    STENCIL(zeta_t2);
    STENCIL(h);
    STENCIL(rhsu);
    STENCIL(rhsv);

    auto D_t1 = zeta_t1 + h;
    auto D_t2 = zeta_t2 + h;
    auto D_t1U = RtoU(D_t1);
    auto D_t1V = RtoV(D_t1);
    auto D_t2U = RtoU(D_t2);
    auto D_t2V = RtoV(D_t2);

    u_t2 = (u_t1(0,0)*D_t1U(0,0) + Dt*rhsu(0,0))/D_t2U(0,0);
    v_t2 = (v_t1(0,0)*D_t1V(0,0) + Dt*rhsv(0,0))/D_t2V(0,0);

}



extern "C"  __global__
void AdamsMoultonCorr3rd(const double Dt, const double *_v_t2, const double *_rhs_t0, const double *_rhs_t1, const double *_rhs_t2)
{
    const unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
    // Adams-Moulton 3rd order coefficients
    constexpr double AM3_2 =  5.0 / 12.0;
    constexpr double AM3_1 =  8.0 / 12.0;
    constexpr double AM3_0 = -1.0 / 12.0;

    if (i >= sz2D) return;

    STENCIL(rhs_t0);
    STENCIL(rhs_t1);
    STENCIL(rhs_t2);
    STENCIL(v_t2);

    v_t2 = v_t2(0,0) + Dt*(AM3_2*rhs_t2(0,0) + AM3_1*rhs_t1(0,0) + AM3_0*rhs_t0(0,0));
}



extern "C"  __global__
void AdamsMoultonCorr3rd2(const double Dt, const double *_u_t2, const double *_v_t2,
                          const double *_rhsu_t0, const double *_rhsu_t1, const double *_rhsu_t2,
                          const double *_rhsv_t0, const double *_rhsv_t1, const double *_rhsv_t2,
                          const double *_h, const double *_zeta_t1, const double *_zeta_t2)
{
    const unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
    // Adams-Moulton 3rd order coefficients
    constexpr double AM3_2 =  5.0 / 12.0;
    constexpr double AM3_1 =  8.0 / 12.0;
    constexpr double AM3_0 = -1.0 / 12.0;

    if (i >= sz2D) return;

    STENCIL(rhsu_t0);
    STENCIL(rhsu_t1);
    STENCIL(rhsu_t2);
    STENCIL(rhsv_t0);
    STENCIL(rhsv_t1);
    STENCIL(rhsv_t2);
    STENCIL(u_t2);
    STENCIL(v_t2);
    STENCIL(zeta_t1);
    STENCIL(zeta_t2);
    STENCIL(h);

    auto D_t1 = zeta_t1 + h;
    auto D_t2 = zeta_t2 + h;
    auto D_t1U = RtoU(D_t1);
    auto D_t1V = RtoV(D_t1);
    auto D_t2U = RtoU(D_t2);
    auto D_t2V = RtoV(D_t2);

    u_t2 = (u_t2(0,0)*D_t1U + Dt*(AM3_2*rhsu_t2(0,0) + AM3_1*rhsu_t1(0,0) + AM3_0*rhsu_t0(0,0)))/D_t2U;
    v_t2 = (v_t2(0,0)*D_t1V + Dt*(AM3_2*rhsv_t2(0,0) + AM3_1*rhsv_t1(0,0) + AM3_0*rhsv_t0(0,0)))/D_t2V;
}