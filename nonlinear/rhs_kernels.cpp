// ηξΔ

#include "./mod_cppkernels.h"




extern "C"  __global__
void horizontalAdvection(const double *_u,  const double *_v, const double *_Huon, const double *_Hvom,
                         const double *_ru, const double *_rv, const int *BC, cont double Gadv, const int N)
{
    const unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;

    int K = 0;  // All stencils share this variable.
    STENCIL3D(u,    K);
    STENCIL3D(v,    K);
    STENCIL3D(Huon, K);
    STENCIL3D(Hvom, K);
    STENCIL3D(ru,   K);
    STENCIL3D(rv,   K);


    if (i >= sz2D || BC[i] >= 0)
    {
        ru = 0.0;
        rv = 0.0;
        return;
    };

//    if (((i % szI) == 0 || ((i % szI) == (szI - 1) || (i/szI) == 0) || (i/szI) == (szJ - 1)) || i >= sz2D) return;

    auto uR = UtoR(u);
    auto uP = UtoP(u);
    auto vR = VtoR(v);
    auto vP = VtoP(v);

    auto HuonR = UtoR(Huon);
    auto HuonP = UtoP(Huon);
    auto HvomP = VtoP(Hvom);
    auto HvomR = VtoR(Hvom);

    auto uξξR = upwindUtoR(DξξUtoU(u), uR);
    auto uηηP = upwindUtoP(DηηUtoU(u), HvomP);
    auto vξξP = upwindVtoP(DξξVtoV(v), HuonP);
    auto vηηR = upwindVtoR(DηηVtoV(v), vR);

    auto HuξξR = UtoR(DξξUtoU(Huon));
    auto HuηηP = UtoP(DηηUtoU(Huon));
    auto HvξξP = VtoP(DξξVtoV(Hvom));
    auto HvηηR = VtoR(DηηVtoV(Hvom));

    // See equations in section "Advection scheme" of the accompanying document.
    auto UFξR = (uR + Gadv*uξξR*(HuonR + Gadv*HuξξR));
    auto UFηP = (uP + Gadv*uηηP*(HvomP + Gadv*HvξξP));

    auto VFξP = (vP + Gadv*uηηP*(HuonP + Gadv*HuηηP));  // TODO: CHECK all ξξ and ηη
    auto VFηR = (vR + Gadv*uξξR*(HvomR + Gadv*HvηηR));



    for (int k=0; k<N; k++)
    {
        K = k;  // Sets the current k index.

        // Add horizontal advection of momentum to the RHS vector.
        ru -= DξRtoU(UFξR) + DηPtoU(UFηP);
        rv -= DξPtoV(VFξP) + DηRtoV(VFηR);
    }

}

//--------------------------------------------------------------------------

extern "C"  __global__
void verticalAdvection(double const *_u, double const *_v, double const *_W, double const *_ru, double const *_rv, int const *BC, const int N)
{
    const unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;

//    if (((i % szI) == 0 || ((i % szI) >= (szI - 2) || (i/szI) == 0) || (i/szI) >= (szJ - 2)) || i >= sz2D) return;

    if (i >= sz2D || BC[i] >= 0)
    {
        return;
    }

    int K = 0;  // All stencils share this variable.
    STENCIL3D(u,  K);
    STENCIL3D(v,  K);
    STENCIL3D(W,  K);
    STENCIL3D(ru, K);
    STENCIL3D(rv, K);


    for (int k=2; k<N-2; k++)
    {
        // Product of the fourth order centered interpolations (at R points) of u (or v) and W  (internal vertical nodes,
        // the ones near the surface or botton are dealt with as BC)
        auto FsigmaUW = UtoUW_4th(u)*WtoUW_4th(W);
        auto FsigmaVW = VtoVW_4th(v)*WtoVW_4th(W);

        K = k;  // Sets current k index.

        ru -= DσUWtoU(FsigmaUW);
        rv -= DσVWtoV(FsigmaVW);
    }

    // vertical boundary conditions (TODO: can it be done inside DσVWtoV?)
    K = 1;
    ru -= (9.0/16.0)*(ru(  0,0,0) + ru(  1,0,0)) - (1.0/16.0)*(ru(  0,0,0) + ru(  2,0,0));
    rv -= (9.0/16.0)*(rv(  0,0,0) + rv(  1,0,0)) - (1.0/16.0)*(rv(  0,0,0) + rv(  2,0,0));

    K = N-2;
    ru -= (9.0/16.0)*(ru(N-1,0,0) + ru(N-2,0,0)) - (1.0/16.0)*(ru(N-1,0,0) + ru(N-3,0,0));
    rv -= (9.0/16.0)*(rv(N-1,0,0) + rv(N-2,0,0)) - (1.0/16.0)*(rv(N-1,0,0) + rv(N-3,0,0));

    K = 0;
    ru = 0.0;
    rv = 0.0;

    K = N-1;
    ru = 0.0;
    rv = 0.0;

}



extern "C"  __global__
void verticalHomogeneousBC()
{
//
//    setK(1);
//    res = (9.0/16.0)*(var(  0,0,0) + var(  1,0,0)) - (1.0/16.0)*(var(  0,0,0) + var(2,0,0));
//    setK(N-2);
//    res = (9.0/16.0)*(var(N-1,0,0) + var(N-2,0,0)) - (1.0/16.0)*(var(N-1,0,0) + var(N-3,0,0));
//
//    setK(0);
//    res = 0.0;
//    setK(N-1);
//    res = 0.0;

}