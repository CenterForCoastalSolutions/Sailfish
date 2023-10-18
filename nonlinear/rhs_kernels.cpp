// ηξΔ

#include "./mod_cppkernels.h"


extern "C"  __global__
void horizontalAdvection(const double *_u,  const double *_v, const double *_Huon, const double *_Hvom,
                         const double *_ru, const double *_rv, const int N)
{
    const unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;

    STENCIL3D(u);
    STENCIL3D(v);
    STENCIL3D(Huon);
    STENCIL3D(Hvom);
    STENCIL3D(ru);
    STENCIL3D(rv);

    if (((i % szI) == 0 || ((i % szI) == (szI - 1) || (i/szI) == 0) || (i/szI) == (szJ - 1)) || i >= sz2D) return;


//    ru = ru[nrhs,:,:,:]
//    rv = rv[nrhs,:,:,:]
    double Gadv = 0.0;

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
        setK(k);  // Sets the current k index.


        // Add horizontal advection of momentum to the RHS vector.
        ru -= DξRtoU(UFξR) + DηPtoU(UFηP);
        rv -= DξPtoV(VFξP) + DηRtoV(VFηR);
    }

}

//--------------------------------------------------------------------------

extern "C"  __global__
void HorizontalHomogeneousBC()
{


}


//--------------------------------------------------------------------------

extern "C"  __global__
void verticalAdvection(double const *_u, double const *_W, double const *_ru)
{
    STENCIL3D(u);
    STENCIL3D(W);
    STENCIL3D(ru);

    // Product of the fourth order centered interpolations (at R points) of u and W  (internal vertical nodes)
    for (int k=1; k<N-2; k++)
    {
        auto Fsigma = UtoUW_4th(u)*WtoUW_4th(W);
    }

//    ru[nrhs,:,:,:] -= Dsigma(Fsigma)
}


void UtoUW_4th(const double *_var)
{

//    for (int k=1; k<N-2; k++)
//    {
//        setK(k);
//        res = (9.0/16.0)*(var(k,0,0) + var(k+1,0,0)) - (1.0/16.0)*(var(k-1,0,0) + var(k+2,0,0));
//    }
//
//    // vertical boundary conditions
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
