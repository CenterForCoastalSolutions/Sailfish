// ηξΔ

#include "./mod_cppkernels.h"


extern "C"  __global__
void addCoriolis(const double *_fomn, const double *_u, const double *_v, const double *_ru, const double *_rv, const double *_Hz)
{
    const unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
    const int N = szK;

    if (!isRNode(i)) return;

    int K = 0;
    STENCIL(fomn);
    STENCIL3D(u,  K);
    STENCIL3D(v,  K);
    STENCIL3D(ru, K);
    STENCIL3D(rv, K);
    STENCIL3D(Hz, K);


    for (K = 0; K  <N; K++)
    {
        printf("> %i, %i \n", K, i);
        auto cff = (Hz*fomn);

        auto UF = cff*(VtoR(v));
        auto VF = cff*(UtoR(u));

        ru += RtoU(UF);
        rv -= RtoV(VF);
    }

}


extern "C"  __global__
void addCurvedGridTerms(void)
{
// TODO: implement
//    DO j=JstrV-1,Jend
//        DO i=IstrU-1,Iend
//             cff1=0.5*(v(i,j  ,k,nrhs) + v(i,j+1,k,nrhs))
//             cff2=0.5*(u(i  ,j,k,nrhs) + u(i+1,j,k,nrhs))
//             cff3=cff1*dndx(i,j)
//             cff4=cff2*dmde(i,j)
//
//             cff=Hz(i,j,k)*(cff3-cff4)
//             UFx(i,j)=cff*cff1
//             VFe(i,j)=cff*cff2
//        END DO
//    END DO
//
//    DO j=Jstr,Jend
//        DO i=IstrU,Iend
//            cff1=0.5*(UFx(i,j)+UFx(i-1,j))
//            ru(i,j,k,nrhs) = ru(i,j,k,nrhs) + cff1
//        END DO
//    END DO
//
//    DO j=JstrV,Jend
//        DO i=Istr,Iend
//             cff1=0.5*(VFe(i,j)+VFe(i,j-1))
//             rv(i,j,k,nrhs) = rv(i,j,k,nrhs) - cff1
//        END DO
//    END DO
}

extern "C"  __global__
void horizontalAdvection(const double *_u,  const double *_v, const double *_Huon, const double *_Hvom,
                         const double *_ru, const double *_rv, const int *BC, const double Gadv, const int N)
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


    for (K=2; K<N-2; K++)
    {
        // Product of the fourth order centered interpolations (at R points) of u (or v) and W  (internal vertical nodes,
        // the ones near the surface or botton are dealt with as BC)
        auto FsigmaUW = UtoUW_4th(u)*WtoUW_4th(W);
        auto FsigmaVW = VtoVW_4th(v)*WtoVW_4th(W);

        ru -= DσUWtoU(FsigmaUW);
        rv -= DσVWtoV(FsigmaVW);
    }

    // vertical boundary conditions (TODO: can it be done inside DσVWtoV?, also revise, clearly there is an error with K and the indices used for k)
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