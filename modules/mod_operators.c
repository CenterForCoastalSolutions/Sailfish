#define STENCIL(var) Stencil<double>  var((double *)(_##var + i))
ptrdiff_t strideJ, strideI;

template<typename T>
class Stencil
{
    T * const p;
public:
    Stencil(T * _p): p(_p)
    {
    }
    T &operator()(int const j, int const i) const
    {
        return *(p + j*strideJ + i);
    }
    T operator=(T val) const
    {
        return (*p = val);
    }
};



extern "C" {

    // Global variables.
    unsigned int sz2D, sz3D, szI, szJ, szK;

    double RtoU(const Stencil<double> &R, const unsigned int i)
    {
        return (R(0, 0) + R(-1, 0)) * 0.5;
    }

    double RtoV(const Stencil<double> &R, const unsigned int i)
    {
        return (R(0, 0) + R(0, -1)) * 0.5;
    }

    double DERtoU(const Stencil<double> &R, const Stencil<double> &on_u, const unsigned int i)
    {
        return on_u(0,0)*(R(0, 0) - R(0, -1));
    }

    double DNRtoV(const Stencil<double> &R, const Stencil<double> &om_v, const unsigned int i)
    {
        return om_v(0,0)*(R(0, 0) - R(-1, 0));
    }




    double _RtoU(const double *_R, const unsigned int i)
    {
        STENCIL(R);

        return (R(0, 0) + R(-1, 0)) * 0.5;
    }

    double _RtoV(const double *_R, const unsigned int i)
    {
        STENCIL(R);

        return (R(0, 0) + R(0, -1)) * 0.5;
    }

    double _DERtoU(const double *_R, const double *_on_u, const unsigned int i)
    {
        STENCIL(R);
        STENCIL(on_u);

        return on_u(0,0)*(R(0, 0) - R(0, -1));
    }

    double _DNRtoV(const double *_R, const double *_om_v, const unsigned int i)
    {
        STENCIL(R);
        STENCIL(om_v);

        return om_v(0,0)*(R(0, 0) - R(-1, 0));
    }


//----------------------------------------------------------------------------------

    __global__ void initOperators(unsigned int sizeK, unsigned int sizeJ, unsigned int sizeI)
    {
        // Stores some global variables that are used a lot (instead of passing the as parameters on each function).
        szI  = sizeI;
        szJ  = sizeJ;
        szK  = sizeK;
        sz2D = szI*szJ;
        sz3D = szI*szJ*szK;
    }


    __global__ void computeMomentumRHS(const double *_h, const double *_gzeta, double *_gzeta2, const double *_on_u, const double *_om_v,
                                       const double *_rhs_ubar, const double *_rhs_vbar, const double g)
    {
        const unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;

        if (((i % szJ) == 0 || (i/szJ) == 0) && i >= sz2D) return;

        STENCIL(h);
        STENCIL(gzeta);
        STENCIL(gzeta2);
        STENCIL(rhs_ubar);
        STENCIL(rhs_vbar);
        STENCIL(on_u);
        STENCIL(om_v);

        rhs_ubar = 0.5*g*(RtoU(h,i)*DERtoU(gzeta,on_u,i) + DERtoU(gzeta2,on_u,i));
        rhs_vbar = 0.5*g*(RtoV(h,i)*DNRtoV(gzeta,om_v,i) + DNRtoV(gzeta2,om_v,i));
    }


//    __global__ void computeMomentumRHS(const double *_D, const double *_gzeta, double *_gzeta2, const double *_on_u, const double *_om_v,
//                                       const double *_rhs_ubar, const double *_rhs_vbar)
//    {
//        const unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
//        strideJ = 502;
//
//
//        # compute the water column depth
//        D = zeta + h
//
//        DU = ubar*RtoU(D)
//        DV = vbar*RtoV(D)
//
//        return divUVtoR(DU, DV)
//    }

}