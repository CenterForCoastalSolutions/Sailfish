extern "C" {


#define STENCIL(var) Stencil<double>  var((double *)(&(_##var[i])))
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

__global__ double RtoU(const double *_R, const unsigned int i)
{
    strideJ = _strideJ;

    // Initializes the stencil variables that allow to access to relative elements.
    STENCIL(R);

    if ((i % strideJ)==0) return 0.0;

    return (R(0, 0) + R(-1, 0)) * 0.5;

}

__global__ double RtoV(const double *_R, const unsigned int i)
{

    STENCIL(R);
    STENCIL(V);

    if ((i/strideJ)==0) return 0.0;

    return (R(0, 0) + R(0, -1)) * 0.5;
}

__global__ double DERtoU(const double *_R, const double *_on_u, const unsigned int i)
{
    STENCIL(R);
    STENCIL(on_u);

    if ((i % strideJ)==0) return 0.0;

    return on_u(0,0)*(R(0, 0) - R(0, -1));
}

__global__ double DξRtoV(const double *_R, const double *_om_, const unsigned int i)v
{
    STENCIL(R);
    STENCIL(om_v);

    if ((i % strideJ)==0) return 0.0;

    return om_v(0,0)*(R(0, 0) - R(-1, 0));
}

__global__ void computeMomentumRHS(const double *_h, const double *_gzeta, const double *_rhs_ubar, const double *_rhs_vbar)
{
    unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;

    STENCIL(h);
    STENCIL(gzeta);
    STENCIL(rhs_ubar);
    STENCIL(rhs_vbar);


//    rhs_ubar = 0.5*g*(RtoU(h)*DERtoU(gzeta) + DERtoU(gzeta*gzeta))
//    rhs_vbar = 0.5*g*(RtoV(h)*DNRtoV(gzeta) + DNRtoV(gzeta*gzeta))
}

}