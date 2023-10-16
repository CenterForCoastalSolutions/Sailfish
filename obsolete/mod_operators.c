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

    double DηRtoU(const Stencil<double> &R, const Stencil<double> &on_u, const unsigned int i)
    {
        return on_u(0,0)*(R(0, 0) - R(0, -1));
    }

    double DξRtoV(const Stencil<double> &R, const Stencil<double> &om_v, const unsigned int i)
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


