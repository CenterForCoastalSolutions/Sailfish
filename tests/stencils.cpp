#define STENCIL(var) Stencil<double>  var((double *)(&(_##var[i])))


ptrdiff_t strideI, strideJ, strideK;


void initializeStrides(const ptrdiff_t *strides)
{
    strideI = 1;
    strideJ = strides[1]/strides[2];
    strideK = strides[0]/strides[2];
}

constexpr int iwest  = 1;
constexpr int isouth = 2;
constexpr int ieast  = 3;
constexpr int inorth = 4;



template<typename T>
class Stencil
{
    T *p;

public:
    Stencil(T *_p): p(_p)
    {
    }

    size_t relIdx(size_t const k, size_t const j, size_t const i) const
    {
        return i*strideI + j*strideJ + k*strideK;
    }


    size_t relIdx(size_t const BCType) const
    {
        const int j = ((BCType-1) & 2) - 1;
        const int i = ((BCType-1) & 1)*j;
        return i*strideI + j*strideJ;
    }



//    T &operator()(size_t const k, size_t const j, size_t const i) const
//    {
//        return *(p + i*strideI + j*strideJ + k*strideK);
//    }

    Stencil &operator()(size_t const k, size_t const j, size_t const i) const
    {
        return Stencil(p + i*strideI + j*strideJ + k*strideK);
    }

    T &operator T(Stencil) {return *p}


    T operator=(T val) const
    {
        return (*p = val);
    }
};