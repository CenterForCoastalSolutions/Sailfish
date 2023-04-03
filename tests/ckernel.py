

preamble =
'''
template<class T>
class Stencil
{
    T *p;
    size_t strideI, strideJ, strideK;

public:
    Stencil(T *_p): p(_p)
    {
      INITIALIZE strideI, strideJ, strideK
    }

    T &operator()(size_t const k, size_t const j, size_t const i) const
    {
        return *(p + i*strideI + j*strideJ + k*strideK);
    }

    T operator=(T val) const
    {
        return (*p = val);
    }
}
'''