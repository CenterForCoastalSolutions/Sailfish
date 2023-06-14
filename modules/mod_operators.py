import mod_grid
import cupy as cp


# import rmm
# pool = rmm.mr.PoolMemoryResource(
#     rmm.mr.ManagedMemoryResource(),
#     initial_pool_size=5*(2**30),
#     maximum_pool_size=5*(2**30)
# )
# rmm.mr.set_current_device_resource(pool)
# cp.cuda.set_allocator(rmm.rmm_cupy_allocator)
#

G = None   # GRID
shp = None


def initModule(GRID):
    # Initialize some global variables that are used extensively in the module.

    global G, shp
    G = GRID
    shp = G.on_u.shape   # I chosen on_u, but it could've been any other array.



preamble2D = r'''
            #define STENCIL(var) Stencil<double>  var((double *)(&(_##var[i])))
            ptrdiff_t strideJ, strideI;

            void initializeStrides(const ptrdiff_t *strides)
            {
                strideI = 1;
                strideJ = strides[0]/strides[1];
            }

            template<typename T>
            class Stencil
            {
                T *p;    
            public:
                Stencil(T *_p): p(_p)
                {
                }
                T &operator()(size_t const j, size_t const i) const
                {
                    return *(p + j*strideJ + i*strideI);
                }
                T operator=(T val) const
                {
                    return (*p = val);
                }
            };
            '''

preamble3D = r'''
            #define STENCIL(var) Stencil<double>  var((double *)(&(_##var[i])))
            ptrdiff_t strideI, strideJ, strideK;

            void initializeStrides(const ptrdiff_t *strides)
            {
                strideI = 1;
                strideJ = strides[1]/strides[2];
                strideK = strides[0]/strides[2];
            }

            template<typename T>
            class Stencil
            {
                T *p;    
            public:
                Stencil(T *_p): p(_p)
                {
                }
                T &operator()(size_t const k, size_t const j, size_t const i) const
                {
                    return *(p + k*strideK + j*strideJ + i*strideI);
                }
                T operator=(T val) const
                {
                    return (*p = val);
                }
            };
            '''

RtoU_CUDA = cp.ElementwiseKernel(
    '''raw float64 _R, raw float64 varForStrides''',
    'raw float64 _U',
    preamble=preamble2D,
    operation=r'''
        // Initializes the stencil variables that allow to access to relative elements.
        initializeStrides(varForStrides.strides());
        STENCIL(R);
        STENCIL(U);

        U = (R(0, 0) + R(-1, 0)) * 0.5;

       ''',
    name='RtoU_CUDA',
    options=('-default-device',))







def RtoU(R):
    mempool = cp.get_default_memory_pool()
    pinned_mempool = cp.get_default_pinned_memory_pool()

    # Create an array on CPU.
    # NumPy allocates 400 bytes in CPU (not managed by CuPy memory pool).
    a_cpu = cp.ndarray(100, dtype=cp.float32)
    print('nbytes', a_cpu.nbytes)  # 400

    # You can access statistics of these memory pools.
    print(mempool.used_bytes())  # 0
    print(mempool.total_bytes())  # 0
    print(pinned_mempool.n_free_blocks())  # 0
    print('max:', mempool.get_limit())
    # print('max2:', pinned_mempool.get_limit())
    res = cp.zeros(shp)

    a = R.reshape(shp)[:, 1:]
    print(mempool.used_bytes())
    b = G.on_u[:, 1:]
    print(mempool.used_bytes())
    c = RtoU_CUDA(a, b, size=G.on_u[:, 1:].size)
    print(mempool.used_bytes())

    res[:, 1:] = c.reshape(G.on_u[:, 1:].shape)
    print(mempool.used_bytes())
    # res[:,1:] = RtoU_CUDA(R.reshape(shp)[:,1:], G.on_u[:,1:], size=G.on_u[:,1:].size).reshape(G.on_u[:,1:].shape)
    return res.ravel()



RtoV_CUDA = cp.ElementwiseKernel(
    '''raw float64 _R, raw float64 varForStrides''',
    'raw float64 _V',
    preamble=preamble2D,
    operation=r'''
        // Initializes the stencil variables that allow to access to relative elements.
        initializeStrides(varForStrides.strides());
        STENCIL(R);
        STENCIL(V);

        V = (R(0, 0) + R(0, -1)) * 0.5;

       ''',
    name='RtoV_CUDA',
    options=('-default-device',))



def RtoV(R):

    res = cp.zeros(shp)
    res[1:,:] = RtoV_CUDA(R.reshape(shp)[1:,:], G.on_u[1:,:], size=G.on_u[1:,:].size).reshape(G.on_u[1:,:].shape)
    return res.ravel()



divUVtoR_CUDA = cp.ElementwiseKernel(
    '''raw float64 _U, raw float64 _V, raw float64 _pm, raw float64 _pn, raw float64 _on_u, raw float64 _om_v, raw float64 varForStrides''',
    'raw float64 _R',
    preamble=preamble2D,
    operation=r'''
        // Initializes the stencil variables that allow to access to relative elements.
        initializeStrides(varForStrides.strides());
        STENCIL(U);
        STENCIL(V);
        STENCIL(R);
        STENCIL(pm);
        STENCIL(pn);
        STENCIL(on_u);
        STENCIL(om_v);
        
        R = ( (U(0, 1)*on_u(0, 1) - U(0, 0)*on_u(0, 0)) + (V(1, 0)*om_v(1, 0) - V(0, 0)*om_v(0, 0)) )*pm(0, 0)*pn(0, 0);
        

       ''',
    name='RtoV_CUDA',
    options=('-default-device',))



def divUVtoR(U, V):
    res = cp.zeros(shp)
    res[1:,1:] = divUVtoR_CUDA(U.reshape(shp)[1:,1:], V.reshape(shp)[1:,1:], G.pm[:1,:1], G.pn[1:,1:], G.on_u[1:,1:], G.om_v[1:,1:],  G.on_u, size=G.om_v[1:,1:].size).reshape(G.om_v[1:,1:].shape)
    return res.ravel()





DξRtoU_CUDA = cp.ElementwiseKernel(
    '''raw float64 _R, raw float64 _on_u, raw float64 varForStrides''',
    'raw float64 _U',
    preamble = preamble2D,
    operation = r'''
            // Initializes the stencil variables that allow to access to relative elements.
            initializeStrides(varForStrides.strides());
            STENCIL(U);
            STENCIL(R);
            STENCIL(on_u);

            U = on_u(0,0)*(R(0, 0) - R(0, -1));


           ''',
    name = 'DXRtoU_CUDA',
    options = ('-default-device',))


def DξRtoU(R):
    res = cp.zeros(shp)
    res[:,1:] =  DξRtoU_CUDA(R.reshape(shp)[:,1:], G.on_u[:,1:], G.on_u, size=G.on_u[:,1:].size).reshape(G.on_u[:,1:].shape)
    return res.ravel()





DηRtoV_CUDA = cp.ElementwiseKernel(
    '''raw float64 _R, raw float64 _om_v, raw float64 varForStrides''',
    'raw float64 _V',
    preamble = preamble2D,
    operation = r'''
            // Initializes the stencil variables that allow to access to relative elements.
            initializeStrides(varForStrides.strides());
            STENCIL(V);
            STENCIL(R);
            STENCIL(om_v);

            V = om_v(0,0)*(R(0, 0) - R(-1, 0));


           ''',
    name = 'DERtoV_CUDA',
    options = ('-default-device',))


def DηRtoV(R):
    res = cp.zeros(shp)
    res[1:,:] = DηRtoV_CUDA(R.reshape(shp)[1:,:], G.om_v[1:,:], G.om_v, size=G.om_v[1:,:].size).reshape(G.om_v[1:,:].shape)
    return res.ravel()