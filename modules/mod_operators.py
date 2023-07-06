import os.path

import mod_grid
import cupy as cp
from misc import *


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
om_v = None
on_u = None
sz = None


def initModule(GRID):
    # Initialize some global variables that are used extensively in the module.

    global G, shp, on_u, om_v, sz
    G = GRID
    on_u = G.on_u
    om_v = G.om_v
    shp = G.on_u.shape   # I chosen on_u, but it could've been any other array.
    sz  = G.on_u.size

    initializeCPPKernels((1,), (1,), (1, shp[0], shp[1]))


preamble2D = r'''
            #define STENCIL(var) Stencil<double>  var((double *)(&(_##var[i])))
            ptrdiff_t strideJ, strideI;

            void initializeStrides(const ptrdiff_t *strides)
            {
                //strideI = 1;
                strideJ = strides[0]/strides[1];
            }

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
                T * const p;    
            public:
                Stencil(T * _p): p(_p)
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

# -------------------------------------------------------------








loaded_from_source = r'''#define STENCIL(var) Stencil<double>  var((double *)(&(_##var[i])))
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

double RtoU(const double *_R, const unsigned int i)
{
    
    // Initializes the stencil variables that allow to access to relative elements.
    STENCIL(R);
    
    if ((i % strideJ)==0) return 0.0;       
                   
    return (R(0, 0) + R(-1, 0)) * 0.5;

}

double RtoV(const double *_R, const unsigned int i)
{
       
    STENCIL(R);
    
    if ((i/strideJ)==0) return 0.0;
            
    return (R(0, 0) + R(0, -1)) * 0.5;
}

double DERtoU(const double *_R, const double *_on_u, const unsigned int i)
{
    STENCIL(R);
    STENCIL(on_u);
    
    if ((i % strideJ)==0) return 0.0; 

    return on_u(0,0)*(R(0, 0) - R(0, -1));
}

double DNRtoV(const double *_R, const double *_om_v, const unsigned int i)
{
    STENCIL(R);
    STENCIL(om_v);
    
    if ((i / strideJ)==0) return 0.0; 

    return om_v(0,0)*(R(0, 0) - R(-1, 0));    
}

__global__ void computeMomentumRHS(const double *_h, const double *_gzeta, double *_gzeta2, const double *_on_u, const double *_om_v,
                                   const double *_rhs_ubar, const double *_rhs_vbar)
{
    unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
    strideJ = 502;
    
    if ((i % strideJ)==0 || (i/strideJ)==0) return;
    
    //STENCIL(h);
    //STENCIL(gzeta);
    STENCIL(rhs_ubar);
    STENCIL(rhs_vbar);
    
    
    double g = 9.8; 
    // printf("%li %li %li %li %p***\n", i, blockDim.x, blockIdx.x, threadIdx.x, _h);
    // rhs_ubar = RtoU(_h,i); 
    rhs_ubar = 0.5*g*(RtoU(_h,i)*DERtoU(_gzeta,_on_u,i) + DERtoU(_gzeta2,_on_u,i));
    rhs_vbar = 0.5*g*(RtoV(_h,i)*DNRtoV(_gzeta,_om_v,i) + DNRtoV(_gzeta2,_om_v,i));
}   

}
'''

# module = cp.RawModule(code=loaded_from_source, options=('-default-device',))
filename = os.path.join(exePath, r'modules/mod_operators.c')
with open(filename, 'r') as file:
    code = file.read()
module = cp.RawModule(code=code, options=('-default-device',))
computeMomentumRHS = module.get_function('computeMomentumRHS')
initOperators = module.get_function('initOperators')



filename = os.path.join(exePath, r'modules/mod_cppkernels.cpp')
with open(filename, 'r') as file:
    code = file.read()
moduleCPPKernels = cp.RawModule(code=code, options=('-default-device', '--std=c++17'))
initializeCPPKernels = moduleCPPKernels.get_function('initialize')

computeMomentumRHS3 = moduleCPPKernels.get_function('computeMomentumRHS')
computeZetaRHS3     = moduleCPPKernels.get_function('computeZetaRHS')
aaa     = moduleCPPKernels.get_function('aaa')
AdamsMoultonCorr3rd = moduleCPPKernels.get_function('AdamsMoultonCorr3rd')
Pred = moduleCPPKernels.get_function('Pred')



# initOperators = moduleCPPKernels.get_function('initOperators')


# -------------------------------------------------------------
copyBC = cp.ElementwiseKernel(
    '''raw float64 _R, raw int32 i1, raw int32 i2, raw float64 varForStrides''',
    '',
    preamble='',
    operation=r'''
        // Initializes the stencil variables that allow to access to relative elements.
        // initializeStrides(varForStrides.strides());

        // STENCIL(R);

        //printf("44444 - %i\n",  i);
        int ii1 = i1[i];
        int ii2 = i2[i];
        
        //printf("44444 - %i  %i  %i\n", ii1, ii2, i);

        *((double *)&(_R[ii1])) = _R[ii2];

       ''',
    name='copyBC',
    options=('-default-device',))

setBC = cp.ElementwiseKernel(
    '''raw float64 _R, raw int32 i1''',
    '',
    preamble='',
    operation=r'''
        int ii1 = i1[i];

        *((double *)&(_R[ii1])) = 0.0;

       ''',
    name='setBC',
    options=('-default-device',))


# -------------------------------------------------------------



RtoU_CUDA = cp.ElementwiseKernel(
    '''raw float64 _R, raw float64 varForStrides''',
    'raw float64 _U',
    preamble=preamble2D,
    operation=r'''
        // Initializes the stencil variables that allow to access to relative elements.
        initializeStrides(varForStrides.strides());
        
        STENCIL(R);
        STENCIL(U);
        
        if ((i % strideJ)==0) { U = 0.0; return; };        
                       
        U = (R(0, 0) + R(0, -1)) * 0.5;

       ''',
    name='RtoU_CUDA',
    options=('-default-device',))






def RtoU(R):

    return RtoU_CUDA.__call__(R, on_u, size=sz)



RtoV_CUDA = cp.ElementwiseKernel(
    '''raw float64 _R, raw float64 varForStrides''',
    'raw float64 _V',
    preamble=preamble2D,
    operation=r'''
        // Initializes the stencil variables that allow to access to relative elements.
        initializeStrides(varForStrides.strides());
        
        STENCIL(R);
        STENCIL(V);
        
        if ((i/strideJ)==0) { V = 0.0; return; };
                
        V = (R(0, 0) + R(-1, 0)) * 0.5;

       ''',
    name='RtoV_CUDA',
    options=('-default-device',))



def RtoV(R):

    # res = cp.zeros(shp[0]*shp[1])
    # res[:,1:] = RtoV_CUDA(R.reshape(shp)[:,1:], G.on_u[:,1:], size=G.on_u[:,1:].size).reshape(G.on_u[:,1:].shape)
    return RtoV_CUDA(R, on_u, size=sz)



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

        if ((i % strideJ)==0 || (i/strideJ)==0) { R = 0.0; return; };
        
        R = ( (U(0, 1)*on_u(0, 1) - U(0, 0)*on_u(0, 0)) + (V(1, 0)*om_v(1, 0) - V(0, 0)*om_v(0, 0)) )*pm(0, 0)*pn(0, 0);
        

       ''',
    name='RtoV_CUDA',
    options=('-default-device',))

import numba.cuda
# @numba.cuda.jit
# @numba.cuda.jit('float64[:](float64[:], float64[:])', device=True, inline=True)
def divUVtoR(U, V):
    res: cp.ndarray = cp.zeros(U.shape, dtype = cp.float64)
    res[:] = divUVtoR_CUDA(U, V, G.pm, G.pn, on_u, om_v, on_u, size=sz)
    return res





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
            
            if ((i % strideJ)==0) { U = 0.0; return; };

            U = on_u(0,0)*(R(0, 0) - R(0, -1));

           ''',
    name = 'DXRtoU_CUDA',
    options = ('-default-device',))



def DξRtoU(R):
    return DξRtoU_CUDA(R, on_u, on_u, size=sz)




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
            
            if ((i/strideJ)==0) { V = 0.0; return; };

            V = om_v(0,0)*(R(0, 0) - R(-1, 0));


           ''',
    name = 'DERtoV_CUDA',
    options = ('-default-device',))


def DηRtoV(R):
    return DηRtoV_CUDA(R, om_v, om_v, size=sz)