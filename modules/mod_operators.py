import os.path
import codecs
import unicodedata

import mod_grid
import cupy as cp
from misc import *


import rmm
pool = rmm.mr.PoolMemoryResource(
    rmm.mr.ManagedMemoryResource(),
    initial_pool_size=5*(2**35),
    maximum_pool_size=5*(2**35)
)
rmm.mr.set_current_device_resource(pool)
cp.cuda.set_allocator(rmm.rmm_cupy_allocator)


G    = None   # GRID
shp  = None
om_v = None
on_u = None
sz   = None
bksz = None
grsz = None


def initModule(GRID):
    # Initialize some global variables that are used extensively in the module.

    global G, shp, on_u, om_v, sz, bksz, grsz
    G = GRID
    on_u = G.on_u
    om_v = G.om_v
    shp = G.on_u.shape   # I chosen on_u, but it could've been any other array.
    sz  = G.on_u.size
    bksz = (blockSize,)
    grsz = (sz//blockSize + 1,)

    initializeCPPKernels   ((1,), (1,), (GRID.N, shp[0], shp[1], G.on_u, G.om_v, G.pn, G.pm))
    initializeRHSKernels   ((1,), (1,), (GRID.N, shp[0], shp[1], G.on_u, G.om_v, G.pn, G.pm))
    initializeStep3DKernels((1,), (1,), (GRID.N, shp[0], shp[1], G.on_u, G.om_v, G.pn, G.pm))


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








# filename = os.path.join(exePath, r'modules/mod_operators.c')
# with open(filename, 'r') as file:
#     code = file.read()
# module = cp.RawModule(code=code, options=('-default-device',))
# computeMomentumRHS = module.get_function('computeMomentumRHS')
# initOperators = module.get_function('initOperators')



filePath = os.path.dirname(os.path.abspath(__file__))

filename = os.path.join(exePath, r'modules/mod_cppkernels.cpp')
with open(filename, 'r') as file:
    code = file.read()
with codecs.open(filename, encoding='utf-8') as file:
    code = file.read()
code = code.replace('η', 'E')
code = code.replace('ξ', 'X')
code = code.replace('Δ', 'D')
code = code.replace('σ', 'sig')
code = unicodedata.normalize('NFKD', code).encode('ascii', 'ignore').decode('ascii')
moduleCPPKernels = cp.RawModule(code=code, options=('-default-device', '--restrict', '--std=c++17', r'-I%s' % filePath))
initializeCPPKernels = moduleCPPKernels.get_function('initialize')

# computeMomentumRHS3 = moduleCPPKernels.get_function('computeMomentumRHS')
computeMomentumRHSCorr = moduleCPPKernels.get_function('computeMomentumRHSCorr')
computeMomentumRHSPred = moduleCPPKernels.get_function('computeMomentumRHSPred')
computeZetaRHS         = moduleCPPKernels.get_function('computeZetaRHS')
computeZetaPred        = moduleCPPKernels.get_function('computeZetaPred')
AdamsMoultonCorr3rd    = moduleCPPKernels.get_function('AdamsMoultonCorr3rd')
AdamsMoultonCorr3rd2   = moduleCPPKernels.get_function('AdamsMoultonCorr3rd')
computeMomentumPred    = moduleCPPKernels.get_function('computeMomentumPred')
computeMomentumRHSCorr = moduleCPPKernels.get_function('computeMomentumRHSCorr')

filename = os.path.join(exePath, r'nonlinear/rhs_kernels.cpp')
with codecs.open(filename, encoding='utf-8') as file:
    code = file.read()
code = code.replace('η', 'E')
code = code.replace('ξ', 'X')
code = code.replace('Δ', 'D')
code = code.replace('σ', 'sig')
code = unicodedata.normalize('NFKD', code).encode('ascii', 'ignore').decode('ascii')
moduleRHSKernels = cp.RawModule(code=code, options=('-default-device', '--restrict', '--std=c++17', r'-I%s' % filePath))
initializeRHSKernels = moduleRHSKernels.get_function('initialize')
horizontalAdvection  = moduleRHSKernels.get_function('horizontalAdvection')
verticalAdvection    = moduleRHSKernels.get_function('verticalAdvection')
addCoriolis          = moduleRHSKernels.get_function('addCoriolis')
vertIntegral         = moduleRHSKernels.get_function('vertIntegral')


filename = os.path.join(exePath, r'nonlinear/step3D_uv.kernels.cpp')
with codecs.open(filename, encoding='utf-8') as file:
    code = file.read()
code = code.replace('η', 'E')
code = code.replace('ξ', 'X')
code = code.replace('Δ', 'D')
code = code.replace('σ', 'sig')
code = unicodedata.normalize('NFKD', code).encode('ascii', 'ignore').decode('ascii')
moduleStep3DKernels = cp.RawModule(code=code, options=('-default-device', '--restrict', '--std=c++17', r'-I%s' % filePath))
initializeStep3DKernels  = moduleStep3DKernels.get_function('initialize')
# adjustBarotropicVelocity = moduleStep3DKernels.get_function('adjustBarotropicVelocity')
# correctBaroclinicVel     = moduleStep3DKernels.get_function('correctBaroclinicVel')
setVerticalVelEq         = moduleStep3DKernels.get_function('setVerticalVelEq')
set_maxflux              = moduleStep3DKernels.get_function('set_maxflux')
omega                    = moduleStep3DKernels.get_function('omega')
set_zeta                 = moduleStep3DKernels.get_function('set_zeta')
set_depth                = moduleStep3DKernels.get_function('set_depth')
step3d_UV                = moduleStep3DKernels.get_function('step3d_UV')
setLateralUVBCs          = moduleStep3DKernels.get_function('setLateralUVBCs')




# initOperators = moduleCPPKernels.get_function('initOperators')


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
    name = 'DNRtoU_CUDA',
    options = ('-default-device',))



def DξRtoU(R):
    return DξRtoU_CUDA(R, on_u, on_u, size=sz)




DξRtoU_CUDA = cp.ElementwiseKernel(
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


def DξRtoU(R):
    return DξRtoU_CUDA(R, om_v, om_v, size=sz)