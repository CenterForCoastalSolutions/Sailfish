import mod_grid
import numpy as cp


#   IJwaterR   Water points IJ couter for RHO-points masked variables.
#   IJwaterU   Water points IJ couter for   U-  points masked variables.
#   IJwaterV   Water points IJ couter for   V-points masked variables.
#   Hz         Thicknesses (m) of vertical RHO-points.
#   Huon       Total U-momentum flux term, Hz*u/pn.
#   Hvom       Total V-momentum flux term, Hz*v/pm.
#   IcePress   Pressure under the ice shelf at RHO-points.
#   Rscope     Adjoint sensitivity spatial scope mask at RHO-points.
#   Tcline     Width (m) of surface or bottom boundary layer where
#                higher vertical resolution is required during
#                stretching.
#   Uscope     Adjoint sensitivity spatial scope mask at U-points.
#   Vscope     Adjoint sensitivity spatial scope mask at V-points.
#   angler     Angle (radians) between XI-axis and true EAST at
#                RHO-points.
#   CosAngler  Cosine of curvilinear angle, angler.
#   SinAngler  Sine of curvilinear angle, angler.
#   dmde       ETA-derivative of inverse metric factor pm,
#                d(1/pm)/d(ETA).
#   dndx       XI-derivative  of inverse metric factor pn,
#                d(1/pn)/d(XI).
#   f          Coriolis parameter (1/s).
#   fomn       Compound term, f/(pm*pn) at RHO points.
#   grdscl     Grid scale used to adjust horizontal mixing according
#                to grid area.
#   h          Bottom depth (m) at RHO-points.
#   latp       Latitude (degrees_north) at PSI-points.
#   latr       Latitude (degrees_north) at RHO-points.
#   latu       Latitude (degrees_north) at U-points.
#   latv       Latitude (degrees_north) at V-points.
#   lonp       Longitude (degrees_east) at PSI-points.
#   lonr       Longitude (degrees_east) at RHO-points.
#   lonu       Longitude (degrees_east) at U-points.
#   lonv       Longitude (degrees_east) at V-points.
#   Mylon      Longitude work array for regridding.
#   omm        RHO-grid area (meters2).
#   om_p       PSI-grid spacing (meters) in the XI-direction.
#   om_r       RHO-grid spacing (meters) in the XI-direction.
#   om_u       U-grid spacing (meters) in the XI-direction.
#   om_v       V-grid spacing (meters) in the XI-direction.
#   on_p       PSI-grid spacing (meters) in the ETA-direction.
#   on_r       RHO-grid spacing (meters) in the ETA-direction.
#   on_u       U-grid spacing (meters) in the ETA-direction.
#   on_v       V-grid spacing (meters) in the ETA-direction.
#   pm         Coordinate transformation metric "m" (1/meters)
#                associated with the differential distances in XI.
#   pmon_p     Compound term, pm/pn at PSI-points.
#   pmon_r     Compound term, pm/pn at RHO-points.
#   pmon_u     Compound term, pm/pn at U-points.
#   pmon_v     Compound term, pm/pn at V-points.
#   pn         Coordinate transformation metric "n" (1/meters)
#                associated with the differential distances in ETA.
#   pnom_p     Compound term, pn/pm at PSI-points.
#   pnom_r     Compound term, pn/pm at RHO-points.
#   pnom_u     Compound term, pn/pm at U-points.
#   pnom_v     Compound term, pn/pm at V-points.
#   pmask      Slipperiness time-independent mask at PSI-points:
#                (0=Land, 1=Sea, 2=no-slip).
#   rmask      Time-independent mask at RHO-points (0=Land, 1=Sea).
#   umask      Time-independent mask at U-points (0=Land, 1=Sea).
#   vmask      Time-independent mask at V-points (0=Land, 1=Sea).
#
#   pmask_full Full mask at PSI-points (0=dry, 1=wet, 2=no-slip).
#   rmask_full Full mask at RHO-points (0=dry, 1=wet).
#   rmask_full Full mask at   U-points (0=dry, 1=wet).
#   rmask_full Full mask at   V-points (0=dry, 1=wet).

#   xp         XI-coordinates (m) at PSI-points.
#   xr         XI-coordinates (m) at RHO-points.
#   xu         XI-coordinates (m) at U-points.
#   xv         XI-coordinates (m) at V-points.
#   yp         ETA-coordinates (m) at PSI-points.
#   yr         ETA-coordinates (m) at RHO-points.
#   yu         ETA-coordinates (m) at U-points.
#   yv         ETA-coordinates (m) at V-points.
#   zice       Depth of ice shelf cavity (m, negative) at RHO-points.
#   z0_r       Time independent depths (m) at horizontal RHO-points and vertical RHO-points.
#   z0_w       Time independent depths (m) at horizontal RHO-points and  vertical W-points.
#   z_r        Actual depths (m) at horizontal RHO-points and vertical RHO-points.
#   z_w        Actual depths (m) at horizontal RHO-points and vertical W-points.





grid = None




preamble2D = r'''
            #define STENCIL(var) Stencil<double>  var((double *)(&(_##var[i])))
            ptrdiff_t strideJ, strideI;

            void initializeStrides(const ptrdiff_t *strides)
            {
                strideJ = 1;
                strideI = strides[0]/strides[1];
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
                strideK = 1;
                strideJ = strides[1]/strides[2];
                strideI = strides[0]/strides[2];
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

def RtoU(R, U):
    return RtoU_CUDA(R, U, size=R.size)



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

def RtoV(R, V):
    return RtoV_CUDA(R, V, size=R.size)



divUVtoR_CUDA = cp.ElementwiseKernel(
    '''raw float64 _U, raw float64 _V, raw float64 pm, raw float64 pn, raw float64 on_u, raw float64 om_v, raw float64 varForStrides''',
    'raw float64 _R',
    preamble=preamble2D,
    operation=r'''
        // Initializes the stencil variables that allow to access to relative elements.
        initializeStrides(varForStrides.strides());
        STENCIL(U);
        STENCIL(V);
        STENCIL(R);
        
        R = (((U(0, 1)*on_u(0, 1) - U(0, 0)*on_u(0, 0)) + ((V(0, 1)*om_v(0, 1) - V(0, 0)*om_v(0, 0))) * pm * pn;
        

       ''',
    name='RtoV_CUDA',
    options=('-default-device',))



def divUVtoR(U, V, R):
    return divUVtoR_CUDA(U, V, R, size=R.size)



rhs_ubar(i,j)=*on_u(i,j)* ((h(i-1,j)+ h(i ,j))* (gzeta(i-1,j)-   gzeta(i  ,j))+ (gzeta2(i - 1, j) -  gzeta2(i, j)))


((h(j, i-1)-



cff1 = g






R = cp.arange(0, 25, dtype=cp.float64)
R = R.reshape(5, 5)
V = cp.zeros_like(R)
U = cp.zeros_like(R)
V = RtoV(R, V)
U = RtoU(R, U)
print('this is R:')
print(R)

print('\nthis is U:')
print(U)


print('this is R:')
print(R)

print('\nthis is V:')
print(V)

#def divUVtoR(u, v):
