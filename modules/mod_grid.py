import cupy as cp
import mod_param
                                                                      !
# !  IJwaterR   Water points IJ couter for RHO-points masked variables.  !
# !  IJwaterU   Water points IJ couter for   U-points masked variables.  !
# !  IJwaterV   Water points IJ couter for   V-points masked variables.  !
# !  Hz         Thicknesses (m) of vertical RHO-points.                  !
# !  Huon       Total U-momentum flux term, Hz*u/pn.                     !
# !  Hvom       Total V-momentum flux term, Hz*v/pm.                     !
# !  IcePress   Pressure under the ice shelf at RHO-points.              !
# !  Rscope     Adjoint sensitivity spatial scope mask at RHO-points.    !
# !  Tcline     Width (m) of surface or bottom boundary layer where      !
# !               higher vertical resolution is required during          !
# !               stretching.                                            !
# !  Uscope     Adjoint sensitivity spatial scope mask at U-points.      !
# !  Vscope     Adjoint sensitivity spatial scope mask at V-points.      !
# !  angler     Angle (radians) between XI-axis and true EAST at         !
# !               RHO-points.                                            !
# !  CosAngler  Cosine of curvilinear angle, angler.                     !
# !  SinAngler  Sine of curvilinear angle, angler.                       !
# !  dmde       ETA-derivative of inverse metric factor pm,              !
# !               d(1/pm)/d(ETA).                                        !
# !  dndx       XI-derivative  of inverse metric factor pn,              !
# !               d(1/pn)/d(XI).                                         !
# !  f          Coriolis parameter (1/s).                                !
# !  fomn       Compound term, f/(pm*pn) at RHO points.                  !
# !  grdscl     Grid scale used to adjust horizontal mixing according    !
# !               to grid area.                                          !
# !  h          Bottom depth (m) at RHO-points.                          !
# !  latp       Latitude (degrees_north) at PSI-points.                  !
# !  latr       Latitude (degrees_north) at RHO-points.                  !
# !  latu       Latitude (degrees_north) at U-points.                    !
# !  latv       Latitude (degrees_north) at V-points.                    !
# !  lonp       Longitude (degrees_east) at PSI-points.                  !
# !  lonr       Longitude (degrees_east) at RHO-points.                  !
# !  lonu       Longitude (degrees_east) at U-points.                    !
# !  lonv       Longitude (degrees_east) at V-points.                    !
# !  Mylon      Longitude work array for regridding.                     !
# !  omm        RHO-grid area (meters2).                                 !
# !  om_p       PSI-grid spacing (meters) in the XI-direction.           !
# !  om_r       RHO-grid spacing (meters) in the XI-direction.           !
# !  om_u       U-grid spacing (meters) in the XI-direction.             !
# !  om_v       V-grid spacing (meters) in the XI-direction.             !
# !  on_p       PSI-grid spacing (meters) in the ETA-direction.          !
# !  on_r       RHO-grid spacing (meters) in the ETA-direction.          !
# !  on_u       U-grid spacing (meters) in the ETA-direction.            !
# !  on_v       V-grid spacing (meters) in the ETA-direction.            !
# !  pm         Coordinate transformation metric "m" (1/meters)          !
# !               associated with the differential distances in XI.      !
# !  pmon_p     Compound term, pm/pn at PSI-points.                      !
# !  pmon_r     Compound term, pm/pn at RHO-points.                      !
# !  pmon_u     Compound term, pm/pn at U-points.                        !
# !  pmon_v     Compound term, pm/pn at V-points.                        !
# !  pn         Coordinate transformation metric "n" (1/meters)          !
# !               associated with the differential distances in ETA.     !
# !  pnom_p     Compound term, pn/pm at PSI-points.                      !
# !  pnom_r     Compound term, pn/pm at RHO-points.                      !
# !  pnom_u     Compound term, pn/pm at U-points.                        !
# !  pnom_v     Compound term, pn/pm at V-points.                        !
# !  pmask      Slipperiness time-independent mask at PSI-points:        !
# !               (0=Land, 1=Sea, 2=no-slip).                            !
# !  rmask      Time-independent mask at RHO-points (0=Land, 1=Sea).     !
# !  umask      Time-independent mask at U-points (0=Land, 1=Sea).       !
# !  vmask      Time-independent mask at V-points (0=Land, 1=Sea).       !
# !                                                                      !
# !  pmask_full Full mask at PSI-points (0=dry, 1=wet, 2=no-slip).       !
# !  rmask_full Full mask at RHO-points (0=dry, 1=wet).                  !
# !  rmask_full Full mask at   U-points (0=dry, 1=wet).                  !
# !  rmask_full Full mask at   V-points (0=dry, 1=wet).                  !

# !  xp         XI-coordinates (m) at PSI-points.                        !
# !  xr         XI-coordinates (m) at RHO-points.                        !
# !  xu         XI-coordinates (m) at U-points.                          !
# !  xv         XI-coordinates (m) at V-points.                          !
# !  yp         ETA-coordinates (m) at PSI-points.                       !
# !  yr         ETA-coordinates (m) at RHO-points.                       !
# !  yu         ETA-coordinates (m) at U-points.                         !
# !  yv         ETA-coordinates (m) at V-points.                         !
# !  zice       Depth of ice shelf cavity (m, negative) at               !
# !               RHO-points.                                            !
# !  z0_r       Time independent depths (m) at horizontal RHO-points and !
# !               vertical RHO-points.                                   !
# !  z0_w       Time independent depths (m) at horizontal RHO-points and !
# !               vertical W-points.                                     !
# !  z_r        Actual depths (m) at horizontal RHO-points and           !
# !               vertical RHO-points.                                   !
# !  z_w        Actual depths (m) at horizontal RHO-points and           !
# !               vertical W-points.                                     !
# !                                                                      !

class Grid:

    def __init__(self, LBi, UBi, LBj, UBj, LBij, UBij):
        
        # Set horizontal array size.
        size2D = (UBi-LBi+1)*(UBj-LBj+1)
    
        shape2D = (UBi - LBi + 1, UBj - LBj)
    
        # Nonlinear model state.
        self.angler    = cp.zeros(shape2D, dtype = cp.float64)
        self.CosAngler = cp.zeros(shape2D, dtype = cp.float64)
        self.SinAngler = cp.zeros(shape2D, dtype = cp.float64)
    
    
        if CURVGRID && UV_ADV:
            self.dmde = cp.zeros(shape2D, dtype = cp.float64)
            self.dndx = cp.zeros(shape2D, dtype = cp.float64)
    
        self.f      = cp.zeros(shape2D, dtype = cp.float64)
        self.fomn   = cp.zeros(shape2D, dtype = cp.float64)
        self.grdscl = cp.zeros(shape2D, dtype = cp.float64)
        self.h      = cp.zeros(shape2D, dtype = cp.float64)
        self.latp   = cp.zeros(shape2D, dtype = cp.float64)
        self.latr   = cp.zeros(shape2D, dtype = cp.float64)
        self.latu   = cp.zeros(shape2D, dtype = cp.float64)
        self.latv   = cp.zeros(shape2D, dtype = cp.float64)
        self.lonp   = cp.zeros(shape2D, dtype = cp.float64)
        self.lonr   = cp.zeros(shape2D, dtype = cp.float64)
        self.lonu   = cp.zeros(shape2D, dtype = cp.float64)
        self.lonv   = cp.zeros(shape2D, dtype = cp.float64)
        self.Mylon  = cp.zeros(shape2D, dtype = cp.float64)
        self.omn    = cp.zeros(shape2D, dtype = cp.float64)
        self.om_p   = cp.zeros(shape2D, dtype = cp.float64)
        self.om_r   = cp.zeros(shape2D, dtype = cp.float64)
        self.om_u   = cp.zeros(shape2D, dtype = cp.float64)
        self.om_v   = cp.zeros(shape2D, dtype = cp.float64)
        self.on_p   = cp.zeros(shape2D, dtype = cp.float64)
        self.on_r   = cp.zeros(shape2D, dtype = cp.float64)
        self.on_u   = cp.zeros(shape2D, dtype = cp.float64)
        self.on_v   = cp.zeros(shape2D, dtype = cp.float64)
        self.pm     = cp.zeros(shape2D, dtype = cp.float64)
        self.pn     = cp.zeros(shape2D, dtype = cp.float64)
        self.pmon_p = cp.zeros(shape2D, dtype = cp.float64)
        self.pmon_r = cp.zeros(shape2D, dtype = cp.float64)
        self.pmon_u = cp.zeros(shape2D, dtype = cp.float64)
        self.pmon_v = cp.zeros(shape2D, dtype = cp.float64)
        self.pnom_p = cp.zeros(shape2D, dtype = cp.float64)
        self.pnom_r = cp.zeros(shape2D, dtype = cp.float64)
        self.pnom_u = cp.zeros(shape2D, dtype = cp.float64)
        self.pnom_v = cp.zeros(shape2D, dtype = cp.float64)
        self.xp     = cp.zeros(shape2D, dtype = cp.float64)
        self.xr     = cp.zeros(shape2D, dtype = cp.float64)
        self.xu     = cp.zeros(shape2D, dtype = cp.float64)
        self.xv     = cp.zeros(shape2D, dtype = cp.float64)
        self.yp     = cp.zeros(shape2D, dtype = cp.float64)
        self.yr     = cp.zeros(shape2D, dtype = cp.float64)
        self.yu     = cp.zeros(shape2D, dtype = cp.float64)
        self.yv     = cp.zeros(shape2D, dtype = cp.float64)
        self.pmask  = cp.zeros(shape2D, dtype = cp.bool)
        self.rmask  = cp.zeros(shape2D, dtype = cp.bool)
        self.umask  = cp.zeros(shape2D, dtype = cp.bool)
        self.vmask  = cp.zeros(shape2D, dtype = cp.bool)
        self.pmask_full = cp.zeros(shape2D, dtype = cp.bool)
        self.rmask_full = cp.zeros(shape2D, dtype = cp.bool)
        self.umask_full = cp.zeros(shape2D, dtype = cp.bool)
        self.vmask_full = cp.zeros(shape2D, dtype = cp.bool)


  