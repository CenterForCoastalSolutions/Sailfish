#!/usr/bin/python
# -*- coding: utf-8 -*-


import cupy as cp
# import mod_param

#   IJwaterR   Water points IJ couter for RHO-points masked variables.
#   IJwaterU   Water points IJ couter for   U-points masked variables.
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

class Grid:

    def __init__(self, input):

        # Read sizes from the input file.
        self.Lm = input.getVal('Lm', dtype = int, minVal = 2)
        self.Mm = input.getVal('Mm', dtype = int, minVal = 2)
        self.L = self.Lm + 1
        self.M = self.Mm + 1

        # Set horizontal array size.
        shape2D = (self.L + 1, self.M + 1)
    
        # Nonlinear model state.
        self.angler    = cp.zeros(shape2D, dtype = cp.float64)
        self.CosAngler = cp.zeros(shape2D, dtype = cp.float64)
        self.SinAngler = cp.zeros(shape2D, dtype = cp.float64)


        self.Hz     = cp.zeros(shape2D, dtype = cp.float64)
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

        # if CURVGRID && UV_ADV: TODO: Implement the if
        self.dmdη   = cp.zeros(shape2D, dtype = cp.float64)
        self.dndξ   = cp.zeros(shape2D, dtype = cp.float64)

        self.pmon_p = cp.zeros(shape2D, dtype = cp.float64)
        self.pmon_r = cp.zeros(shape2D, dtype = cp.float64)
        self.pmon_u = cp.zeros(shape2D, dtype = cp.float64)
        self.pmon_v = cp.zeros(shape2D, dtype = cp.float64)
        self.pnom_p = cp.zeros(shape2D, dtype = cp.float64)
        self.pnom_r = cp.zeros(shape2D, dtype = cp.float64)
        self.pnom_u = cp.zeros(shape2D, dtype = cp.float64)
        self.pnom_v = cp.zeros(shape2D, dtype = cp.float64)
        # self.xp     = cp.zeros(shape2D, dtype = cp.float64)
        # self.xr     = cp.zeros(shape2D, dtype = cp.float64)
        # self.xu     = cp.zeros(shape2D, dtype = cp.float64)
        # self.xv     = cp.zeros(shape2D, dtype = cp.float64)
        # self.yp     = cp.zeros(shape2D, dtype = cp.float64)
        # self.yr     = cp.zeros(shape2D, dtype = cp.float64)
        # self.yu     = cp.zeros(shape2D, dtype = cp.float64)
        # self.yv     = cp.zeros(shape2D, dtype = cp.float64)
        self.pmask  = cp.zeros(shape2D, dtype = bool)
        self.rmask  = cp.zeros(shape2D, dtype = bool)
        self.umask  = cp.zeros(shape2D, dtype = bool)
        self.vmask  = cp.zeros(shape2D, dtype = bool)
        self.pmask_full = cp.zeros(shape2D, dtype = bool)
        self.rmask_full = cp.zeros(shape2D, dtype = bool)
        self.umask_full = cp.zeros(shape2D, dtype = bool)
        self.vmask_full = cp.zeros(shape2D, dtype = bool)


    def updateMetrics(self):
        """This routine computes various horizontal metric terms."""

        # Compute 1/m, 1/n, 1/mn, and f/mn at horizontal RHO-points.
        # -----------------------------------------------------------------------


        self.om_r = 1.0/self.pm
        self.on_r = 1.0/self.pn
        self.omn  = 1.0/(self.pm*self.pn)
        self.fomn = self.f*self.omn


        # Compute n/m, and m/n at horizontal RHO-points.
        # -----------------------------------------------------------------------

        self.pnom_r = self.pn/self.pm
        self.pmon_r = self.pm/self.pn


        # Compute m/n, 1/m, and 1/n at horizontal U-points.
        # -----------------------------------------------------------------------

        self.pmon_u = (self.pm[:-1,:] + self.pm[1:,:])/(self.pn[:-1,:] + self.pn[1:,:])
        self.pnom_u = (self.pn[:-1,:] + self.pn[1:,:])/(self.pm[:-1,:] + self.pm[1:,:])
        self.om_u[1:,:] = 2.0/(self.pm[:-1,:] + self.pm[1:,:])
        self.on_u[1:,:] = 2.0/(self.pn[:-1,:] + self.pn[1:,:])



        # Compute m/n, 1/m, and 1/n at horizontal U-points.
        # -----------------------------------------------------------------------

        self.pmon_v = (self.pm[:,:-1] + self.pm[:,1:]) / (self.pn[:,:-1] + self.pn[:,1:])
        self.pnom_v = (self.pn[:,:-1] + self.pn[:,1:]) / (self.pm[:,:-1] + self.pm[:,1:])
        self.om_v[:,1:] = 2.0 / (self.pm[:,:-1] + self.pm[:,1:])
        self.on_v[:,1:] = 2.0 / (self.pn[:,:-1] + self.pn[:,1:])


        # Compute m/n, 1/m, and 1/n at horizontal PSI-points.
        # -----------------------------------------------------------------------

        self.pmon_p = (self.pm[:-1,:-1] + self.pm[:-1,1:] + self.pm[1:,:-1] + self.pm[1:,1:]) / (self.pn[:-1,:-1] + self.pn[:-1,1:] + self.pn[1:,:-1] + self.pn[1:,1:])
        self.pnom_p = (self.pn[:-1,:-1] + self.pn[:-1,1:] + self.pn[1:,:-1] + self.pn[1:,1:]) / (self.pm[:-1,:-1] + self.pm[:-1,1:] + self.pm[1:,:-1] + self.pm[1:,1:])
        self.om_p = 4.0 / (self.pm[:-1,:-1] + self.pm[:-1,1:] + self.pm[1:,:-1] + self.pm[1:,1:])
        self.on_p = 4.0 / (self.pn[:-1,:-1] + self.pn[:-1,1:] + self.pn[1:,:-1] + self.pn[1:,1:])




    # #ifdef MASKING
    # !
    # !-----------------------------------------------------------------------
    # !  Set slipperiness (no-slip) mask at PSI-points.
    # !-----------------------------------------------------------------------
    # !
    # ! Set no-slip boundary conditions on land-mask boundaries regardless of
    # ! supplied value of gamma2.
    # !
    #       cff1=1.0_r8       ! computation of off-diagonal nonlinear terms
    #       cff2=2.0_r8
    #       DO j=JstrP,JendP
    #         DO i=IstrP,IendP
    #           IF ((rmask(i-1,j  ).gt.0.5_r8).and.                           &
    #      &        (rmask(i  ,j  ).gt.0.5_r8).and.                           &
    #      &        (rmask(i-1,j-1).gt.0.5_r8).and.                           &
    #      &        (rmask(i  ,j-1).gt.0.5_r8)) THEN
    #             pmask(i,j)=1.0_r8
    #           ELSE IF ((rmask(i-1,j  ).lt.0.5_r8).and.                      &
    #      &             (rmask(i  ,j  ).gt.0.5_r8).and.                      &
    #      &             (rmask(i-1,j-1).gt.0.5_r8).and.                      &
    #      &             (rmask(i  ,j-1).gt.0.5_r8)) THEN
    #             pmask(i,j)=cff1
    #           ELSE IF ((rmask(i-1,j  ).gt.0.5_r8).and.                      &
    #      &             (rmask(i  ,j  ).lt.0.5_r8).and.                      &
    #      &             (rmask(i-1,j-1).gt.0.5_r8).and.                      &
    #      &             (rmask(i  ,j-1).gt.0.5_r8)) THEN
    #             pmask(i,j)=cff1
    #           ELSE IF ((rmask(i-1,j  ).gt.0.5_r8).and.                      &
    #      &             (rmask(i  ,j  ).gt.0.5_r8).and.                      &
    #      &             (rmask(i-1,j-1).lt.0.5_r8).and.                      &
    #      &             (rmask(i  ,j-1).gt.0.5_r8)) THEN
    #             pmask(i,j)=cff1
    #           ELSE IF ((rmask(i-1,j  ).gt.0.5_r8).and.                      &
    #      &             (rmask(i  ,j  ).gt.0.5_r8).and.                      &
    #      &             (rmask(i-1,j-1).gt.0.5_r8).and.                      &
    #      &             (rmask(i  ,j-1).lt.0.5_r8)) THEN
    #             pmask(i,j)=cff1
    #           ELSE IF ((rmask(i-1,j  ).gt.0.5_r8).and.                      &
    #      &             (rmask(i  ,j  ).lt.0.5_r8).and.                      &
    #      &             (rmask(i-1,j-1).gt.0.5_r8).and.                      &
    #      &             (rmask(i  ,j-1).lt.0.5_r8)) THEN
    #             pmask(i,j)=cff2
    #           ELSE IF ((rmask(i-1,j  ).lt.0.5_r8).and.                      &
    #      &             (rmask(i  ,j  ).gt.0.5_r8).and.                      &
    #      &             (rmask(i-1,j-1).lt.0.5_r8).and.                      &
    #      &             (rmask(i  ,j-1).gt.0.5_r8)) THEN
    #             pmask(i,j)=cff2
    #           ELSE IF ((rmask(i-1,j  ).gt.0.5_r8).and.                      &
    #      &             (rmask(i  ,j  ).gt.0.5_r8).and.                      &
    #      &             (rmask(i-1,j-1).lt.0.5_r8).and.                      &
    #      &             (rmask(i  ,j-1).lt.0.5_r8)) THEN
    #             pmask(i,j)=cff2
    #           ELSE IF ((rmask(i-1,j  ).lt.0.5_r8).and.                      &
    #      &             (rmask(i  ,j  ).lt.0.5_r8).and.                      &
    #      &             (rmask(i-1,j-1).gt.0.5_r8).and.                      &
    #      &             (rmask(i  ,j-1).gt.0.5_r8)) THEN
    #             pmask(i,j)=cff2
    #           ELSE
    #             pmask(i,j)=0.0_r8
    #           END IF
    #         END DO
    #       END DO





        # Compute cosine and sine of grid rotation angle.
        # -----------------------------------------------------------------------

        self.CosAngler = cp.cos(self.angler)
        self.SinAngler = cp.sin(self.angler)



        # Compute minimum and maximum grid spacing.
        # -----------------------------------------------------------------------

        # Compute grid spacing range.


        return      # TODO: remove this


        rmask = self.rmask
        # dtfast = compTimes.dtfast

        grd = cp.sqrt(self.om_r[rmask]*self.on_r[rmask])

        grdmax = cp.max(grd)


        # my_DXmin= Large
        #       my_DXmax=-Large
        #       my_DYmin= Large
        #       my_DYmax=-Large
        # #ifdef MASKING
        #       my_DXminW= Large
        #       my_DXmaxW=-Large
        #       my_DYminW= Large
        #       my_DYmaxW=-Large
        # #endif
        # #ifdef SOLVE3D
        #       my_DZmin= Large
        #       my_DZmax=-Large
        # # ifdef MASKING
        #       my_DZminW= Large
        #       my_DZmaxW=-Large
        # # endif
        # #endif
        #       my_grdmax=-Large
        #       DO j=JstrT,JendT
        #         DO i=IstrT,IendT
        #
        #           cff=SQRT(om_r(i,j)*on_r(i,j))
        #
        #           my_DXmin=MIN(my_DXmin,om_r(i,j))
        #           my_DXmax=MAX(my_DXmax,om_r(i,j))
        #           my_DYmin=MIN(my_DYmin,on_r(i,j))
        #           my_DYmax=MAX(my_DYmax,on_r(i,j))
        # #ifdef MASKING
        #           IF (rmask(i,j).gt.0.0_r8) THEN
        #             my_grdmax=MAX(my_grdmax,cff)
        #             my_DXminW=MIN(my_DXminW,om_r(i,j))
        #             my_DXmaxW=MAX(my_DXmaxW,om_r(i,j))
        #             my_DYminW=MIN(my_DYminW,on_r(i,j))
        #             my_DYmaxW=MAX(my_DYmaxW,on_r(i,j))
        #           END IF
        # #else
        #           my_grdmax=MAX(my_grdmax,cff)
        # #endif



        # Compute gravity waves Courant number.
        # -----------------------------------------------------------------------

        # The 2D Courant number is defined as:

        #     Cg = c * dt * SQRT (1/dx^2 + 1/dy^2)

        # where c=SQRT(g*h) is gravity wave speed, and dx, dy are grid spacing in each direction.

        # Cg = dtfast * cp.sqrt(g * cp.abs(h) * (self.pm*self.pm + self.pn*self.pn))
        #
        # Cg_min = cp.min(Cg[rmask])
        # Cg_max = cp.max(Cg[rmask])
        #
        # Cg_Coriolis_min = dt*cp.min(cp.abs(f[rmask]))
        # Cg_Coriolis_max = dt*cp.max(cp.abs(f[rmask]))
        #
        # Cg_min = min(Cg_min, Cg_Coriolis_min)
        # Cg_max = max(Cg_max, Cg_Coriolis_max)


    def printReport(self):
        pass # XXXXX implement this


    #
    #
    # !$OMP CRITICAL (REDUCTIONS)
    #       Cg_min(ng)=MIN(Cg_min(ng),my_Cg_min)
    #       Cg_max(ng)=MAX(Cg_max(ng),my_Cg_max)
    #       Cg_Cor(ng)=MAX(Cg_Cor(ng),my_Cg_Cor)
    #       grdmax(ng)=MAX(grdmax(ng),my_grdmax)
    #       DXmin(ng)=MIN(DXmin(ng),my_DXmin)
    #       DXmax(ng)=MAX(DXmax(ng),my_DXmax)
    #       DYmin(ng)=MIN(DYmin(ng),my_DYmin)
    #       DYmax(ng)=MAX(DYmax(ng),my_DYmax)
    # #ifdef MASKING
    #       DXminW(ng)=MIN(DXminW(ng),my_DXminW)
    #       DXmaxW(ng)=MAX(DXmaxW(ng),my_DXmaxW)
    #       DYminW(ng)=MIN(DYminW(ng),my_DYminW)
    #       DYmaxW(ng)=MAX(DYmaxW(ng),my_DYmaxW)
    #endif
    #         buffer(9)=0.0_r8
    #         op_handle(9)='MIN'
    #         buffer(10)=0.0_r8
    #         op_handle(10)='MAX'
    #         buffer(11)=0.0_dp
    #         op_handle(11)='MIN'
    #         buffer(12)=0.0_dp
    #         op_handle(12)='MAX'
    #         buffer(13)=0.0_dp
    #         op_handle(13)='MIN'
    #         buffer(14)=0.0_dp
    #         op_handle(14)='MAX'
    # # endif
    # # ifdef MASKING
    #         buffer(15)=DXminW(ng)
    #         op_handle(15)='MIN'
    #         buffer(16)=DXmaxW(ng)
    #         op_handle(16)='MAX'
    #         buffer(17)=DYminW(ng)
    #         op_handle(17)='MIN'
    #         buffer(18)=DYmaxW(ng)
    #         op_handle(18)='MAX'
    # #  ifdef SOLVE3D
    #         buffer(19)=DZminW(ng)
    #         op_handle(19)='MIN'
    #         buffer(20)=DZmaxW(ng)
    #         op_handle(20)='MAX'
    # #  else
    #         buffer(19)=0.0_dp
    #         op_handle(19)='MIN'
    #         buffer(20)=0.0_dp
    #         op_handle(20)='MAX'
    # #  endif
    #         CALL mp_reduce (ng, model, 20, buffer, op_handle)
    # # else
    #         CALL mp_reduce (ng, model, 14, buffer, op_handle)
    # # endif
    #         Cg_min(ng)=buffer(1)
    #         Cg_max(ng)=buffer(2)
    #         Cg_Cor(ng)=buffer(3)
    #         grdmax(ng)=buffer(4)
    #         DXmin(ng)=buffer(5)
    #         DXmax(ng)=buffer(6)
    #         DYmin(ng)=buffer(7)
    #         DYmax(ng)=buffer(8)
    # # ifdef SOLVE3D
    #         DZmin(ng)=buffer(9)
    #         DZmax(ng)=buffer(10)
    # # endif
    # # ifdef DIFF_3DCOEF
    #         DiffMin(ng)=buffer(11)
    #         DiffMax(ng)=buffer(12)
    # # endif
    # # ifdef VISC_3DCOEF
    #         ViscMin(ng)=buffer(13)
    #         ViscMax(ng)=buffer(14)
    # # endif
    # # ifdef MASKING
    #         DXminW(ng)=buffer(15)
    #         DXmaxW(ng)=buffer(16)
    #         DYminW(ng)=buffer(17)
    #         DYmaxW(ng)=buffer(18)
    #
    # # endif
    # #endif
    #         IF (Master.and.LwrtInfo(ng)) THEN
    #           WRITE (stdout,10) ng,                                         &
    # #ifdef MASKING
    #      &                      DXmin(ng)*0.001_dp, DXminW(ng)*0.001_dp,    &
    #      &                      DXmax(ng)*0.001_dp, DXmaxW(ng)*0.001_dp,    &
    #      &                      DYmin(ng)*0.001_dp, DYminW(ng)*0.001_dp,    &
    #      &                      DYmax(ng)*0.001_dp, DYmaxW(ng)*0.001_dp
    #   10      FORMAT (/,' Metrics information for Grid ',i2.2,':',          &
    #      &            /,' ===============================',/,               &
    #      &            /,' Minimum X-grid spacing, DXmin = ',1pe15.8,' km',  &
    #      &              4x,'Water points = ',1pe15.8,' km',                 &
    #      &            /,' Maximum X-grid spacing, DXmax = ',1pe15.8,' km',  &
    #      &              4x,'Water points = ',1pe15.8,' km',                 &
    #      &            /,' Minimum Y-grid spacing, DYmin = ',1pe15.8,' km',  &
    #      &              4x,'Water points = ',1pe15.8,' km',                 &
    #      &            /,' Maximum Y-grid spacing, DYmax = ',1pe15.8,' km',  &
    #      &              4x,'Water points = ',1pe15.8,' km')
    # #else
    #      &                      DXmin(ng)*0.001_dp,                         &
    #      &                      DXmax(ng)*0.001_dp,                         &
    #      &                      DYmin(ng)*0.001_dp,                         &
    #      &                      DYmax(ng)*0.001_dp
    #   10      FORMAT (/,' Metrics information for Grid ',i2.2,':',          &
    #      &            /,' ===============================',/,               &
    #      &            /,' Minimum X-grid spacing, DXmin = ',1pe15.8,' km',  &
    #      &            /,' Maximum X-grid spacing, DXmax = ',1pe15.8,' km',  &
    #      &            /,' Minimum Y-grid spacing, DYmin = ',1pe15.8,' km',  &
    #      &            /,' Maximum Y-grid spacing, DYmax = ',1pe15.8,' km')
    # #endif
    # #ifdef SOLVE3D
    # # ifdef MASKING
    #           WRITE (stdout,20) DZmin(ng), DZminW(ng), DZmax(ng), DZmaxW(ng)
    #   20      FORMAT (' Minimum Z-grid spacing, DZmin = ',1pe15.8,' m',     &
    #      &            5x,'Water points = ',1pe15.8,' m',/,                  &
    #      &            ' Maximum Z-grid spacing, DZmax = ',1pe15.8,' m',     &
    #      &            5x,'Water points = ',1pe15.8,' m')
    # # else
    #           WRITE (stdout,20) DZmin(ng), DZmax(ng)
    #   20      FORMAT (' Minimum Z-grid spacing, DZmin = ',1pe15.8,' m',/,   &
    #      &            ' Maximum Z-grid spacing, DZmax = ',1pe15.8,' m')
    # # endif
    # #endif
    #           WRITE (stdout,30) Cg_min(ng), Cg_max(ng), Cg_Cor(ng)
    #   30      FORMAT (/,' Minimum barotropic Courant Number = ', 1pe15.8,/, &
    #      &              ' Maximum barotropic Courant Number = ', 1pe15.8,/, &
    #      &              ' Maximum Coriolis   Courant Number = ', 1pe15.8,/)
    # # if defined VISC_GRID || defined DIFF_GRID
    #           WRITE (stdout,40) grdmax(ng)/1000.0_r8
    #   40      FORMAT (' Horizontal mixing scaled by grid area squared root',&
    # #  ifdef MASKING
    #      &            ', MAXVAL(grdscl) = ',1pe15.8,' km',                  &
    #      &            2x,'(Water points)')
    # #  else
    #      &            ', MAXVAL(grdscl) = ',1pe15.8,' km')
    # #  endif
    # # endif
    # # ifdef DIFF_3DCOEF
    # #  ifdef TS_DIF2
    #           units='m2/s'
    # #  elif defined TS_DIF4
    #           units='m4/s'
    # #  endif
    #           WRITE (stdout,50) DiffMin(ng), TRIM(units),                   &
    #      &                      DiffMax(ng), TRIM(units)
    #   50      FORMAT (/,' Minimum horizontal diffusion coefficient = ',     &
    #      &             1pe15.8,1x,a,                                        &
    #      &            /,' Maximum horizontal diffusion coefficient = ',     &
    #      &             1pe15.8,1x,a)
    # # endif
    # # ifdef VISC_3DCOEF
    # #  ifdef UV_VIS2
    #           units='m2/s'
    # #  elif defined UV_VIS4
    #           units='m4/s'
    # #  endif
    #           WRITE (stdout,60) ViscMin(ng), TRIM(units),                   &
    #      &                      ViscMax(ng), TRIM(units)
    #   60      FORMAT (/,' Minimum horizontal viscosity coefficient = ',     &
    #      &             1pe15.8,1x,a,                                        &
    #      &            /,' Maximum horizontal viscosity coefficient = ',     &
    #      &             1pe15.8,1x,a)
    # # endif
    #         END IF
    #       END IF





# Vertical grid structure. This structure used to be called T_SCALARS and been in mod_scalars, but I don't think it
# made sense there. In fact, I would argue is should be part of Grid.
# -----------------------------------------------------------------------
# Cs_r          Set of S-curves used to stretch the vertical grid that follows the bathymetry at vertical RHO-points.
# Cs_w          Set of S-curves used to stretch the vertical grid that follows the bathymetry at vertical W-points.
# sc_r          S-coordinate independent variable, [-1 < sc < 0] at vertical RHO-points.
# sc_w          S-coordinate independent variable, [-1 < sc < 0] at  vertical W-points.

class VertGrid:
    def __init__(self, input):

        N = input.getVal('N' , minVal = 0)               # Number of vertical levels.


        # sigma coordinates.
        Vtransform  = input.getVal('Vtransform',  validValues = (1, 2), dtype = int)
        Vstretching = input.getVal('Vstretching', validValues = (1, 2, 3, 4, 5), dtype=int)


        theta_s = input.getVal('THETA_S', dtype=float)       # S-coordinate surface control parameter. The range of optimal values depends on the vertical stretching function, C(s).
        theta_b = input.getVal('THETA_B', dtype=float)       # S-coordinate bottom  control parameter. The range of optimal values depends on the vertical stretching function, C(s).
        Tcline  = input.getVal('TCLINE',  dtype=float)       # Critical depth (hc) in meters (positive) controlling the stretching. It can be interpreted as the width of surface or
                                                            #   bottom boundary layer in which higher vertical resolution (levels) is required during stretching.
        hc = Tcline



#         # real(dp), pointer :: Cs_r(:)
#         # real(dp), pointer :: Cs_w(:)
#         # real(dp), pointer :: sc_r(:)
#         # real(dp), pointer :: sc_w(:)


  