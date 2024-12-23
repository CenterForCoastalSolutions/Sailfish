# from get_data import get_data
# from set_data import set_data
from rhs import rhs3d
from step2d import step2dPredictor, step2dCorrector
from mod_operators import step3d_UV, set_maxflux, omega, set_zeta, set_depth, setVerticalVelEq, setLateralUVBCs
import cupy as cp
from misc import *
import time

import matplotlib.pyplot as plt
import matplotlib as mpl
# mpl.use("TkAgg")

doSaveFile = False
doPlot = True


if doSaveFile:
    from netCDF4 import Dataset


def output():
    pass

def main3d(compTimes, GRID, OCEAN, BOUNDARY):
    """ This subroutine is the main driver for nonlinear ROMS/TOMS when configurated as a full 3D baroclinic ocean
    model only.
    """

    if doSaveFile:
        from adhocNetcdfOutput import create_roms_netcdf
        from netCDF4 import date2num
        from datetime import datetime, timedelta

        previousSaveTime = None
        outputFile = create_roms_netcdf('/blue/olabarrieta/jo.gonzalez/projects/sailfish/ROMS_tests/output.nc', GRID.L+1, GRID.M+1, GRID.N)
        outputFile['x_rho'][:,:] = GRID.lonr.get()
        outputFile['y_rho'][:,:] = GRID.latr.get()
        outputFile['x_u'  ][:,:] = GRID.lonu.get()[:,:-1]
        outputFile['y_u'  ][:,:] = GRID.latu.get()[:,:-1]
        outputFile['x_v'  ][:,:] = GRID.lonv.get()[:-1,:]
        outputFile['y_v'  ][:,:] = GRID.latv.get()[:-1,:]
        outputFile['h'    ][:,:] = GRID.h.get()
        outputFile['Cs_r' ][:]   = GRID.Cs_r.get()
        outputFile['s_rho'][:]   = GRID.sc_r.get()


    eqSD   = cp.zeros(GRID.shape3DW, dtype = cp.float64)
    eqD    = cp.zeros(GRID.shape3DW, dtype = cp.float64)
    eqRHS  = cp.zeros(GRID.shape3DW, dtype = cp.float64)
    tmpU   = cp.zeros(OCEAN.u_t2.shape, dtype = cp.float64)
    tmpV   = cp.zeros(OCEAN.v_t2.shape, dtype = cp.float64)
    from mod_operators import grsz, bksz

    BC = BOUNDARY.zetaBC.bcIdxFieldIdx2

    # Initialize the arrays for the matrix equations.
    setVerticalVelEq((1,), (1,), (eqSD, eqD, eqRHS))


    # Recompute depths and thicknesses using the new time filtered free-surface.
    set_depth(grsz, bksz, (GRID.Vtransform, OCEAN.Zt_avg1, GRID.z_w, GRID.z_r, GRID.h, GRID.hc, GRID.Hz,
                                 GRID.sc_r,  GRID.sc_w, GRID.Cs_r, GRID.Cs_w))


    # Time-step nonlinear 3D primitive equations by the specified time.
    # =======================================================================



    # Initialize all time levels and compute other initial fields.
    # -----------------------------------------------------------------------

    # if compTimes.isInitialTimeStep:
        # Initialize free-surface and compute initial level thicknesses and depths.

        # TODO: Write these.
        # ini_zeta()
        # set_depth(GRID.Vtransform, GRID.Zt_avg1, GRID.z, GRID.z_r, GRID.z_w, GRID.h, GRID.hc, GRID.Hz, GRID.sc_r,  GRID.sc_w, GRID.Cs_r, GRID.Cs_w)
        #
        # # Initialize other state variables.
        # ini_fields()

    imgIdx = 0
    elapsed2d = 0.0

    while not compTimes.isFinalTimeStep():

        # Notice that RunInterval is set in the calling driver. Its value may span the full period of the
        # simulation, a multi-model coupling interval (RunInterval > ifac*dt), or just a single step (RunInterval = 0).


        # Time-step governing equations for Nsteps.  TODO: CHECK, I changed Nsteps by ndtfast
        # for istep in range(compTimes.ndtfast):

        # Set time indices and time clock.

        # Gets next time step and cycles variables.
        # TODO: HERE OR BELOW????

        # OCEAN.cycleTimes3D()



        # Read in required data, if any, from input NetCDF files.
        # ----------------------------------------------------------------------
        # get_data()

        # If applicable, process input data: time interpolate between data snapshots.
        # -----------------------------------------------------------------------
        # set_data()



        # Compute horizontal mass fluxes (Hz*u/n and Hz*v/m), density related
        # quatities and report global diagnostics.
        # -----------------------------------------------------------------------
        set_maxflux(grsz, bksz, (OCEAN.u_t2, OCEAN.v_t2, OCEAN.Huon, OCEAN.Hvom, GRID.Hz))



        # Set fields for vertical boundary conditions. Process tidal forcing,
        # if any.
        # ---------------------------------------------------------------------

        # TODO: Recover!!
        # set_vbc()


# Compute time-dependent vertical/horizontal mixing coefficients for
# momentum and tracers. Compute S-coordinate vertical velocity,
# diagnostically from horizontal mass divergence.
# -----------------------------------------------------------------------

        # XXX todo: put other mixings.
        # TODO: Implement
        # ana_vmix()



        omega(grsz, bksz, (OCEAN.W, OCEAN.u_t2, OCEAN.v_t2, OCEAN.Huon.ravel(), OCEAN.Hvom.ravel(), GRID.z_w, BC))
        # bc_w3d()

        # This is needed only for output. I will not implement it just yet.
        # wvelocity (nstp)
        # bc_w3d()  # TODO: twice, really?


        # Set free-surface to its time-averaged value.  If applicable, accumulate time-averaged output data which
        # needs a irreversible loop in shared-memory jobs.
        set_zeta(grsz, bksz, (OCEAN.zeta_t1, OCEAN.zeta_t2, OCEAN.Zt_avg1))

        # If appropriate, write out fields into output NetCDF files.  Notice that IO data is written in delayed
        # and serial mode.  Exit if last time step.
        output()


        # Compute right-hand-side terms for 3D equations.
        rhs3d(GRID, OCEAN, BOUNDARY)


        # Solve the vertically integrated primitive equations for the
        # free-surface and barotropic momentum components.
        # -----------------------------------------------------------------------

        compTimes.first2DTimeStep()
        t2d_1 = time.time()
        for iif in range(compTimes.nfast+1):

            if compTimes.isFirst2DStep():
                OCEAN.Zt_avg1[:] = 0.0
                OCEAN.DU_avg1[:] = 0.0
                OCEAN.DV_avg1[:] = 0.0

            # Set time indices for predictor step. The PREDICTOR_2D_STEP switch
            # it is assumed to be false before the first time-step.
            compTimes.next2DTimeStep()
            OCEAN.cycleTimes2D()


            step2dPredictor(compTimes, GRID, OCEAN, BOUNDARY)


            # Corrector step - Apply 2D time-step corrector scheme.
            # Original Fortran code says: Notice that there is not need for a corrector step during the auxiliary (nfast+1) time-step.
            step2dCorrector(compTimes, GRID, OCEAN, BOUNDARY)

            # TODO: This seems to be necessary, but I'm not sure why. Maybe it is hiding a deeper bug?
            # OCEAN.vbar_t2[:]=0.0


        t2d_2 = time.time()
        elapsed2d += t2d_2 - t2d_1

 # TODO: I'm not sure if I have done this correctly.
 #            # Predictor step - Advance barotropic equations using 2D time-step
 #            # ==============   predictor scheme.  No actual time-stepping is
 #            # performed during the auxiliary (nfast+1) time-step. It is needed
 #            # to finalize the fast-time averaging of 2D fields, if any, and
 #            # compute the new time-evolving depths.

        if doSaveFile and (previousSaveTime is None or compTimes.time >= previousSaveTime + 30*60.0 -1e-4):
            if previousSaveTime is None:
                previousSaveTime = 0.0
            previousSaveTime += 30*60.0
            try:
                idxTime += 1
            except:
                idxTime = 0

            outputFile['ocean_time'][idxTime] = date2num(timedelta(seconds=compTimes.time) + datetime(2001,1,1), units = "seconds since 2001-01-01 00:00:00")
            outputFile['zeta'][idxTime,:,:] = OCEAN.zeta_t2.get()
            outputFile['u'][idxTime,:,:,:] = OCEAN.u_t2.get().reshape(GRID.N+1, GRID.M+1, GRID.L+1)[:-1,:,:-1]
            outputFile['v'][idxTime,:,:,:] = OCEAN.v_t2.get().reshape(GRID.N+1, GRID.M+1, GRID.L+1)[:-1,:-1,:]
            outputFile['ubar'][idxTime,:,:] = OCEAN.ubar_t2.get().reshape(GRID.M+1, GRID.L+1)[:,:-1]
            outputFile['vbar'][idxTime,:,:] = OCEAN.vbar_t2.get().reshape(GRID.M+1, GRID.L+1)[:-1,:]

        if doPlot and (compTimes.iic % 1000==0):
            # plt.close(True)
            try:
                fig.clf()
            except:
                # fig = plt.figure(figsize=(20,12))
                fig = plt.figure(figsize=(16,10))
            ax1 = fig.add_subplot(221)
            fig.suptitle('time = %.1f s' % compTimes.time)


            # cmap = mpl.cm.RdBu
            # norm = mpl.colors.Normalize(vmin=-.03, vmax=.03)
            # ax1.imshow(OCEAN.zeta_t2.reshape(GRID.M+1, GRID.M+1).get(), cmap = cmap, vmin = -.03, vmax = 0.03)
            # ax1.pcolormesh(GRID.xr[0,:].get()/1000, GRID.yr[:,0].get()/1000, OCEAN.zeta_t2.reshape(GRID.M+1, GRID.M+1).get(), cmap = cmap, vmin = -.03, vmax = 0.03)
            # ax1.title.set_text('Free surface (m)')
            # # plt.imshow(OCEAN.u_t2.get().reshape(17,402,402)[0,:,:])
            # fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),ax=ax1)
            # plt.pause(1)
            # print('.... %.2f s   ' % (compTimes.time))


            ax1.title.set_text('zeta (m) x-transect')
            ax1.plot(GRID.xp[0,2:].get(), -OCEAN.zeta_t2.reshape(GRID.M+1, GRID.M+1)[100,2:].get())
            ax1.set_ylim((-.13,0.13))

            ax2 = fig.add_subplot(222)
            cmap = mpl.cm.RdBu
            norm = mpl.colors.Normalize(vmin=-.03, vmax=.03)
            # ax2.imshow(OCEAN.ubar_t2.reshape(GRID.M+1, GRID.M+1).get(), cmap = cmap, vmin = -.02, vmax = 0.02)
            ax2.pcolormesh(GRID.xr[0,:].get()/1000, GRID.yr[:,0].get()/1000, OCEAN.ubar_t2.reshape(GRID.M+1, GRID.M+1).get(), cmap = cmap, vmin = -.03, vmax = 0.03)
            ax2.title.set_text('barotropic u (m/s)')
            # plt.imshow(OCEAN.u_t2.get().reshape(17,402,402)[0,:,:])
            fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),ax=ax2)


            ax3 = fig.add_subplot(223)
            ax3.title.set_text('barotropic u (m/s) x-transect')
            ax3.plot(GRID.xp[0,2:].get(), OCEAN.ubar_t2.reshape(GRID.M+1, GRID.M+1)[100,2:].get())
            ax3.set_ylim((-.1,0.1))

            ax4 = fig.add_subplot(224)
            ax4.title.set_text('u  (m/s) z-profile')
            ax4.plot(OCEAN.u_t2.reshape(GRID.N+1, GRID.M+1, GRID.M+1)[:-1,100,10].get(), GRID.z_w[:-1,0,0].get())
            ax4.plot(OCEAN.u_t2.reshape(GRID.N+1, GRID.M+1, GRID.M+1)[:-1,100,20].get(), GRID.z_w[:-1,0,0].get())
            ax4.plot(OCEAN.u_t2.reshape(GRID.N+1, GRID.M+1, GRID.M+1)[:-1,100,30].get(), GRID.z_w[:-1,0,0].get())
            ax4.plot(OCEAN.u_t2.reshape(GRID.N+1, GRID.M+1, GRID.M+1)[:-1,100,100].get(), GRID.z_w[:-1,0,0].get())
            ax4.plot(OCEAN.u_t2.reshape(GRID.N+1, GRID.M+1, GRID.M+1)[:-1,100,200].get(), GRID.z_w[:-1,0,0].get())
            ax4.plot(OCEAN.ubar_t2.reshape(GRID.M+1, GRID.M+1)[100,100].get() + 0*GRID.z_w[:-1,0,0].get(), GRID.z_w[:-1,0,0].get())
            ax4.legend(['x=%.3f' % (GRID.xr[0,10].get()/1000), 'x=%.3f' % (GRID.xr[0,20].get()/1000), 'x=%.3f' % (GRID.xr[0,30].get()/1000),
                                    'x=%.3f' % (GRID.xr[0,100].get()/1000), 'x=%.3f' % (GRID.xr[0,200].get()/1000)])

            ax4.set_xlim((-.05,0.05))

            imgIdx += 1

            # plt.subplot(1,2,1)
            # plt.imshow((OCEAN.u_t2).reshape(GRID.N+1,GRID.M+1,GRID.L+1)[-2,:190,:190].get())
            # plt.subplot(1,2,2)
            # plt.imshow((OCEAN.v_t2).reshape(GRID.N+1,GRID.M+1,GRID.L+1)[-2,:190,:190].get())

            # plt.savefig(r'D:\projects\src\oceangpu\img%.4i' % imgIdx)
            plt.savefig(r'img%.4i' % imgIdx)

            plt.pause(3)
            plt.ion()
            print('.... %.2f s    %.2f' % (compTimes.time, OCEAN.zeta_t2.sum()))




        # Recompute depths and thicknesses using the new time filtered free-surface.
        set_depth(grsz, bksz, (GRID.Vtransform, OCEAN.Zt_avg1, GRID.z_w, GRID.z_r, GRID.h, GRID.hc, GRID.Hz,
                               GRID.sc_r,  GRID.sc_w, GRID.Cs_r, GRID.Cs_w))


        # Time-step 3D momentum equations and couple with vertically integrated equations.
        λ = 1.0

        # TODO: Explain this better. This is to avoid a "sweeingp" error:
        tmpU[:] = OCEAN.u_t2
        tmpV[:] = OCEAN.v_t2



        # XXX: I had to reduce the number of threads because an error related to GPU's limited resources. This has to be done in a better way.
        step3d_UV((grsz[0]*GPUMUL,), (bksz[0]//GPUMUL,), (OCEAN.u_t2, OCEAN.v_t2, OCEAN.ru_t2, OCEAN.rv_t2, OCEAN.ubar_t1, OCEAN.vbar_t1,
                                                OCEAN.ubar_t2, OCEAN.vbar_t2, GRID.Hz, OCEAN.AKv, GRID.z_r, OCEAN.DU_avg1, OCEAN.DV_avg1, tmpU, tmpV, GRID.h,
                                                compTimes.iic, compTimes.ntfirst, λ, OCEAN.AK, compTimes.dt))




        setLateralUVBCs((grsz[0],), (bksz[0],), (compTimes.time, OCEAN.u_t2, OCEAN.v_t2, BOUNDARY.uvelBC.bcIdxFieldIdx2, BOUNDARY.vvelBC.bcIdxFieldIdx2, BOUNDARY.uvelBC.bcIdxFieldType, BOUNDARY.vvelBC.bcIdxFieldType))

        # OCEAN.v[:,:,0,:] = 0.0
        # OCEAN.vbar[:,0,:] = 0.0
        # OCEAN.u[:,:,0,:] = OCEAN.u[:,:,2,:]
        # OCEAN.ubar[:,0,:] = OCEAN.ubar[:,2,:]
        # OCEAN.u[:,:,1,:] = OCEAN.u[:,:,2,:]
        # OCEAN.ubar[:,1,:] = OCEAN.ubar[:,2,:]

        #   Time-step vertical mixing turbulent equations and passive tracer source and sink terms, if applicable.
        omega(grsz, bksz, (OCEAN.W, OCEAN.u_t2, OCEAN.v_t2, OCEAN.Huon, OCEAN.Hvom, GRID.z_w, BC))
        # bc_w3d()

        compTimes.nextTimeStep()
        print('time: %.2f s' % compTimes.time)

    print('TIME 2D ', elapsed2d)

    outputFile.close()