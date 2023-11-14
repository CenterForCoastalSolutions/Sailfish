from misc import *
# from get_data import get_data
# from set_data import set_data
from step2d import step2dPredictor, step2dCorrector
from mod_operators import *

import matplotlib.pyplot as plt


def output():
    pass

def main3d(compTimes, GRID, OCEAN, BOUNDARY):
    """ This subroutine is the main driver for nonlinear ROMS/TOMS when configurated as a full 3D baroclinic ocean
    model only.
    """

    VerticalVelEq(velEq, BC, FC, RHS)

    # Time-step nonlinear 3D primitive equations by the specified time.
    # =======================================================================



    # Initialize all time levels and compute other initial fields.
    # -----------------------------------------------------------------------

    if compTimes.isInitialTimeStep:
        # Initialize free-surface and compute initial level thicknesses and depths.
        ini_zeta()
        set_depth(GRID.Vtransform, GRID.Zt_avg1, GRID.z, GRID.z_r, GRID.z_w, GRID.h, GRID.hc, GRID.Hz, GRID.sc_r,  GRID.sc_w, GRID.Cs_r, GRID.Cs_w)

        # Initialize other state variables.
        ini_fields()


    while not compTimes.isFinalTimeStep():

        # Notice that RunInterval is set in the calling driver. Its value may span the full period of the
        # simulation, a multi-model coupling interval (RunInterval > ifac*dt), or just a single step (RunInterval = 0).


        # Time-step governing equations for Nsteps.
        for istep in range(compTimes.Nsteps):

        # Set time indices and time clock.

            # Gets next time step and cycles variables.
            # TODO: HERE OR BELOW????

            OCEAN.cycleTimes3D()

#         iic += 1
#         nstp(ng)=1+MOD(iic(ng)-ntstart(ng),2)
#         nnew(ng)=3-nstp(ng)
#         nrhs(ng)=nstp(ng)
# !       time(ng)=time(ng)+dt(ng)
#               tdays(ng)=time(ng)*sec2day
#               IF (step_counter(ng).eq.Rsteps) Time_Step=.FALSE.


            # Read in required data, if any, from input NetCDF files.
            # ----------------------------------------------------------------------
            # get_data()

            # If applicable, process input data: time interpolate between data snapshots.
            # -----------------------------------------------------------------------
            # set_data()



            # Compute horizontal mass fluxes (Hz*u/n and Hz*v/m), density related
            # quatities and report global diagnostics.
            # -----------------------------------------------------------------------
            set_maxflux(u, v, Huon, Hvom)



            # Set fields for vertical boundary conditions. Process tidal forcing,
            # if any.
            # ---------------------------------------------------------------------


            set_vbc()


# Compute time-dependent vertical/horizontal mixing coefficients for
# momentum and tracers. Compute S-coordinate vertical velocity,
# diagnostically from horizontal mass divergence.
# -----------------------------------------------------------------------

            # XXX todo: put other mixings.
            ana_vmix()


            omega(BC)
            # bc_w3d()

            # This is needed only for output. I will not implement it just yet.
            # wvelocity (nstp)
            # bc_w3d()  # TODO: twice, really?


            # Set free-surface to its time-averaged value.  If applicable, accumulate time-averaged output data which
            # needs a irreversible loop in shared-memory jobs.
            set_zeta(zeta[1,:,:], zeta[2,:,:], ZtAvg1)

            # If appropriate, write out fields into output NetCDF files.  Notice that IO data is written in delayed
            # and serial mode.  Exit if last time step.
            output()


            # Compute right-hand-side terms for 3D equations.
            rhs3d()


            # Solve the vertically integrated primitive equations for the
            # free-surface and barotropic momentum components.
            # -----------------------------------------------------------------------

            for iif in range(compTimes.nfast + 1):

                # Set time indices for predictor step. The PREDICTOR_2D_STEP switch
                # it is assumed to be false before the first time-step.
                compTimes.nextTimeStep()
                OCEAN.cycleTimes2D()

                Δt = compTimes.dtfast

                step2dPredictor(compTimes, GRID, OCEAN, BOUNDARY)

                if compTimes.iic % 500 == 0:
                    plt.clf()
                    plt.imshow(OCEAN.zeta_t2.reshape(GRID.M+1, GRID.M+1).get())
                    plt.colorbar()
                    plt.pause(1)
                    print('.... %.2f s    %.2f' % (compTimes.time, OCEAN.zeta_t2.sum()))

                # Corrector step - Apply 2D time-step corrector scheme.  XXXX-> Notice that there is not need for a corrector step
                # during the auxiliary (nfast+1) time-step.
                step2dCorrector(compTimes, GRID, OCEAN, BOUNDARY)

     # TODO: THIS*****
     #            # Predictor step - Advance barotropic equations using 2D time-step
     #            # ==============   predictor scheme.  No actual time-stepping is
     #            # performed during the auxiliary (nfast+1) time-step. It is needed
     #            # to finalize the fast-time averaging of 2D fields, if any, and
     #            # compute the new time-evolving depths.
     #
     #            IF (my_iif.le.(nfast(ng)+1)) THEN
     #                CALL step2d()
     #            END IF
     #




            # Recompute depths and thicknesses using the new time filtered free-surface.
            set_depth(GRID.Vtransform, GRID.Zt_avg1, GRID.z, GRID.z_r, GRID.z_w, GRID.h, GRID.hc, GRID.Hz, GRID.sc_r,  GRID.sc_w, GRID.Cs_r, GRID.Cs_w)


            # Time-step 3D momentum equations and couple with vertically integrated equations.
            step3d_uv(OCEAN.u_t2, OCEAN.v_t2, OCEAN.ru_t1, OCEAN.rv_t1, VerticalVelEq)


            setLateralUVBCs(OCEAN.u_t2, OCEAN.v_t2)

            #   Time-step vertical mixing turbulent equations and passive tracer
            #   source and sink terms, if applicable.
            omega(BC)
            # bc_w3d()

            cycleTimes3D()
            compTimes.nextTimeStep()


