from misc import *
# from get_data import get_data
# from set_data import set_data
from step2d import step2dPredictor, step2dCorrector

import matplotlib.pyplot as plt


def main2d(compTimes, GRID, OCEAN, BOUNDARY):
""" This subroutine is the main driver for nonlinear ROMS/TOMS when configurated as a full 3D baroclinic ocean
model only. It advances forward  the vertically integrated primitive equations by thespecified  time interval
(seconds), RunInterval.
"""


def main3d (RunInterval):
# This subroutine is the main driver for nonlinear ROMS/TOMS when     !
# configurated as a full 3D baroclinic ocean model.  It  advances     !
# forward the primitive equations for all  nested  grids, if any,     !
# for the specified time interval (seconds), RunInterval.             !


    # Time-step nonlinear 3D primitive equations by the specified time.
    # =======================================================================

    # Time-step the 3D kernel for the specified time interval (seconds),
    # RunInterval.

    Time_Step = True


    # Initialize all time levels and compute other initial fields.
    # -----------------------------------------------------------------------

    if (iic == ntstart):
        # Initialize free-surface and compute initial level thicknesses and depths.
        ini_zeta()
        set_depth()

        # Initialize other state variables.
        ini_fields()

    while Time_Step:


        # Notice that RunInterval is
        # set in the calling driver. Its value may span the full period of the
        # simulation, a multi-model coupling interval (RunInterval > ifac*dt),
        # or just a single step (RunInterval=0).

        ntimesteps (iNLM, RunInterval, nl, Nsteps, Rsteps)


        # Time-step governing equations for Nsteps.
        STEP_LOOP : DO istep=1,Nsteps

        # Set time indices and time clock.

            # Gets next time step and cycles variables.
        # TODO: HERE OR BELOW????
            compTimes.nextTimeStep()
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
            set_massflux()



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


            omega(iNLM)
            wvelocity (nstp)


            # Set free-surface to it time-averaged value.  If applicable, accumulate time-averaged output data which
            # needs a irreversible loop in shared-memory jobs.
            set_zeta()

            # If appropriate, write out fields into output NetCDF files.  Notice that IO data is written in delayed
            # and serial mode.  Exit if last time step.
            output()


            # Compute right-hand-side terms for 3D equations.
            rhs3d()


# ifdef GLS_MIXING
#             gls_prestep()
# endif


            # Solve the vertically integrated primitive equations for the
            # free-surface and barotropic momentum components.
            # -----------------------------------------------------------------------

            LOOP_2D : DO my_iif=1, MAXVAL(nfast)+1

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
     # #            next_indx1=3-indx1(ng)
     # #            IF (.not.PREDICTOR_2D_STEP(ng).and.                     &
     # # &              my_iif.le.(nfast(ng)+1)) THEN
     # #              PREDICTOR_2D_STEP(ng)=.TRUE.
     # #              iif(ng)=my_iif
     # #              IF (FIRST_2D_STEP) THEN
     # #                kstp(ng)=indx1(ng)
     # #              ELSE
     # #                kstp(ng)=3-indx1(ng)
     # #              END IF
     # #              knew(ng)=3
     # #              krhs(ng)=indx1(ng)
     # #            END IF
     #
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
     #
     #
     #         # Set time indices for corrector step.
     #            IF (PREDICTOR_2D_STEP(ng)) THEN
     #              PREDICTOR_2D_STEP(ng)=.FALSE.
     #              knew(ng)=next_indx1
     #              kstp(ng)=3-knew(ng)
     #              krhs(ng)=3
     #              IF (iif(ng).lt.(nfast(ng)+1)) indx1(ng)=next_indx1
     #            END IF
     #
     #        # Corrector step - Apply 2D time-step corrector scheme.  Notice that
     #        # ==============   there is not need for a corrector step during the
     #        # auxiliary (nfast+1) time-step.
     #            IF (iif(ng).lt.(nfast(ng)+1)) THEN
     #                CALL step2d ()
     #            END IF

            END DO LOOP_2D



            # Recompute depths and thicknesses using the new time filtered free-surface.
            set_depth()


            # Time-step 3D momentum equations and couple with vertically integrated equations.
            step3d_uv()


            #   Time-step vertical mixing turbulent equations and passive tracer
            #   source and sink terms, if applicable.
            omega(iNLM)
            # gls_corstep()


            # -----------------------------------------------------------------------
            #  Advance time index and time clock.
            # -----------------------------------------------------------------------
            iic(ng)=iic(ng)+1
            time(ng)=time(ng)+dt(ng)
            step_counter(ng)=step_counter(ng)-1



