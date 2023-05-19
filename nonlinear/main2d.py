def main2d(RunInterval):
""" This subroutine is the main driver for nonlinear ROMS/TOMS when configurated as shallow water (barotropic) ocean
model only. It advances forward  the vertically integrated primitive equations for all nested grids,  if any,  by the
specified  time interval (seconds), RunInterval.
"""


    import mod_param
    import mod_scalars
    from ntimestep import ntimestep
    # import mod_stepping
    # import dateclock_mod
    # import step2d_mod
    #from ini_fields_mod import ini_zeta, ini_fields, set_vbc
    # import mod_iounits
    #use set_vbc_mod
    #use step_floats_mod
    # use strings_mod





    # Time-step vertically integrated equations.
    # =======================================================================
    # !
    # !  Time-step the 3D kernel for the specified time interval (seconds),
    # !  RunInterval.
    # !
    Time_Step = True

    # kernel loop
    while Time_Step:

        # Determine number of time steps to compute in each nested grid layer
        # based on the specified time interval (seconds), RunInterval. Non
        # nesting applications have NestLayers=1. Notice that RunInterval is
        # set in the calling driver. Its value may span the full period of the
        # simulation, a multi-model coupling interval (RunInterval > ifac*dt),
        # or just a single step (RunInterval=0).
        Nsteps, Rsteps = ntimesteps(iNLM, RunInterval, Nsteps, Rsteps)


        # Time-step governing equations for Nsteps.
        for istep in range(Nsteps):

            # Set time indices and time clock.
            iic += 1
            time += dt
            tdays = time*sec2day
            time_string(time, time_code)
            if (step_counter == Rsteps):
                Time_Step = False

            # Read in required data, if any, from input NetCDF files.
            # ----------------------------------------------------------------------
            get_data()

            # If applicable, process input data: time interpolate between data snapshots.
            # -----------------------------------------------------------------------
            set_data()


            # Initialize all time levels and compute other initial fields.
            # -----------------------------------------------------------------------
            if iic == ntstart:

                # Initialize free-surface.
                ini_zeta(iNLM)

                # Initialize other state variables.
                ini_fields (iNLM)



                # Set vertical boundary conditions. Process tidal forcing.
                # -----------------------------------------------------------------------
                set_vbc()



            # If appropriate, write out fields into output NetCDF files.  Notice
            # that IO data is written in delayed and serial mode.  Exit if last
            # time step.
            # -----------------------------------------------------------------------

            output()
            if (iic == ntend + 1):
                exit()


            # Solve the vertically integrated primitive equations for the
            # free-surface and momentum components.
            # -----------------------------------------------------------------------

            # Set time indices for predictor step. The PREDICTOR_2D_STEP switch
            # it is assumed to be false before the first time-step.

            iif = 1
            nfast = 1
            next_indx1 = 3 - indx1
            if not PREDICTOR_2D_STEP:
                PREDICTOR_2D_STEP = True
                if FIRST_2D_STEP:
                    kstp = indx1
                else:
                    kstp = 3 - indx1

                knew = 3
                krhs = indx1

            # Predictor step - Advance barotropic equations using 2D time-step
            # ==============   predictor scheme.
            call step2d_LF_AM3()


            # Set time indices for corrector step.
            if PREDICTOR_2D_STEP:
                PREDICTOR_2D_STEP = False
                knew = next_indx1
                kstp = 3 - knew
                krhs = 3
                if iif < nfast + 1:
                    indx1 = next_indx1

            # Corrector step - Apply 2D time-step corrector scheme.  Notice that
            # ==============   there is not need for a corrector step during the
            # auxiliary (nfast+1) time-step.

            if iif < nfast + 1:
                step2d()





