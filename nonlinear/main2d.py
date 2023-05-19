import mod_param
import mod_scalars
from ntimestep import ntimestep
from step2d import step2d


def main2d(RunInterval, compTimes,  GRID, OCEAN):
    """ This subroutine is the main driver for nonlinear ROMS/TOMS when configurated as shallow water (barotropic) ocean
    model only. It advances forward  the vertically integrated primitive equations for all nested grids,  if any,  by the
    specified  time interval (seconds), RunInterval.
    """


    # Time-step vertically integrated equations.
    # =======================================================================

    Time_Step = True


    # Initialize all time levels and compute other initial fields.
    # -----------------------------------------------------------------------
    if compTimes.isInitialTime():

        # Initialize free-surface.
        ini_zeta()

        # Initialize other state variables.
        ini_fields ()


        # Set vertical boundary conditions. Process tidal forcing.
        set_vbc()


    # Main loop
    while Time_Step:

        # Determine number of time steps to compute in each nested grid layer
        # based on the specified time interval (seconds), RunInterval. Non
        # nesting applications have NestLayers=1. Notice that RunInterval is
        # set in the calling driver. Its value may span the full period of the
        # simulation, a multi-model coupling interval (RunInterval > ifac*dt),
        # or just a single step (RunInterval=0).
        Nsteps, Rsteps = ntimesteps(RunInterval, Nsteps, Rsteps)


        # Time-step governing equations for Nsteps.
        for istep in range(Nsteps):


            # Set time indices and time clock. Get's the index of the first dimension of variables like OCEAN.zeta.
            compTimes.nextTimeStep(dt)
            compTimes.updateIndices(True)



            # Read in required data, if any, from input NetCDF files.
            # ----------------------------------------------------------------------
            get_data()

            # If applicable, process input data: time interpolate between data snapshots.
            # -----------------------------------------------------------------------
            set_data()



            # If appropriate, write out fields into output NetCDF files.
            output()


            if compTime.isFinalTimeStep():
                exitProgram()


            # Solve the vertically integrated primitive equations for the free-surface and momentum components.
            # -----------------------------------------------------------------------

            nfast = 1    # TODO: What is this?????? Can line 88 ever be true??

            # Predictor step - Advance barotropic equations using 2D time-step
            # ==============   predictor scheme.
            step2d(True, compTimes, GRID, OCEAN)


            # Computes the indices (of the first dimensions of variables like OCEAN.zeta) for the corrector step.
            compTimes.updateIndices(False)

            if compTimes.iif < nfast + 1:

                # Corrector step - Apply 2D time-step corrector scheme.  Notice that there is not need for a corrector step
                # during the auxiliary (nfast+1) time-step.
                step2d(False)





