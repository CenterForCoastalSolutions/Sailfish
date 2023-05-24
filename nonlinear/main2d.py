import mod_param
import mod_scalars
from misc import *
# from get_data import get_data
# from set_data import set_data
from step2d import step2dPredictor, step2dCorrector


def main2d(RunInterval, compTimes,  GRID, OCEAN):
    """ This subroutine is the main driver for nonlinear ROMS/TOMS when configurated as shallow water (barotropic) ocean
    model only. It advances forward  the vertically integrated primitive equations for all nested grids,  if any,  by the
    specified  time interval (seconds), RunInterval.
    """


    # TODO: Recover this
    # # Initialize all time levels and compute other initial fields.
    # # -----------------------------------------------------------------------
    # if compTimes.isInitialTime():   # Not sure if this if is necessary.
    #
    #     # Initialize free-surface.
    #     ini_zeta()
    #
    #     # Initialize other state variables.
    #     ini_fields ()
    #
    #
    #     # Set vertical boundary conditions. Process tidal forcing.
    #     set_vbc()


    # Main loop
    while compTimes.keepRunning:


        # Set time indices and time clock. Get's the index of the first dimension of variables like OCEAN.zeta.
        compTimes.nextTimeStep(dt)
        OCEAN.cycleTimes()
        # compTimes.updateIndices(True)



        # Read in required data, if any, from input NetCDF files.
        # ----------------------------------------------------------------------
        # get_data()

        # If applicable, process input data: time interpolate between data snapshots.
        # -----------------------------------------------------------------------
        # set_data()



        # If appropriate, write out fields into output NetCDF files.
        # output()


        if compTimes.isFinalTimeStep():
            exitProgram()


        # Solve the vertically integrated primitive equations for the free-surface and momentum components.
        # -----------------------------------------------------------------------

        nfast = 1    # TODO: What is this?????? Can line 88 ever be true??

        # Predictor step - Advance barotropic equations using 2D time-step
        # ==============   predictor scheme.
        step2dPredictor(compTimes, GRID, OCEAN)


        # # Computes the indices (of the first dimensions of variables like OCEAN.zeta) for the corrector step.
        # compTimes.updateIndices(False)


        # if compTimes.iif < nfast + 1: TODO: Keep this extra step???

        # Corrector step - Apply 2D time-step corrector scheme.  Notice that there is not need for a corrector step
        # during the auxiliary (nfast+1) time-step.
        step2dCorrector(compTimes, GRID, OCEAN, BOUNDARY)





