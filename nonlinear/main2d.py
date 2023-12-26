from misc import *
# from get_data import get_data
# from set_data import set_data
from step2d import step2dPredictor, step2dCorrector

import matplotlib.pyplot as plt


def main2d(compTimes, GRID, OCEAN, BOUNDARY):
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


    import time
    t1 = time.time()
    print('TTTT', t1)

    # Main loop
    # ----------
    while not compTimes.isFinalTimeStep():


        # Gets next time step and cycles variables.
        compTimes.next2DTimeStep()
        OCEAN.cycleTimes2D()   # Fast time cycle (barotropic)


        # Read in required data, if any, from input NetCDF files.
        # ----------------------------------------------------------------------
        # get_data()

        # If applicable, process input data: time interpolate between data snapshots.
        # -----------------------------------------------------------------------
        # set_data()



        # If appropriate, write out fields into output NetCDF files.
        # output()



        # Solve the vertically integrated primitive equations for the free-surface and momentum components.
        # -----------------------------------------------------------------------

        nfast = 1    # TODO: What is this??????

        # Predictor step - Advance barotropic equations using 2D time-step
        # ==============   predictor scheme.
        step2dPredictor(compTimes, GRID, OCEAN, BOUNDARY)

        if compTimes.iif % 50000 == 0:
            plt.clf()
            plt.imshow(OCEAN.zeta_t2.reshape(GRID.M+1, GRID.M+1).get())
            plt.colorbar()
            plt.pause(1)
            print('.... %.2f s   ' % (compTimes.time2D))



        # if compTimes.iif < nfast + 1: TODO: Keep this extra step???

        # Corrector step - Apply 2D time-step corrector scheme.  XXXX-> Notice that there is not need for a corrector step
        # during the auxiliary (nfast+1) time-step.
        step2dCorrector(compTimes, GRID, OCEAN, BOUNDARY)


    t2 = time.time()
    print('Time ->', t2 - t1)
    exitProgram()




