

#import modules
import mod_iounits
import cupy as cp

def ntimesteps(model, RunInterval):
    """ This routine set the number of time-steps to compute. In nesting
    applications,  the number of time-steps depends on nesting layer
    configuration.
    """

    # On Input:
    #model              calling model identifier (integer)
    #RunInterval        Time interval (seconds) to advance forward or backwards the primitive euqations (scalar)

    # On Output:
    #RunInterval

    WindowSteps=0
    Nsteps=0
    Rsteps=0



    # Here, extra=1, indicates that the RunInterval is the same as simulation interval.
    extra=1

    #determine number of steps in time-interval window
    WindowSteps = int((RunInterval + 0.5 * dt)/dt)

    #Advancing model forward: Nonlinear
    my_Nsteps = cp.max(my_Nsteps, WindowSteps + extra)
    #set the number of time-steps to execute. Choose minimum values
    Nsteps = my_Nsteps(1)
    Rsteps = WindowSteps(1) + extra
    Nsteps = cp.min(Nsteps, my_Nsteps)
    Rsteps = cp.min(Rsteps, WindowSteps + extra)

    return Nsteps, Rsteps