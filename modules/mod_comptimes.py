from misc import *


# Data about computational times.
# -----------------------------------------------------------------------
class CompTimes:
    def __init__(self, input):

        self.ntimes  = input.getVal('NTIMES',  dtype = int)            # Total number timesteps in current run. In 3D configurations, "ntimes"
                                                                  # is the total of baroclinic timesteps. In 2D configuration, "ntimes"
                                                                  # is the total of barotropic timesteps.
        self.dt      = input.getVal('DT',      dtype = float)          # Baroclinic timestep (s)
        self.ndtfast = input.getVal('NDTFAST', dtype = int)            # Number of barotropic timesteps between each baroclinic timestep.
        self.dtfast  = self.dt/self.ndtfast                                      # Barotropic timestep (s)

        self.LastRec = (input.getVal('NRREC', dtype = int) < 0)

        self.PerfectRST        = False
        self.PREDICTOR_2D_STEP = False


        self.indx1 = 0       # 2D timestep rolling counter.
        self.iic   = 0       # Timestep counter for 3D primitive equations
        self.iif   = 0       # Timestep counter for 2D primitive equations

        self.nfast   = 0     # Number of barotropic timesteps needed to compute time-averaged barotropic variables centered at time level n+1.

        self.tdays = 0.0     # (days)    Model time clock.
        self.time  = 0.0     # (seconds) Model time clock

        self.dtfast = 0.0

        self.TimeEnd  = 0    # (s)
        self.IMPtime  = 0    # (s) Impulse forcing time to process.
        self.INItime  = 0    # (s) Nonlinear model initial conditions time.
        self.INItimeS = 0    # (s)

        self.dstart        = input.getVal('DSTART',     dtype=float)    # (days)    Time stamp assigned to model initialization (usually a Calendar day, like modified Julian Day).
        self.tide_start    = input.getVal('TIDE_START', dtype=float)    # (days)    Reference time for tidal forcing

        self.F_code = ''     # Final time string for simulation.
        self.I_code = ''     # Initial time string for simulation.

        self.time_code = ''  # Model time clock (string, YYYY-MM-DD hh:mm:ss.ss).



        self.knew = 1        # Barotropic (fast) time-step index corresponding to the newest values for 2D primitive equation variables.
        self.krhs = 1        # Barotropic (fast) time-step index used to compute the right-hand-terms of 2D primitive equation variables.
        self.kstp = 1        # Barotropic (fast) time-step index to which the current changes are added to compute new 2D primitive equation variables.

        # Time-step counter for current execution time-window.
        self. step_counter = 0.0


        # First, starting, and ending timestepping parameters
        self.ntfirst = 0      # Forward-Euler step
        self.ntstart = 0      # Start step
        self.ntend   = 0      # End step


    def isFirst2DStep(self):
        # This works only for 2D XXX
        return self.iic == self.ntfirst


    def get2DTimes(self):
        # This function is used in the 2D functions.

        if self.isFirst2DStep():
            know = self.krhs
            dt2d = self.dtfast

        elif self.PREDICTOR_2D_STEP:
            know = self.krhs
            dt2d = 2.0*self.dtfast

        else:
            know = self.kstp
            dt2d = self.dtfast


        return know, dt2d