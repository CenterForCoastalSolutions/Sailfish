from misc import *
import cupy as cp


# Data about computational times.
# -----------------------------------------------------------------------
class CompTimes:
    def __init__(self, input):

        self.ntimes  = input.getVal('NTIMES',  dtype = int)            # Total number timesteps in current run. In 3D configurations, "ntimes"
                                                                  # is the total of baroclinic timesteps. In 2D configuration, "ntimes"
                                                                  # is the total of barotropic timesteps.
        self.dt      = input.getVal('DT',      dtype = float)          # Baroclinic timestep (s)
        self.ndtfast = input.getVal('NDTFAST', dtype = int)            # Number of barotropic timesteps between each baroclinic timestep.
        self.dtfast  = self.dt/self.ndtfast
        self.nfast   = 0

        self.LastRec = (input.getVal('NRREC', dtype = int) < 0)

        self.PerfectRST        = False
        self.is2DPredictorStep = False
        self.is2DCorrectorStep = True


        self.indx1 = 0       # 2D timestep rolling counter.
        self.iic   = 0       # Timestep counter for 3D primitive equations
        self.iif   = 0       # Timestep counter for 2D primitive equations

        self.nfast   = 0     # Number of barotropic timesteps needed to compute time-averaged barotropic variables centered at time level n+1.

        self.tdays  = 0.0     # (days)    Model time clock.
        self.time   = 0.0     # (seconds) Model time clock
        self.time2D = 0.0    # (seconds) Model 2D time clock


        self.TimeEnd  = 0    # (s)
        self.IMPtime  = 0    # (s) Impulse forcing time to process.
        self.INItime  = 0    # (s) Nonlinear model initial conditions time.
        self.INItimeS = 0    # (s)

        self.dstart        = input.getVal('DSTART',     dtype=float)    # (days)    Time stamp assigned to model initialization (usually a Calendar day, like modified Julian Day).
        self.tide_start    = input.getVal('TIDE_START', dtype=float)    # (days)    Reference time for tidal forcing

        self.F_code = ''     # Final time string for simulation.
        self.I_code = ''     # Initial time string for simulation.

        self.time_code = ''  # Model time clock (string, YYYY-MM-DD hh:mm:ss.ss).


        # Time-step counter for current execution time-window.
        self.step_counter = 0.0


        # First, starting, and ending timestepping parameters
        self.ntfirst = 0      # Forward-Euler step
        self.ntstart = 0      # Start step
        self.ntend   = 500      # End step  TODO: READ this from the file


        # Weights for time filtering.
        self.weight1 = cp.zeros((257), dtype=cp.float64)
        self.weight2 = cp.zeros((257), dtype=cp.float64)


    def isInitialTimeStep(self):
        return self.iic == self.ntstart


    def isFinalTimeStep(self):
        # return False     # Remove
        return self.iic >= self.ntend + 1


    def isFirst2DStep(self):
        # Not sure what is the difference between initial and First (or ntstart and ntfirst)
        return self.iif == 0


    def first2DTimeStep(self):
        self.iif = 0

    def next2DTimeStep(self):
        """Advances one time step"""

        self.iif += 1
        self.time2D += self.dtfast


    def nextTimeStep(self):
        """Advances one time step"""

        self.iic += 1
        self.time += self.dt
        self.time2D = self.time
        # self.tdays = self.time * sec2day    TODO : Recover this




    def getDeltaTime(self):
        # This function is used in the 2D functions. It takes into account that all predictor except for time = 0 do leapfrog and hence, the

        if self.is2DPredictorStep:
            return 2.0*self.dtfast

        else:
            return self.dtfast


