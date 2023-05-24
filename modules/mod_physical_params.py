# This routine reads and reports physical model input parameters.
# =======================================================================

# import mod_param
# import mod_iounits
# import mod_ncparam
# import mod_netcdf
# import mod_scalars
# import mod_strings
import numpy as cp

# import dateclock

from misc import *


class PhysicalParams:

    def __init__(self, input, GRID):


        self.itracer  = 0             # LBC tracer counter


        self.title    = input.getVal('TITLE')
        self.MyAppCPP = input.getVal('MyAppCPP')
        self.varname  = input.getVal('VARNAME')

        self.rho0    = input.getVal('RHO0',      dtype = float)
        self.bvf_bak = input.getVal('BVF_BAK',   dtype = float)


        time_refStr = input.getVal('TIME_REF')
        time_ref    = input.getVal('TIME_REF', dtype = float)
        msgInfo('Recover this.')
        # dateclock.ref_clock(time_ref, time_refStr) TODO: Recover this.

        # MIXING of MOMENTUM and TRACERS
        # nl_tnu2    = input.getVal('TNU2', count = NAT + NPT, dtype = float)
        # nl_tnu4    = input.getVal('TNU4', count = NAT + NPT, dtype = float)

        self.nl_visc2   = input.getVal('VISC2',     dtype = float)
        self.nl_visc4   = input.getVal('VISC4',     dtype = float)
        self.LuvSponge  = input.getVal('LuvSponge', dtype = bool)

        # Akt_bak    = input.getVal('AKT_BAK',   count = NAT + NPT, dtype = float)
        # Akt_limit  = input.getVal('AKT_LIMIT', count = NAT + NPT, dtype = float)

        self.Akv_bak    = input.getVal('AKV_BAK',    dtype = float)
        self.Akv_limit  = input.getVal('AKV_LIMIT',  dtype = float)
        self.ad_Akv_fac = input.getVal('ad_AKV_fac', dtype = float)
        self.Akk_bak    = input.getVal('AKK_BAK',    dtype = float)
        self.Akp_bak    = input.getVal('AKP_BAK',    dtype = float)
        self.tkenu2     = input.getVal('TKENU2',     dtype = float)
        self.tkenu4     = input.getVal('TKENU4',     dtype = float)


        # STATE EQUATION
        self.R0         = input.getVal('R0',    dtype = float, minVal = 0.0)
        self.T0         = input.getVal('T0',    dtype = float)
        self.S0         = input.getVal('S0',    dtype = float)
        self.Tcoef      = cp.abs(input.getVal('TCOEF', dtype = float))
        self.Scoef      = cp.abs(input.getVal('SCOEF', dtype = float))

        if (self.R0 < 100):
            self.R0 += 1000.0  # Apparently, one can write 24 instead of 1024



        self.gamma2    = input.getVal('GAMMA2', dtype = float)

        self.LuvSrc    = input.getVal('LuvSrc', dtype = bool)
        self.LwSrc     = input.getVal('LwSrc',  dtype = bool)


        # LsshCLM     = input.getVal('LsshCLM',     dtype=float)
        # Lm2CLM      = input.getVal('Lm2CLM',      dtype=float)
        # LnudgeM2CLM = input.getVal('LnudgeM2CLM', dtype=float)
        # nl_tnu4     = input.getVal('TNU4',        dtype=float)
        # nl_tnu4     = input.getVal('TNU4',        dtype=float)


        #  Check if both point sources methodologies are activated. Only one method is allowed for a particular grid.  Otherwise, the point
        #  source will be applied twice.
        if self.LuvSrc and self.LwSrc:
            msgError('Both LuvSrc and LuvSrc are selected. Only one method is legal', 4)






