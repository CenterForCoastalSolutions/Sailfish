import cupy as cp

# !  2D Primitive Variables.                                             !
# !                                                                      !
# !  rubar        Right-hand-side of 2D U-momentum equation (m4/s2).     !
# !  rvbar        Right-hand-side of 2D V-momentum equation (m4/s2).     !
# !  rzeta        Right-hand-side of free surface equation (m3/s).       !
# !  ubar         Vertically integrated U-momentum component (m/s).      !
# !  vbar         Vertically integrated V-momentum component (m/s).      !
# !  zeta         Free surface (m).                                      !


# Nonlinear model state.
class T_OCEAN:
    def __init__(self, input, GRID):
        M = GRID.M
        L = GRID.L

        # RHS of the differential equations (like dzeta/dt = rzeta).
        self.rubar = cp.zeros((2, M + 1, L + 1), dtype = cp.float64)
        self.rvbar = cp.zeros((2, M + 1, L + 1), dtype = cp.float64)
        self.rzeta = cp.zeros((2, M + 1, L + 1), dtype = cp.float64)

        # Nonlinear model state.
        self.ubar  = cp.zeros((3, M + 1, L + 1), dtype = cp.float64)
        self.vbar  = cp.zeros((3, M + 1, L + 1), dtype = cp.float64)
        self.zeta  = cp.zeros((3, M + 1, L + 1), dtype = cp.float64)


        # These aliases are used in the predictor/corrector steps. t2, t1 and t0 always mean t+2Δt, t+Δt and t,
        # but their location in the arrays changes (as given by cycleTimes) to avoid copying lots of data.
        self.zeta_t2 = self.zeta[2, :, :].ravel()
        self.zeta_t1 = self.zeta[1, :, :].ravel()
        self.zeta_t0 = self.zeta[0, :, :].ravel()

        self.ubar_t2 = self.ubar[2, :, :].ravel()
        self.ubar_t1 = self.ubar[1, :, :].ravel()
        self.ubar_t0 = self.ubar[0, :, :].ravel()

        self.vbar_t2 = self.vbar[2, :, :].ravel()
        self.vbar_t1 = self.vbar[1, :, :].ravel()
        self.vbar_t0 = self.vbar[0, :, :].ravel()

        self.rzeta_t1 = self.rzeta[1, :, :].ravel()
        self.rzeta_t0 = self.rzeta[0, :, :].ravel()

        self.rubar_t1 = self.rubar[1, :, :].ravel()
        self.rubar_t0 = self.rubar[0, :, :].ravel()

        self.rvbar_t1 = self.rvbar[1, :, :].ravel()
        self.rvbar_t0 = self.rvbar[0, :, :].ravel()



    def cycleTimes(self):
        # Cycle the variables that contain the previous times. In reality this is a shift, where variables are moved
        # one barotropic time-step in such a way that t2 always point to time = t + 2Δt and so on.
        # The reason for cycling is to reuse the same memory. the t2 variables contain "garbage" at this point. They
        # must be filled with newly computed values.

        temp = self.zeta_t0
        self.zeta_t0 = self.zeta_t1
        self.zeta_t1 = self.zeta_t2
        self.zeta_t2 = temp

        temp = self.ubar_t0
        self.ubar_t0 = self.ubar_t1
        self.ubar_t1 = self.ubar_t2
        self.ubar_t2 = temp

        temp = self.vbar_t0
        self.vbar_t0 = self.vbar_t1
        self.vbar_t1 = self.vbar_t2
        self.vbar_t2 = temp

        temp = self.rzeta_t0
        self.rzeta_t0 = self.rzeta_t1
        self.rzeta_t1 = temp

        temp = self.rubar_t0
        self.rubar_t0 = self.rubar_t1
        self.rubar_t1 = temp

        temp = self.rvbar_t0
        self.rvbar_t0 = self.rvbar_t1
        self.rvbar_t1 = temp


    def getVars(self):
        return (self.zeta_t2, self.zeta_t1, self.zeta_t0,
                self.ubar_t2, self.ubar_t1, self.ubar_t0,
                self.vbar_t2, self.vbar_t1, self.vbar_t0,
                self.rzeta_t1, self.rzeta_t0,
                self.rubar_t1, self.rubar_t0,
                self.rvbar_t1, self.rvbar_t0)





# Set array initialization range.

    # IF (DOMAIN(ng)%Western_Edge(tile)) THEN
    #     Imin=BOUNDS(ng)%LBi(tile)
    # ELSE
    #     Imin=Istr
    # END IF
    # IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
    #     Imax=BOUNDS(ng)%UBi(tile)
    # ELSE
    #     Imax=Iend
    # END IF
    # IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
    #     Jmin=BOUNDS(ng)%LBj(tile)
    # ELSE
    #     Jmin=Jstr
    # END IF
    # IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
    #     Jmax=BOUNDS(ng)%UBj(tile)
    # ELSE
    #     Jmax=Jend
    # END IF






