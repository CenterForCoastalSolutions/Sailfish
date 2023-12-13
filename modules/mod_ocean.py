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
        N = GRID.N
        M = GRID.M
        L = GRID.L


        # 2D momentum
        # -----------

        # RHS of the differential equations (like dzeta/dt = rzeta).
        self.rubar = cp.zeros((3, M + 1, L + 1), dtype = cp.float64)
        self.rvbar = cp.zeros((3, M + 1, L + 1), dtype = cp.float64)
        self.rzeta = cp.zeros((3, M + 1, L + 1), dtype = cp.float64)

        # Nonlinear model state.
        self.ubar  = cp.zeros((3, M + 1, L + 1), dtype = cp.float64)
        self.vbar  = cp.zeros((3, M + 1, L + 1), dtype = cp.float64)
        self.zeta  = cp.zeros((3, M + 1, L + 1), dtype = cp.float64)


        # These aliases are used in the predictor/corrector steps. t2, t1 and t0 always mean t+2Δt, t+Δt and t,
        # but their location in the arrays changes (as given by cycleTimes2D) to avoid copying lots of data.
        self.zeta_t2 = self.zeta[2, :, :].ravel()
        self.zeta_t1 = self.zeta[1, :, :].ravel()
        self.zeta_t0 = self.zeta[0, :, :].ravel()

        self.ubar_t2 = self.ubar[2, :, :].ravel()
        self.ubar_t1 = self.ubar[1, :, :].ravel()
        self.ubar_t0 = self.ubar[0, :, :].ravel()

        self.vbar_t2 = self.vbar[2, :, :].ravel()
        self.vbar_t1 = self.vbar[1, :, :].ravel()
        self.vbar_t0 = self.vbar[0, :, :].ravel()

        self.rzeta_t2 = self.rzeta[2, :, :].ravel()
        self.rzeta_t1 = self.rzeta[1, :, :].ravel()
        self.rzeta_t0 = self.rzeta[0, :, :].ravel()

        self.rubar_t2 = self.rubar[2, :, :].ravel()
        self.rubar_t1 = self.rubar[1, :, :].ravel()
        self.rubar_t0 = self.rubar[0, :, :].ravel()

        self.rvbar_t2 = self.rvbar[2, :, :].ravel()
        self.rvbar_t1 = self.rvbar[1, :, :].ravel()
        self.rvbar_t0 = self.rvbar[0, :, :].ravel()



        # 3D momentum
        # -----------
        self.ru = cp.zeros((3, N + 1, M + 1, L + 1), dtype = cp.float64)
        self.rv = cp.zeros((3, N + 1, M + 1, L + 1), dtype = cp.float64)

        # Nonlinear model state.
        self.u  = cp.zeros((3, N + 1, M + 1, L + 1), dtype = cp.float64)
        self.v  = cp.zeros((3, N + 1, M + 1, L + 1), dtype = cp.float64)
        self.W  = cp.zeros((3, N + 1, M + 1, L + 1), dtype = cp.float64)


        self.Huon  = cp.zeros((N + 1, M + 1, L + 1), dtype = cp.float64)    # Total U-momentum flux term, Hz*u/pn.
        self.Hvom  = cp.zeros((N + 1, M + 1, L + 1), dtype = cp.float64)    # Total V-momentum flux term, Hz*v/pm.

        self.DUon  = cp.zeros((M + 1, L + 1), dtype = cp.float64).ravel()
        self.DVom  = cp.zeros((M + 1, L + 1), dtype = cp.float64).ravel()


        # Coupling averages.
        # TODO: This needs to be filled somewhere.
        self.Zt_avg1 = cp.zeros((M + 1, L + 1), dtype = cp.float64).ravel()
        self.U_avg1  = cp.zeros((M + 1, L + 1), dtype = cp.float64).ravel()
        self.V_avg1  = cp.zeros((M + 1, L + 1), dtype = cp.float64).ravel()
        self.DU_avg1  = cp.zeros((M + 1, L + 1), dtype = cp.float64).ravel()
        self.DV_avg1  = cp.zeros((M + 1, L + 1), dtype = cp.float64).ravel()

        # Turbulence.
        # TODO: This needs to be filled somewhere.
        self.AKv  = 0.00001 + cp.zeros((N + 1, M + 1, L + 1), dtype = cp.float64)
        self.AK   = 0.00001



        # These aliases are used in the predictor/corrector steps. t2, t1 and t0 always mean t+2LΔt, t+LΔt and t,
        # but their location in the arrays changes (as given by cycleTimes) to avoid copying lots of data.
        self.u_t2 = self.u[2, :, :, :].ravel()
        self.u_t1 = self.u[1, :, :, :].ravel()
        self.u_t0 = self.u[0, :, :, :].ravel()

        self.v_t2 = self.v[2, :, :, :].ravel()
        self.v_t1 = self.v[1, :, :, :].ravel()
        self.v_t0 = self.v[0, :, :, :].ravel()



        self.ru_t2 = self.ru[2, :, :, :].ravel()
        self.ru_t1 = self.ru[1, :, :, :].ravel()
        self.ru_t0 = self.ru[0, :, :, :].ravel()

        self.rv_t2 = self.rv[2, :, :, :].ravel()
        self.rv_t1 = self.rv[1, :, :, :].ravel()
        self.rv_t0 = self.rv[0, :, :, :].ravel()







    def cycleTimes2D(self):
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
        self.rzeta_t1 = self.rzeta_t2
        self.rzeta_t2 = temp

        temp = self.rubar_t0
        self.rubar_t0 = self.rubar_t1
        self.rubar_t1 = self.rubar_t2
        self.rubar_t2 = temp

        temp = self.rvbar_t0
        self.rvbar_t0 = self.rvbar_t1
        self.rvbar_t1 = self.rvbar_t2
        self.rvbar_t2 = temp


    def cycleTimes3D(self):
        # Cycle the variables that contain the previous times. In reality this is a shift, where variables are moved
        # one barotropic time-step in such a way that t2 always point to time = t + 2LΔt and so on.
        # The reason for cycling is to reuse the same memory. the t2 variables contain "garbage" at this point. They
        # must be filled with newly computed values.

        temp = self.u_t0
        self.u_t0 = self.u_t1
        self.u_t1 = self.u_t2
        self.u_t2 = temp

        temp = self.v_t0
        self.v_t0 = self.v_t1
        self.v_t1 = self.v_t2
        self.v_t2 = temp

        temp = self.ru_t0
        self.ru_t0 = self.ru_t1
        self.ru_t1 = self.ru_t2
        self.ru_t2 = temp

        temp = self.rv_t0
        self.rv_t0 = self.rv_t1
        self.rv_t1 = self.rv_t2
        self.rv_t2 = temp


    def getVars(self):
        return (self.zeta_t2, self.zeta_t1, self.zeta_t0,
                self.ubar_t2, self.ubar_t1, self.ubar_t0,
                self.vbar_t2, self.vbar_t1, self.vbar_t0,
                self.rzeta_t1, self.rzeta_t0,
                self.rubar_t1, self.rubar_t0,
                self.rvbar_t1, self.rvbar_t0)








