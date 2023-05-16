import numpy as cp

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

        # Nonlinear model state.
        self.rubar = cp.zeros((2, M + 1, L + 1), dtype = cp.float64)
        self.rvbar = cp.zeros((2, M + 1, L + 1), dtype = cp.float64)
        self.rzeta = cp.zeros((2, M + 1, L + 1), dtype = cp.float64)

        self.ubar  = cp.zeros((3, M + 1, L + 1), dtype = cp.float64)
        self.vbar  = cp.zeros((3, M + 1, L + 1), dtype = cp.float64)
        self.zeta  = cp.zeros((3, M + 1, L + 1), dtype = cp.float64)



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






