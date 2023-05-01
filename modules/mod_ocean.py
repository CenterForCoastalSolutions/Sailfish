import cupy as cp

# !  2D Primitive Variables.                                             !
# !                                                                      !
# !  rubar        Right-hand-side of 2D U-momentum equation (m4/s2).     !
# !  rvbar        Right-hand-side of 2D V-momentum equation (m4/s2).     !
# !  rzeta        Right-hand-side of free surface equation (m3/s).       !
# !  ubar         Vertically integrated U-momentum component (m/s).      !
# !  vbar         Vertically integrated V-momentum component (m/s).      !
# !  zeta         Free surface (m).                                      !
                                                                   !


# Nonlinear model state.
class T_OCEAN:
    def __init__(self, Ni, Nj):
        # Nonlinear model state.
        rubar = cp.zeros(Ni, Nj, 2)
        rvbar = cp.zeros(Ni, Nj, 2)
        rzeta = cp.zeros(Ni, Nj, 2)
        ubar  = cp.zeros(Ni, Nj, 3)
        vbar  = cp.zeros(Ni, Nj, 3)
        zeta  = cp.zeros(Ni, Nj, 3)



def allocate_ocean(LBi, UBi, LBj, UBj):
# This routine allocates all variables in the module for all nested grids.                                                              !


# Allocate and initialize module variables.
# -----------------------------------------------------------------------

    OCEAN = T_OCEAN(UBi - LBi,UBj - LBj)

    # Set horizontal array size.
    size2d=REAL((UBi-LBi+1)*(UBj-LBj+1),r8)



def initialize_ocean():
# This routine initialize all variables in the module using first     !
# touch distribution policy.


    IniVal = 0.0


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


    # Initialize module variables.
    #-----------------------------------------------------------------------
    self.rubar = createVar(shape2D, 2, IniVal)
    self.rvbar = createVar(shape2D, 2, IniVal)
    self.rzeta = createVar(shape2D, 2, IniVal)


    self.rubar = createVar(shape2D, 3, IniVal)
    self.rvbar = createVar(shape2D, 3, IniVal)
    self.rzeta = createVar(shape2D, 3, IniVal)



