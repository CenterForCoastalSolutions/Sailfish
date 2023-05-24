import cupy as cp

def computeCoordinates(Xsize, Esize, GRID):
    # Compute the (ξ,η) coordinates at PSI- and RHO-points.
    # Set grid spacing (m).
    # -----------------------------------------------------------------------

    # Determine I- and J-ranges for computing grid data.  These ranges
    # are special in periodic boundary conditons since periodicity cannot
    # be imposed in the grid coordinates.

    L  = GRID.L
    M  = GRID.M
    Lm = GRID.Lm
    Mm = GRID.Mm

    I, J = cp.meshgrid(cp.arange(0, M+1), cp.arange(0, L+1))
    GRID.I, GRID.J = (I, J)

    dx = Xsize/Lm
    dy = Esize/Mm
    GRID.dx = dx*cp.ones(I.shape)
    GRID.dy = dy*cp.ones(I.shape)


    GRID.xp = dx * (I - 1)
    GRID.yp = dy * (J - 1)

    GRID.xr = GRID.xp + 0.5*dx
    GRID.yr = GRID.yp + 0.5*dy

    GRID.xu = GRID.xp
    GRID.yu = GRID.yr

    GRID.xv = GRID.xr
    GRID.yv = GRID.yp

    return GRID


def computeAngle():
    return None

def computeCoriolis():
    return None


def computeBathymetry(depth, GRID):

    h = cp.zeros((GRID.M, GRID.L))

    GRID.h[:,:] = depth


def computeStatistics():
    pass


def computeCoordinateTransform(GRID):
    GRID.pm[:,:] = 1/GRID.dx
    GRID.pn[:,:] = 1/GRID.dy

    GRID.dndξ[:, :] = 0.0
    GRID.dmdη[:, :] = 0.0


    return GRID