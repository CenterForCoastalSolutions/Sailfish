import numpy as cp

def computeCoordinates(Xsize, Esize):
    # Compute the (ξ,η) coordinates at PSI- and RHO-points.
    # Set grid spacing (m).
    # -----------------------------------------------------------------------

    # Determine I- and J-ranges for computing grid data.  These ranges
    # are special in periodic boundary conditons since periodicity cannot
    # be imposed in the grid coordinates.

    I, J = cp.meshgrid(cp.arange(0, Mm + 1), cp.arange(0, Lm + 1))

    dx = Xsize/Lm
    dy = Esize/Mm


    xp = dx * (I - 1)
    yp = dy * (J - 1)

    xr = xp + 0.5*dx
    yr = yp + 0.5*dy

    xu = xp
    yu = yr

    xv = xr
    yv = yp

    return I, J, dx, dy, xp, yp, xr, yr, xu, yu, xv, yv


def computeBathymetry(depth):

    h[:.:] = depth