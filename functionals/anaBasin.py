# BASIN analytical case

import anaBASE

def getBasicGridParameters():
    # Gets the basic geometric parameters.

    Xsize = 3600.0E+03
    Esize = 2800.0E+03
    depth = 5000.0
    f0 = 1.0E-04
    beta = 2.0E-11

    return Xsize, Esize, depth, f0, beta


def computeAngle():
    return anaBASE.computeAngle()

def computeCoriolis():
    return anaBASE.computeCoriolis()


def computeCoordinates(Xsize, Esize, GRID):
    # Compute the (ξ,η) coordinates at PSI- and RHO-points.
    # Set grid spacing (m).
    # -----------------------------------------------------------------------

    # Determine I- and J-ranges for computing grid data.  These ranges
    # are special in periodic boundary conditons since periodicity cannot
    # be imposed in the grid coordinates.

    # Uses the basic method.
    return anaBASE.computeCoordinates(Xsize, Esize, GRID)


def computeStatistics():
    return anaBASE.computeStatistics()


def computeCoordinateTransform(GRID):
    return anaBASE.computeCoordinateTransform(GRID)


def computeBathymetry(depth, GRID):
    return anaBASE.computeBathymetry(depth, GRID)

def sources(GRID):