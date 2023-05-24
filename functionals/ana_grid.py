import mod_param
import mod_param
import mod_grid
# import mod_ncparam
import cupy as cp
import importlib
import anaBASE

def ana_grid(analiticalCase, GRID):
    """
    This routine sets model grid using an analytical expressions.

    On Output:  stored in common blocks:

       f        Coriolis parameter (1/seconds) at RHO-points.
       h        Bathymetry (meters; positive) at RHO-points.
       hmin     Minimum depth of bathymetry (m).
       hmax     Maximum depth of bathymetry (m).
       pm       Coordinate transformation metric "m" (1/meters) associated with the differential distances in XI at RHO-points.
       pn       Coordinate transformation metric "n" (1/meters) associated with the differential distances in ETA at RHO-points.
       xp       XI-coordinates (m) at PSI-points.
       xr       XI-coordinates (m) at RHO-points.
       yp       ETA-coordinates (m) at PSI-points.
       yr       ETA-coordinates (m) at RHO-points.
    """

    # Imports the specified analytical function module
    modAnalitical = importlib.import_module('ana%s' % analiticalCase)


    # Set grid parameters:
    #    Xsize    Length (m) of domain box in the XI-direction.
    #    Esize    Length (m) of domain box in the ETA-direction.
    #    depth    Maximum depth of bathymetry (m).
    #    f0       Coriolis parameter, f-plane constant (1/s).
    #    beta     Coriolis parameter, beta-plane constant (1/s/m).
    Xsize, Esize, depth, f0, beta = modAnalitical.getBasicGridParameters()



    # Compute the (ξ,η) coordinates at U, V, PSI and RHO points.
    # Set grid spacing (m).
    # Determine I- and J-ranges for computing grid data.  These ranges
    # are special in periodic boundary conditons since periodicity cannot
    # be imposed in the grid coordinates.
    A = modAnalitical.computeCoordinates(Xsize, Esize, GRID)


    modAnalitical.computeStatistics()


    # Compute coordinate transformation metrics at RHO-points "pm" and
    # "pn" (1/m) associated with the differential distances in ξ and
    # η, respectively.
    # Compute d(1/n)/dξ and d(1/m)/dη at RHO-points.
    # pm, pn, dndξ, dmdη
    GRID = modAnalitical.computeCoordinateTransform(GRID)


    # Compute angle (in radians) between ξ-axis and true EAST at RHO-points.
    # -----------------------------------------------------------------------
    angler = modAnalitical.computeAngle()


    # Compute Coriolis parameter (1/s) at RHO-points.
    f = modAnalitical.computeCoriolis()


    # Compute bathymetry (meters, positive down) at RHO-points.
    h = modAnalitical.computeBathymetry(depth, GRID)

