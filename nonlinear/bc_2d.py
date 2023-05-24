
# These routines apply close, gradient or periodic boundary
# conditions to generic 2D fields.
#
# On Input:
#
#    var                 Processed 2D field.
#
# Routines:
#
#    bc_r:         Boundary conditions for field at RHO-points
#    bc_u2d:       Boundary conditions for field at U-points
#    bc_v2d:       Boundary conditions for field at V-points

import cupy as cp

import mod_boundary


def bc_r2d(vars, BOUNDARY):
    """ BC for rho-type cells: Boundary conditions are "imposed" by setting values to BC nodes"""

    idxZGradBC = ((BOUNDARY.zetaBC.LBC & mod_boundary.bcGradient) != 0)
    idxZGrad1 = BOUNDARY.zetaBC.bcIdx1[idxZGradBC]
    idxZGrad2 = BOUNDARY.zetaBC.bcIdx2[idxZGradBC]

    for var in vars:
        # Gets a flat view of the array.
        var = var.ravel()

        # Sets zero-gradient boundary conditions.
        var[idxZGrad1] = var[idxZGrad2]



def bc_u2d(vars, gamma2):
    """ BC for U-type cells:  Boundary conditions are "imposed" by setting values to ghost  nodes"""

    for var in vars:
        # Gets a flat view of the array.
        var = var.ravel()

        # Sets zero-gradient boundary conditions.
        var[idxDstZGradBCU] = var[idxSrcZGradBCR]

        # Sets free-slip boundary conditions.
        var[idxDstZGradBCU] = gamma2*var[idxSrcZGradBCU]

        # Sets no-slip boundary conditions.
        var[idxDstClosedBCU] = 0.0


def bc_v2d(vars, gamma2):
    """ BC for U-type cells:  Boundary conditions are "imposed" by setting values to ghost  nodes"""

    for var in vars:
        # Gets a flat view of the array.
        var = var.ravel()

        # Sets zero-gradient boundary conditions.
        var[idxDstZGradBCV] = var[idxSrcZGradBCV]

        # Sets free-slip boundary conditions.
        var[idxDstZGradBCV] = gamma2*var[idxSrcZGradBCV]

        # Sets no-slip boundary conditions.
        var[idxDstClosedBCV] = 0.0

