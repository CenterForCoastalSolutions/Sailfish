# !================================================== Hernan G. Arango ===
# !  Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
# !    Licensed under a MIT/X style license                              !
# !    See License_ROMS.txt                                              !
# !=======================================================================
# !                                                                      !
# !  These routines apply close, gradient or periodic boundary           !
# !  conditions to generic 2D fields.                                    !
# !                                                                      !
# !  On Input:                                                           !
# !                                                                      !
# !     ng                Nested grid number.                            !
                                  !
# !                                                                      !
# !  On Output:                                                          !
# !                                                                      !
# !     var                 Processed 2D field.                            !
# !                                                                      !
# !  Routines:                                                           !
# !                                                                      !
# !     bc_r:-            Boundary conditions for field at RHO-points    !
# !     bc_u2d_tile       Boundary conditions for field at U-points      !
# !     bc_v2d_tile       Boundary conditions for field at V-points      !
# !                                                                      !
# !=======================================================================

import cupy as cp


def bc_r2d_tile(ng, vars):
    """ BC for rho-type cells: Boundary conditions are "imposed" by setting values to ghost  nodes"""

    for var in vars:
        # Sets zero-gradient boundary conditions.
        var.flat[idxDstZGradBCR[ng]] = var.flat[idxSrcZGradBCR[ng]]

        # Sets periodic boundary conditions.
        var.flat[idxDstPeriodicBCR[ng]] = var.flat[idxSrcPeriodicBCR[ng]]



def bc_u2d_tile(vars, gamma2):
    """ BC for U-type cells:  Boundary conditions are "imposed" by setting values to ghost  nodes"""

    for var in vars:
        # Gets a linear view of the array.
        var = var.ravel()

        # Sets zero-gradient boundary conditions.
        var[idxDstZGradBCU] = var[idxSrcZGradBCR]

        # Sets free-slip boundary conditions.
        var[idxDstZGradBCU] = gamma2*var[idxSrcZGradBCU]

        # Sets no-slip boundary conditions.
        var[idxDstClosedBCU] = 0.0

        # Sets periodic boundary conditions.
        var[idxDstPeriodicBCU] = var[idxSrcPeriodicBCU]


def bc_v2d_tile(vars, gamma2):
    """ BC for U-type cells:  Boundary conditions are "imposed" by setting values to ghost  nodes"""

    for var in vars:
        # Gets a linear view of the array.
        var = var.ravel()

        # Sets zero-gradient boundary conditions.
        var[idxDstZGradBCV] = var[idxSrcZGradBCV]

        # Sets free-slip boundary conditions.
        var[idxDstZGradBCV] = gamma2*var[idxSrcZGradBCV]

        # Sets no-slip boundary conditions.
        var[idxDstClosedBCV] = 0.0

        # Sets periodic boundary conditions.
        var[idxDstPeriodicBCV] = var[idxSrcPeriodicBCV]

