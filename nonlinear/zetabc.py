from misc import *
import cupy as cp

def zetabc(zeta, compTimes, BOUNDARY):
    """This routine sets lateral boundary conditions for free-surface."""


    # Gradient/closed boundary condition.
    # "This boundary condition is extremely simple and consists of setting the gradient of a field to zero
    # at the edge. The outside value is set equal to the closest interior value"
    zeta[BOUNDARY.zetaClosedOrGradientBCIdx1] = zeta[BOUNDARY.zetaClosedOrGradientBCIdx2]


    # Clamped boundary condition.
    # "Very simple BC that consist in setting the boundary value to a known exterior value"
    # zetaKout[BOUNDARY.zetaClampedBCIdx1] = BOUNDARY.zeta[BOUNDARY.zetaClampedBCIdx1]
    msgInfo('Implement the real Clamped BC, here I am using a fake one!!!')
    omega = 0.01  # s^-1
    zeta[BOUNDARY.zetaClampedBCIdx1] = cp.sin(compTimes.time*omega)
















#     # Implicit upstream radiation condition.
#     # "In realistic domains, open boundary conditions can be extremely difficult to get right. There are
#     # situations in which incoming flow and outgoing flow happen along the same boundary or even at
#     # different depths at the same horizontal location. Orlanski [1976] proposed a radiation scheme in
#     # which a local normal phase velocity is computed and used to radiate things out (if it is indeed going
#     # out). This works well for a wave propagating normal to the boundary, but has problems when waves
#     # approach the boundary at an angle. Raymond and Kuo [1984] have modified the scheme to account
#     # for propagation in all three directions. In ROMS, only the two horizontal directions are accounted
#     # for (with the recommended RADIATION_2D option). (See documentation for more info and formulas)
#     elif lbc.radiation:
#
#         # radiationBCKernel()
#
#         grad(idxLateralBC1) = zeta(know, idxLateralBC1) - zeta(know, idxLateralBC2);
#         grad(idxLateralBC3) = zeta(know, idxLateralBC3) - zeta(know, idxLateralBC4);
#
#
#         dZdt = zeta(idxLateralBC1,know) - zeta(idxLateralBC1,kout)
#         dZdη = zeta(idxLateralBC1,kout) - zeta(idxLateralBC2,kout)
#
#         dZdξ = diffVtoU(idxLateralBC1, idxLateralBC2)
#
#
#
#         if dZdt*dZdη < 0.0:
#             dZdt=0.0
#
#         dZdξ = upwndDiff(idxLateralBC1, idxLateralBC2, dZdt)
#           # # Upwind?
#           # if ((dZdt*(grad(Istr,j) + grad(Istr,j+1))) > 0.0):
#           #   dZdξ =grad(Istr,j  )
#           # else:
#           #   dZdξ =grad(Istr,j+1)
#
#
#
#           dξZ2_dηZ2 = max(dZdξ*dZdξ + dZdη*dZdη, eps)
#           Cη = dZdt*dZdη
#           Cξ = 0.0
#           if RADIATION_2D:
#             Cξ = clamp(dZdt*dZdξ, -dξZ2_dηZ2, dξZ2_dηZ2)
#
#
# # #if defined CELERITY_WRITE && defined FORWARD_WRITE
# #           boundary.zeta_west_Cx(j) = Cη
# #           boundary.zeta_west_Cξ(j) = Cξ
# #           boundary.zeta_west_C2(j) = dξZ2_dηZ2
# # #endif
#           zeta[:, kout] = (dξZ2_dηZ2*zeta[:,know] +
#  &                          Cη*zeta[idx2,kout] -             &
#  &                          max(Cξ, 0.0)*grad(Istr-1,j  ) -     &
#  &                          min(Cξ, 0.0)*grad(Istr-1,j+1))/    &
#  &                          (dξZ2_dηZ2 + Cη)
#
#     if lbc.nudging:
#         # Nudging coefficients (1/s) for passive/active (outflow/inflow) open boundary conditions.
#         if dZdt * dZdx < 0.0:
#             tau = FSobc_in(ng, iwest)
#         else:
#             tau = FSobc_out(ng, iwest)
#         tau = tau * dt2d
#
#         zeta(idx,kout) = zeta(idx, kout) +
#                        tau*(boundary.zeta_west(j) -
#         &                             zeta(idx,know))
#
#
#     # Chapman boundary condition.
#     # "The condition for surface elevation to be used with either the Flather or Shchepetkin momentum
#     # boundary is that of [Chapman, 1985], assuming all outgoing signals leave at the shallow-water wave
#     # speed of sqrt(g*D)". (See documentation for more info.)
#     elif lbc.Chapman:
#         # Dcrit, WET_DRY
#         chapmanBCKernel(grid.h, zeta[:,know], grid.pm, idx, idx2, g, Dcrit, WET_DRY)
#












