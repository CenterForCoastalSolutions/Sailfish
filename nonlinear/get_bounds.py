import mod_param
import mod_parallel
import mod_ncparam
import mod_scalars


# ROMS uses at staggered stencil:
#
#       -------v(i,j+1,k)-------               ------W(i,j,k)-------
#       |                      |               |                   |
#    u(i,j,k)   r(i,j,k)   u(i+1,j,k)          |     r(i,j,k)      |
#       |                      |               |                   |
#       --------v(i,j,k)--------               -----W(i,j,k-1)------
#
#           horizontal stencil                   vertical stencil
#                C-grid
#
#
# M   r..u..r..u..r..u..r..u..r..u..r..u..r..u..r..u..r..u..r..u..r
#     :                                                           :
#  M  v  p++v++p++v++p++v++p++v++p++v++p++v++p++v++p++v++p++v++p  v
#     :  +     |     |     |     |     |     |     |     |     +  :
# Mm  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r
#     :  +     |     |     |     |     |     |     |     |     +  :
#  Mm v  p--v--p--v--p--v--p--v--p--v--p--v--p--v--p--v--p--v--p  v
#     :  +     |     |     |     |     |     |     |     |     +  :
#     r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r
#     :  +     |     |     |     |     |     |     |     |     +  :
#     v  p--v--p--v--p--v--p--v--p--v--p--v--p--v--p--v--p--v--p  v
#     :  +     |     |     |     |     |     |     |     |     +  :
#     r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r
#     :  +     |     |     |     |     |     |     |     |     +  :
#     v  p--v--p--v--p--v--p--v--p--v--p--v--p--v--p--v--p--v--p  v
#     :  +     |     |     |     |     |     |     |     |     +  :
# 2   r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r
#     :  +     |     |     |     |     |     |     |     |     +  :
#  2  v  p--v--p--v--p--v--p--v--p--v--p--v--p--v--p--v--p--v--p  v
#     :  +     |     |     |     |     |     |     |     |     +  :
# 1   r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r
#     :  +     |     |     |     |     |     |     |     |     +  :
#  1  v  p++v++p++v++p++v++p++v++p++v++p++v++p++v++p++v++p++v++p  v
#     :                                                           :
# 0   r..u..r..u..r..u..r..u..r..u..r..u..r..u..r..u..r..u..r..u..r
#        1     2                                         Lm    L
#     0     1     2                                         Lm    L
#
#                          interior       Boundary Conditions
#                        computations     W     E     S     N
#
#   RH0-type variables:  [1:Lm, 1:Mm]   [0,:] [L,:] [:,0] [:,M]
#   PSI-type variables:  [2:Lm, 2:Mm]   [1,:] [L,:] [:,1] [:,M]
#     U-type variables:  [2:Lm, 1:Mm]   [1,:] [L,:] [:,0] [:,M]
#     V-type variables:  [1:Lm, 2:Mm]   [0,:] [L,:] [:,1] [:,M]



def get_bounds(gtype, LBi, UBi, LBj, UBj):
"""
  This routine computes grid bounds in the I- and J-directions.

  On Input:

     gtype      C-grid type. If zero, compute array allocation bounds.
                  Otherwise, compute bounds for IO processing.


  On Output:

     Itile      I-tile coordinate (a value from 0 to NtileI(ng)).
     Jtile      J-tile coordinate (a value from 0 to NtileJ(ng)).
     LBi        I-dimension Lower bound.
     UBi        I-dimension Upper bound.
     LBj        J-dimension Lower bound.
     UBj        J-dimension Upper bound.

"""



      LBi = 0
      LBj = 0
      UBi = Im + 1
      UBj = Jm + 1

      if (EWperiodic(ng)):
          LBi = LBi - NghostPoints
          UBi = UBi + NghostPoints - 1

      if (NSperiodic(ng)):
          LBj = LBj - NghostPoints
          UBj = UBj + NghostPoints - 1


def get_domain (gtype, Nghost, epsilon, Lfullgrid, Xmin, Xmax, Ymin, Ymax):
'''
This routine computes tile minimum and maximum fractional grid
coordinates.

On Input:

   Nghost     Number of ghost-points in the halo region:
                Nghost = 0,  compute non-overlapping coordinates.
                Nghost > 0,  compute overlapping bounds.
   gtype      C-grid type
   epsilon    Small value to add to Xmax and Ymax when the tile
                is lying on the eastern and northern boundaries
                of the grid. This is usefull when processing
                observations.
   Lfullgrid  Switch to include interior and boundaries points
                (TRUE) or just interior points (FALSE).

On Output:

   Xmin       Minimum tile fractional X-coordinate.
   Xmax       Maximum tile fractional X-coordinate.
   Ymin       Minimum tile fractional Y-coordinate.
   Ymax       Maximum tile fractional Y-coordinate.
'''

# Computes tile minimum and maximum fractional-grid coordinates.
# -----------------------------------------------------------------------

def get_bounds(gtype, Nghost, Imin, Imax, Jmin, Jmax):

# Include interior and boundary points.

      IF (Lfullgrid) THEN
        IF ((Itile.eq.0).and.                                           &
     &      ((gtype.eq.r2dvar).or.(gtype.eq.r3dvar).or.                 &
     &       (gtype.eq.v2dvar).or.(gtype.eq.v3dvar))) THEN
          Xmin=REAL(Imin,r8)
        ELSE
          Xmin=REAL(Imin,r8)-0.5_r8
        END IF
        IF (Itile.eq.(NtileI(ng)-1)) THEN
          IF ((gtype.eq.u2dvar).or.(gtype.eq.u3dvar)) THEN
            Xmax=REAL(Imax,r8)-0.5_r8
          ELSE
            Xmax=REAL(Imax,r8)
          END IF
        ELSE
          Xmax=REAL(Imax,r8)+0.5_r8
        END IF
        IF ((Jtile.eq.0).and.                                           &
     &      ((gtype.eq.r2dvar).or.(gtype.eq.r3dvar).or.                 &
     &       (gtype.eq.u2dvar).or.(gtype.eq.u3dvar))) THEN
          Ymin=REAL(Jmin,r8)
        ELSE
          Ymin=REAL(Jmin,r8)-0.5_r8
        END IF
        IF (Jtile.eq.(NtileJ(ng)-1)) THEN
          IF ((gtype.eq.v2dvar).or.(gtype.eq.v3dvar)) THEN
            Ymax=REAL(Jmax,r8)-0.5_r8
          ELSE
            Ymax=REAL(Jmax,r8)
          END IF
        ELSE
          Ymax=REAL(Jmax,r8)+0.5_r8
        END IF



      ELSE
         # Include only interior points.
        IF (Itile.eq.0) THEN
          IF ((gtype.eq.u2dvar).or.(gtype.eq.u3dvar)) THEN
             Xmin=REAL(Imin,r8)
          ELSE
             Xmin=REAL(Imin,r8)+0.5_r8
          END IF
        ELSE
          Xmin=REAL(Imin,r8)-0.5_r8
        END IF
        IF (Itile.eq.(NtileI(ng)-1)) THEN
          IF ((gtype.eq.u2dvar).or.(gtype.eq.u3dvar)) THEN
            Xmax=REAL(Imax,r8)-1.0_r8
          ELSE
            Xmax=REAL(Imax,r8)-0.5_r8
          END IF
        ELSE
          Xmax=REAL(Imax,r8)+0.5_r8
        END IF
        IF (Jtile.eq.0) THEN
          IF ((gtype.eq.v2dvar).or.(gtype.eq.v3dvar)) THEN
            Ymin=REAL(Jmin,r8)
          ELSE
            Ymin=REAL(Jmin,r8)+0.5
          END IF
        ELSE
          Ymin=REAL(Jmin,r8)-0.5_r8
        END IF
        IF (Jtile.eq.(NtileJ(ng)-1)) THEN
          IF ((gtype.eq.v2dvar).or.(gtype.eq.v3dvar)) THEN
            Ymax=REAL(Jmax,r8)-1.0_r8
          ELSE
            Ymax=REAL(Jmax,r8)-0.5_r8
          END IF
        ELSE
          Ymax=REAL(Jmax,r8)+0.5_r8
        END IF
      END IF
!
!  If tile lie at the grid eastern or northen boundary, add provided
!  offset value to allow processing at those boundaries.
!
      IF (Itile.eq.(NtileI(ng)-1)) THEN
        Xmax=Xmax+epsilon
      END IF
      IF (Jtile.eq.(NtileJ(ng)-1)) THEN
        Ymax=Ymax+epsilon
      END IF



def get_domain_edges (Eastern_Edge,     Western_Edge,     Northern_Edge,    Southern_Edge,
                      NorthEast_Corner, NorthWest_Corner, SouthEast_Corner, SouthWest_Corner,
                      NorthEast_Test,   NorthWest_Test,   SouthEast_Test,   SouthWest_Test):
    '''
    This routine sets the logical switches (T/F) needed for processing
    model variables in tiles adjacent to the domain boundary edges. It
    facilitates complicated nesting configurations.

    On Input:

    On Output:

       Eastern_Edge      tile next to the domain eastern  boundary
       Western_Edge      tile next to the domain western  boundary
       Northern_Edge     tile next to the domain northern boundary
       Southern_Edge     tile next to the domain southern boundary

       NorthEast_Corner  tile next to the domain northeastern corner
       NorthWest_Corner  tile next to the domain northwestern corner
       SouthEast_Corner  tile next to the domain southeastern corner
       SouthWest_Corner  tile next to the domain southwestern corner

       NorthEast_Test    test for tiles in the northeastern corner
       NorthWest_Test    test for tiles in the northwestern corner
       SouthEast_Test    test for tiles in the southeastern corner
       SouthWest_Test    test for tiles in the southwestern corner
    '''




  # Set switches for the full grid (tile=-1) to TRUE, since it contains
  # all the boundary edges and corners.  This is a special case use for
  # other purposes and need only in routine "var_bounds".

        # Western_Edge=.TRUE.
        # Eastern_Edge=.TRUE.
        # Southern_Edge=.TRUE.
        # Northern_Edge=.TRUE.
        #
        # SouthWest_Test=.TRUE.
        # SouthEast_Test=.TRUE.
        # NorthWest_Test=.TRUE.
        # NorthEast_Test=.TRUE.
        #
        # SouthWest_Corner=.TRUE.
        # SouthEast_Corner=.TRUE.
        # NorthWest_Corner=.TRUE.
        # NorthEast_Corner=.TRUE.

        WE HAVE TO REMOVE ALL THESE VARIABLES.



def get_iobounds():
'''
!  This routine computes the horizontal lower bound, upper bound, and  !
!  grid size for IO (NetCDF) variables.                               !
!                                                                      !
!                                                                      !
!  On Output, the horizontal lower/upper bounds and grid size for      !
!  each variable type and nested grid number  are loaded into the      !
!  IOBOUNDS structure which is declared in module MOD_PARAM:           !
!                                                                      !
!   IOBOUNDS(ng) % ILB_psi     I-direction lower bound (PSI)           !
!   IOBOUNDS(ng) % IUB_psi     I-direction upper bound (PSI)           !
!   IOBOUNDS(ng) % JLB_psi     J-direction lower bound (PSI)           !
!   IOBOUNDS(ng) % JUB_psi     J-direction upper bound (PSI)           !
!                                                                      !
!   IOBOUNDS(ng) % ILB_rho     I-direction lower bound (RHO)           !
!   IOBOUNDS(ng) % IUB_rho     I-direction upper bound (RHO)           !
!   IOBOUNDS(ng) % JLB_rho     J-direction lower bound (RHO)           !
!   IOBOUNDS(ng) % JUB_rho     J-direction upper bound (RHO)           !
!                                                                      !
!   IOBOUNDS(ng) % ILB_u       I-direction lower bound (U)             !
!   IOBOUNDS(ng) % IUB_u       I-direction upper bound (U)             !
!   IOBOUNDS(ng) % JLB_u       J-direction lower bound (U)             !
!   IOBOUNDS(ng) % JUB_u       J-direction upper bound (U)             !
!                                                                      !
!   IOBOUNDS(ng) % ILB_v       I-direction lower bound (V)             !
!   IOBOUNDS(ng) % IUB_v       I-direction upper bound (V)             !
!   IOBOUNDS(ng) % JLB_v       J-direction lower bound (V)             !
!   IOBOUNDS(ng) % JUB_v       J-direction upper bound (V)             !
!                                                                      !
!   IOBOUNDS(ng) % xi_psi      Number of I-direction points (PSI)      !
!   IOBOUNDS(ng) % xi_rho      Number of I-direction points (RHO)      !
!   IOBOUNDS(ng) % xi_u        Number of I-direction points (U)        !
!   IOBOUNDS(ng) % xi_v        Number of I-direction points (V)        !
!                                                                      !
!   IOBOUNDS(ng) % eta_psi     Number of J-direction points (PSI)      !
!   IOBOUNDS(ng) % eta_rho     Number of J-direction points (RHO)      !
!   IOBOUNDS(ng) % eta_u       Number of I-direction points (U)        !
!   IOBOUNDS(ng) % eta_v       Number of I-direction points (V)        !
'''

# !-----------------------------------------------------------------------
# !  Set IO lower/upper bounds and grid size for each C-grid type
# !  variable.
# !-----------------------------------------------------------------------
# !
# !  Recall that in non-nested applications the horizontal range,
# !  including interior and boundary points, for all variable types
# !  are:
# !
# !    PSI-type      [xi_psi, eta_psi] = [1:Lm(ng)+1, 1:Mm(ng)+1]
# !    RHO-type      [xi_rho, eta_rho] = [0:Lm(ng)+1, 0:Mm(ng)+1]
# !    U-type        [xi_u,   eta_u  ] = [1:Lm(ng)+1, 0:Mm(ng)+1]
# !    V-type        [xi_v,   eta_v  ] = [0:Lm(ng)+1, 1:Mm(ng)+1]
!
    IOBOUNDS.ILB_psi = 1
    IOBOUNDS.IUB_psi = Lm + 1
    IOBOUNDS.JLB_psi = 1
    IOBOUNDS.JUB_psi = Mm + 1

    IOBOUNDS.ILB_rho = 0
    IOBOUNDS.IUB_rho = Lm + 1
    IOBOUNDS.JLB_rho = 0
    IOBOUNDS.JUB_rho = Mm + 1

    IOBOUNDS.ILB_u   = 1
    IOBOUNDS.IUB_u   = Lm + 1
    IOBOUNDS.JLB_u   = 0
    IOBOUNDS.JUB_u   = Mm + 1

    IOBOUNDS.ILB_v   = 0
    IOBOUNDS.IUB_v   = Lm + 1
    IOBOUNDS.JLB_v   = 1
    IOBOUNDS.JUB_v   = Mm + 1


    # Set IO NetCDF files horizontal dimension size. Recall that NetCDF
    # does not support arrays with zero index as an array element.
!
      IOBOUNDS.IorJ    = BOUNDS.UBij - BOUNDS.LBij + 1
!
      IOBOUNDS(ng) % xi_psi  = IOBOUNDS(ng) % IUB_psi -                 &
     &                         IOBOUNDS(ng) % ILB_psi + 1
      IOBOUNDS(ng) % xi_rho  = IOBOUNDS(ng) % IUB_rho -                 &
     &                         IOBOUNDS(ng) % ILB_rho + 1
      IOBOUNDS(ng) % xi_u    = IOBOUNDS(ng) % IUB_u   -                 &
     &                         IOBOUNDS(ng) % ILB_u   + 1
      IOBOUNDS(ng) % xi_v    = IOBOUNDS(ng) % IUB_v   -                 &
     &                         IOBOUNDS(ng) % ILB_v   + 1
!
      IOBOUNDS(ng) % eta_psi = IOBOUNDS(ng) % JUB_psi -                 &
     &                         IOBOUNDS(ng) % JLB_psi + 1
      IOBOUNDS(ng) % eta_rho = IOBOUNDS(ng) % JUB_rho -                 &
     &                         IOBOUNDS(ng) % JLB_rho + 1
      IOBOUNDS(ng) % eta_u   = IOBOUNDS(ng) % JUB_u   -                 &
     &                         IOBOUNDS(ng) % JLB_u   + 1
      IOBOUNDS(ng) % eta_v   = IOBOUNDS(ng) % JUB_v   -                 &
     &                         IOBOUNDS(ng) % JLB_v   + 1



