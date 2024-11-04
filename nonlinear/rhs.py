import mod_ocean
import mod_boundary
from mod_operators import horizontalAdvection, verticalAdvection, addCoriolis, vertIntegral

# This subroutine evaluates right-hand-side terms for 3D momentum and tracers equations.
import mod_grid
import mod_ocean


# ηξ Δ




def rhs3d(GRID, OCEAN, BOUNDARY):
    from mod_operators import grsz, bksz

    # TODO: Check this BC
    BC = BOUNDARY.zetaBC.bcIdxFieldIdx2


    CURVGRID = False
    UV_COR = True
    UV_ADV = True


    # Initialize computations for new time step of the 3D primitive variables.
#    pre_step3d()

    # Compute baroclinic pressure gradient.
#    prsgrd()


    # Compute right-hand-side terms for the 3D momentum equations.
    # -----------------------------------------------------------------------

    OCEAN.ru_t2[:] = 0.0
    OCEAN.rv_t2[:] = 0.0

    # Add in Coriolis terms.
    if UV_COR:
        addCoriolis(grsz, bksz, (GRID.fomn, OCEAN.u_t2, OCEAN.v_t2, OCEAN.ru_t2, OCEAN.rv_t2, GRID.Hz))   # TODO: Check that all must be _t2.


    # These are the terms (u*d_xi(1/n) - v*d_eta(1/m))*Hz*u  and  (v*d_xi(1/n) - u*d_eta(1/m))*Hz*v
    # Where u is computed using an upwind interpolation.
    if CURVGRID and UV_ADV:
        # Add in curvilinear transformation terms. (for the advection terms only)
        addCurvedGridTerms()



    pass
    # TODO: Add in nudging towards 3D momentum climatology.



    # Compute horizontal advection of momentum.
    # -----------------------------------------------------------------------
    # TODO: I had to reduce the number of threads because an error related to GPU's limited resources. This has to be done in a better way.
    horizontalAdvection((grsz[0]*GPUMUL,), (bksz[0]//GPUMUL,), (OCEAN.u_t2, OCEAN.v_t2, OCEAN.Huon.ravel(), OCEAN.Hvom.ravel(), OCEAN.ru_t2, OCEAN.rv_t2, BC))



    if UV_ADV:

        # Add in vertical advection.
        # -----------------------------------------------------------------------

        verticalAdvection((grsz[0]*GPUMUL,), (bksz[0]//GPUMUL,), (OCEAN.u_t2, OCEAN.v_t2, OCEAN.W, OCEAN.ru_t2, OCEAN.rv_t2, BC))

          # TODO: remember


        # Compute forcing term for the 2D momentum equations.
        # -----------------------------------------------------------------------

        # Vertically integrate baroclinic right-hand-side terms. If not ***body force stresses***, add in the difference between surface and
        # bottom stresses.

        # TODO: Bring this back
        # vertIntegral(OCEAN.ru_t2, rufrc)
        # vertIntegral(OCEAN.rv_t2, rvfrc)

        # if BODY_FORCE:
        #     rufrc += om_u*on_u*(sustr - bustr)
        #     rvfrc += om_v*on_v*(svstr - bvstr)





