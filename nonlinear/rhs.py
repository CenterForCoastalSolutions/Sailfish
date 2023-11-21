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



    # Add in Coriolis terms.
    if UV_COR:
        addCoriolis(grsz, bksz, (GRID.fomn, OCEAN.u_t2, OCEAN.v_t2, OCEAN.ru_t2, OCEAN.rv_t2, GRID.Hz))   # TODO: Check that all must be _t2.





    # if defined CURVGRID && defined UV_ADV
    # These are the terms (u*d_xi(1/n) - v*d_eta(1/m))*Hz*u  and  (v*d_xi(1/n) - u*d_eta(1/m))*Hz*v
    # Where u is computed using an upwind interpolation.
    if CURVGRID:
        # Add in curvilinear transformation terms. (for the advection terms only)
        addCurvedGridTerms()





    pass
    # TODO: Add in nudging towards 3D momentum climatology.


    # Compute horizontal advection of momentum.
    # -----------------------------------------------------------------------
    horizontalAdvection(grsz, bksz, (OCEAN.u_t2, OCEAN.v_t2, OCEAN.Huon, OCEAN.Hvom, OCEAN.ru_t2, OCEAN.rv_t2, BC))
    # TODO: remember to put the correct parameters Huvn... instead of ru_t2
    # pass
        # mod_ocean.T_OCEAN.,  const double *_v, const double *_Huon, const double *_Hvom,
        #                  const double *_ru, const double *_rv, const int N)

    # uξξ  = Dξξ(   u[nrhs,:,:,:], BC)
    # Huξξ = Dξξ(Huon[nrhs,:,:,:], BC)
    # # DO j=Jstr,Jend
    # #   DO i=IstrUm1,Iendp1
    # #     uξξ(i,j)  = u(i-1,j,k,nrhs) - 2.0*u(i,j,k,nrhs) + u(i+1,j,k,nrhs)
    # #     Huξξ(i,j) = Huon(i-1,j,k)   - 2.0*Huon(i,j,k)   + Huon(i+1,j,k)
    # #   END DO
    # # END DO
    # #
    # #
    # # # BC's for 2nd derivatives
    # # IF (DOMAIN(ng)%Western_Edge(tile)) THEN
    # #     DO j=Jstr,Jend
    # #       uξξ (Istr,j) = uξξ (Istr+1,j)
    # #       Huξξ(Istr,j) = Huξξ(Istr+1,j)
    # #     END DO
    # # END IF
    # #
    # #
    # # IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
    # #     DO j=Jstr,Jend
    # #       uξξ (Iend+1,j) = uξξ (Iend,j)
    # #       Huξξ(Iend+1,j) = Huξξ(Iend,j)
    # #     END DO
    # # END IF
    pass


    if UV_ADV:

        # Add in vertical advection.
        # -----------------------------------------------------------------------

        verticalAdvection(grsz, bksz, (OCEAN.u_t2, OCEAN.v_t2, OCEAN.W, OCEAN.ru_t2, OCEAN.rv_t2, BC))
        # TODO: remember to put the correct parameters Huvn... instead of 0


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





