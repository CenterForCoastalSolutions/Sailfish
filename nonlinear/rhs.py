import mod_ocean
import mod_boundary
from mod_operators import horizontalAdvection, verticalAdvection, addCoriolis

# This subroutine evaluates right-hand-side terms for 3D momentum and tracers equations.
import mod_grid
import mod_ocean


# ηξ Δ




def rhs3d(GRID, OCEAN, BOUNDARY):
    from mod_operators import grsz, bksz

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
        addCoriolis(grsz, bksz, (GRID.fomn, OCEAN.u_t2, OCEAN.v_t2, OCEAN.ru_t2, OCEAN.rv_t2))   # TODO: Check that all must be _t2.





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
    horizontalAdvection(grsz, bksz, (OCEAN.u_t2, OCEAN.v_t2, OCEAN.u_t2, OCEAN.u_t2, OCEAN.ru_t2, OCEAN.rv_t2, BC, Gadv, 16))
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

        verticalAdvection(grsz, bksz, (0,0,0,0,BC,0))
        # TODO: remember to put the correct parameters Huvn... instead of 0


        # Compute forcing term for the 2D momentum equations.
        # -----------------------------------------------------------------------

        # Vertically integrate baroclinic right-hand-side terms. If not body force stresses, add in the difference between surface and
        # bottom stresses.

        rufrc += vertIntegral(ru)
        rvfrc += vertIntegral(rv)

        # DO i=IstrU,Iend
        #   rufrc(i,j) = ru(i,j,1,nrhs)
        # END DO
        #
        # DO k=2,N(ng)
        #   DO i=IstrU,Iend
        #     rufrc(i,j) = rufrc(i,j) + ru(i,j,k,nrhs)
        #   END DO
        # END DO
        #
        # DO i=IstrU,Iend
        #     rufrc(i,j) += om_v(i,j)*on_v(i,j)*(sustr(i,j) - bustr(i,j))
        #
        #
        # rvfrc(i,j)=rv(i,j,1,nrhs)
        # DO i=Istr,Iend
        #       rvfrc(i,j)=rvfrc(i,j) + rv(i,j,k,nrhs)
        #
        #   DO i=Istr,Iend
        #     rvfrc(i,j) += om_v(i,j)*on_v(i,j)*(svstr(i,j) - bvstr(i,j))
        #
        # END IF


