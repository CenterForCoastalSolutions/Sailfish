import mod_ocean
import mod_boundary
from mod_operators import horizontalAdvection, verticalAdvection

# This subroutine evaluates right-hand-side terms for 3D momentum and tracers equations.
import mod_grid
import mod_ocean


# ηξ Δ




def rhs3d(OCEAN, BOUNDARY):
    from mod_operators import grsz, bksz

    BC = BOUNDARY.zetaBC.bcIdxField


    CURVGRID = False

    # Initialize computations for new time step of the 3D primitive variables.
#    pre_step3d()

    # Compute baroclinic pressure gradient.
#    prsgrd()


    # Compute right-hand-side terms for the 3D momentum equations.
    # -----------------------------------------------------------------------


    # K_LOOP : DO k=1,N
    # K loop should go inside the kernel.




        # Add in Coriolis terms.

#        addCoriolis()
        # ifdef UV_COR
            # for j=JstrV-1,Jend:
            #   for i=IstrU-1,Iend:
            #     cff=Hz(i,j,k)*fomn(i,j)
            #     UFξ(i,j)=cff*0.5*(v(i,j,k,nrhs) + v(i,j+1,k,nrhs))
            #     VFη(i,j)=cff*0.5*(u(i,j,k,nrhs) + u(i+1,j,k,nrhs))

            # for j=Jstr,Jend:
            #   for i=IstrU,Iend:
            #     cff1=0.5*(UFx(i,j) + UFx(i-1,j))
            #     ru(i,j,k,nrhs) = ru(i,j,k,nrhs) + cff1
            #
            #
            # for j=JstrV,Jend:
            #   for i=Istr,Iend:
            #     cff1=0.5*(VFe(i,j) + VFe(i,j-1))
            #     rv(i,j,k,nrhs) = rv(i,j,k,nrhs) - cff1

        # endif



    pass
    # if defined CURVGRID && defined UV_ADV
    # These are the terms (u*d_xi(1/n) - v*d_eta(1/m))*Hz*u  and  (v*d_xi(1/n) - u*d_eta(1/m))*Hz*v
    # Where u is computed using an upwind interpolation.
    if CURVGRID:
        addCurvedGridTerms()
    # Add in curvilinear transformation terms. (for the advection terms only)
    # -----------------------------------------------------------------------



    # DO j=JstrV-1,Jend
    #   DO i=IstrU-1,Iend
    #     cff1=0.5*(v(i,j  ,k,nrhs) + v(i,j+1,k,nrhs))
    #     cff2=0.5*(u(i  ,j,k,nrhs) + u(i+1,j,k,nrhs))
    #     cff3=cff1*dndx(i,j)
    #     cff4=cff2*dmde(i,j)
    #
    #     cff=Hz(i,j,k)*(cff3-cff4)
    #     UFx(i,j)=cff*cff1
    #     VFe(i,j)=cff*cff2
    #   END DO
    # END DO
    #
    # DO j=Jstr,Jend
    #   DO i=IstrU,Iend
    #     cff1=0.5*(UFx(i,j)+UFx(i-1,j))
    #     ru(i,j,k,nrhs) = ru(i,j,k,nrhs) + cff1
    #   END DO
    # END DO
    #
    # DO j=JstrV,Jend
    #   DO i=Istr,Iend
    #     cff1=0.5*(VFe(i,j)+VFe(i,j-1))
    #     rv(i,j,k,nrhs) = rv(i,j,k,nrhs) - cff1
    #   END DO
    # END DO
# endif

    pass
    # Add in nudging towards 3D momentum climatology.
    pass
    # # -----------------------------------------------------------------------
    #     IF (LnudgeM3CLM(ng)) THEN
    #       DO j=Jstr,Jend
    #         DO i=IstrU,Iend
    #           cff=0.25*(CLIMA(ng)%M3nudgcof(i-1,j,k)  CLIMA(ng)%M3nudgcof(i  ,j,k))*om_u(i,j)*on_u(i,j)
    #           ru(i,j,k,nrhs)=ru(i,j,k,nrhs) + cff*(Hz(i-1,j,k)+Hz(i,j,k))*(CLIMA(ng)%uclm(i,j,k) - u(i,j,k,nrhs))
    #         END DO
    #       END DO
    #
    #       DO j=JstrV,Jend
    #         DO i=Istr,Iend
    #           cff=0.25*(CLIMA(ng)%M3nudgcof(i,j-1,k) + CLIMA(ng)%M3nudgcof(i,j  ,k))*om_v(i,j)*on_v(i,j)
    #           rv(i,j,k,nrhs)=rv(i,j,k,nrhs) + cff*(Hz(i,j-1,k)+Hz(i,j,k))*(CLIMA(ng)%vclm(i,j,k) - v(i,j,k,nrhs))
    #         END DO
    #       END DO
    #     END IF


    pass


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


# ifdef UV_ADV


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


#endif
