
# -----------------------------------------------------------------------
#   Compute open segments cross-section area and mass flux.
# -----------------------------------------------------------------------

    my_area = 0.0
    my_flux = 0.0


# west
# cff=0.5_r8*(zeta(Istr-1,j,kinp)+h(Istr-1,j)+zeta(Istr  ,j,kinp)+h(Istr  ,j))*on_u(Istr,j)
# my_flux=my_flux+cff*ubar(Istr,j,kinp)
# east
# cff=0.5_r8*(zeta(Iend  ,j,kinp)+h(Iend  ,j)+zeta(Iend+1,j,kinp)+h(Iend+1,j))*on_u(Iend+1,j)
# my_flux=my_flux-cff*ubar(Iend+1,j,kinp)
# south
# cff=0.5_r8*(zeta(i,Jstr-1,kinp)+h(i,Jstr-1)+zeta(i,Jstr, kinp)+h(i,Jstr   ))*om_v(i,Jstr)
# my_flux=my_flux+cff*vbar(i,JstrV-1,kinp)
# north
# cff=0.5_r8*(zeta(i,Jend  ,kinp)+h(i,Jend  )+zeta(i,Jend+1,kinp)+h(i,Jend+1))*om_v(i,Jend+1)
# my_flux=my_flux-cff*vbar(i,Jend+1,kinp)


    REDUCTION SUM over the whole boundary
        cff = RtoBC(zeta(0,0, kinp) + h(0,0), idx)*compBC(on_u(0,0), om_v(0,0), idx);

        bc_area += cff;
        bc_flux += cff*compBCsign(ubar, vbar, idx);


    ubar_xs = bc_flux/bc_area;

SUBROUTINE set_DUV_bc_tile (ng)
#-----------------------------------------------------------------------
# Set vertically integrated mass fluxes "Duon" and "Dvom" along
#   the open boundaries in such a way that the integral volume is
#   conserved.  This is done by applying "ubar_xs" correction to
#   the velocities.
# -----------------------------------------------------------------------


# define I_RANGE MAX(2,IstrU-1),MIN(Iend+1,Lm(ng))
# define J_RANGE MAX(2,JstrV-1),MIN(Jend+1,Mm(ng))



compBC(Duon, Dvon, idx) = 0.5*RtoBC(Drhs, idx)*BCcomp(ubar(0,0,kimp), vbar(0,0,kimp), idx) + BCsign(idx)*ubar_xs)*compBC(on_u, on_v, idx);
# W: Duon(Istr,j)  =0.5_r8*(Drhs(Istr,j)+Drhs(Istr-1,j))*(ubar(Istr,  j,kinp)-ubar_xs)*on_u(Istr,j)
# E: Duon(Iend+1,j)=0.5_r8*(Drhs(Iend+1,j)+Drhs(Iend,j))*(ubar(Iend+1,j,kinp)+ubar_xs)*on_u(Iend+1,j)
# S: Dvom(i,Jstr)  =0.5_r8*(Drhs(i,Jstr)+Drhs(i,Jstr-1))*(vbar(i,Jstr,  kinp)-ubar_xs)*om_v(i,Jstr)
# N: Dvom(i,Jend+1)=0.5_r8*(Drhs(i,Jend+1)+Drhs(i,Jend))*(vbar(i,Jend+1,kinp)+ubar_xs)*om_v(i,Jend+1)



def conserve_mass_tile (ng):


# -----------------------------------------------------------------------
#   Corrects velocities across the open boundaries to enforce global
#   mass conservation constraint.
# -----------------------------------------------------------------------

    BCcomp(ubar, idx) = BCcomp(ubar, idx) + BCsign(idx)*ubar_xs

    #
    # W: ubar(Istr  ,j,kinp)=(ubar(Istr,j,  kinp)-ubar_xs)
    # E: ubar(Iend+1,j,kinp)=(ubar(Iend+1,j,kinp)+ubar_xs)
    # S: vbar(i,Jstr,  kinp)=(vbar(i,Jstr,kinp)  -ubar_xs)
    # N: vbar(i,Jend+1,kinp)=(vbar(i,Jend+1,kinp)+ubar_xs)



