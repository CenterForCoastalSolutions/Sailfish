import geometry
import cupy as cp

# This subroutines computes vertical velocity (m/s) at W-points       !
# from the vertical mass flux (omega*hz/m*n).  This computation       !
# is done solely for output purposes.                                 !

def wvelocity (Ninp):
    '''
      Compute "true" vertical velocity (m/s).

          In ROMS, the terrain-following vertical velocity, omega, is given by:

                 Hz * omega = w - dz/dt - div(z)

          where w is the "true" vertical velocity and

                 div(z) = pm*u*(dz/dξ) + pn*v*(dz/dη)

          The vertical coordinate is a function of several parameters but only
          the free-surface is time dependent. However, in sediment applications
          with stratigraphy, the bathymetry (h) also evolves in time.
    '''


    # Compute contribution due to quasi-horizontal motions along S-coordinate surfaces:  U·GRADs(z).
    vert = UtoR(u[Ninp,:,:,:]*dξRtoU(z_r)) + VtoR(v[Ninp,:,:,:]*DξRtoU(z_r))


    # Compute contribution due to time tendency of the free-surface, d(zeta)/d(t), which is the vertical velocity at the free-surface
    # and it is expressed in terms of barotropic mass flux divergence.
    # Notice that it is divided by the total depth of the water column.
    # This is needed because this contribution is linearly distributed throughout the water column by multiplying it by the distance from
    # the bottom to the depth at which the vertical velocity is computed.

    wcDepth = z_w[:,:,-1] - z_w[:, :, 0]  # Water Column (wc) depth (2D)
    lcDepth = z_w[:,:,: ] - z_w[:, :, 0]  # Local (lc) depth (3D)

    tmp = W - ηξdivR(DU_avg1)*lcDepth/wcDepth


    wvel = pm*pn*tmp + RtoW(vert)


    krnl(u, v, z_r, pm, pn)


    for (K=1; K<N; K++):
        dξz_rU = dξRtoU(z_r)
        dηz_rV = dηRtoU(z_r)

        vertW = RtoW(UtoR(u*dξz_rU) + VtoR(v*dηz_rV))

    # Compute contribution due to time tendency of the free-surface, d(zeta)/d(t), which is the vertical velocity at the free-surface
    # and it is expressed in terms of barotropic mass flux divergence.
    # Notice that it is divided by the total depth of the water column.
    # This is needed because this contribution is linearly distributed throughout the water column by multiplying it by the distance from
    # the bottom to the depth at which the vertical velocity is computed.

    tmpW = W - UVdivR(DU_avg1)*lcDepth/wcDepth


    wvel = pm*pn*tmpW + vertW



    # Set lateral boundary conditions.
    bc_w3d(N, wvel)

