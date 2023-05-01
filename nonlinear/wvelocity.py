import geometry
import cupy as cp

# This subroutines computes vertical velocity (m/s) at W-points       !
# from the vertical mass flux (omega*hz/m*n).  This computation       !
# is done solely for output purposes.                                 !

def wvelocity (ng, Ninp):
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

    # Exchange time-averaged fields.
    exchange_u2d(ng, DU_avg1)
    exchange_v2d(ng, DV_avg1)

    # Compute contribution due to quasi-horizontal motions along
    # S-coordinate surfaces:  U·GRADs(z) = (u*i + v*j)·GRADs(z).
    vert = UtoR(u[:,:,:,Ninp]*dξUtoR(z_r, pm)) + VtoR(v[:,:,:,Ninp]*dηVtoR(z_r, pn))


    # Compute contribution due to time tendency of the free-surface,
    # d(zeta)/d(t), which is the vertical velocity at the free-surface
    # and it is expressed in terms of barotropic mass flux divergence.
    # Notice that it is divided by the total depth of the water column.
    # This is needed because this contribution is linearly distributed
    # throughout the water column by multiplying it by the distance from
    # the bottom to the depth at which the vertical velocity is computed.

    wcDepth = z_w[:,:,-1] - z_w[:, :, 0]  # Water Column (wc) depth
    lcDepth = z_w[:,:,: ] - z_w[:, :, 0]  # Local (lc) depth

    tmp = W - ηξdivR(DU_avg1)*lcDepth/wcDepth

    if OMEGA_IMPLICIT:
        tmp += Wi

    wvel = pm*pn*tmp + RtoW(vert)

    # ///////////////////////////////////  /

    krnl(u, v, z_r, pm, pn)

    // Depth of the water Column (wc)
    const double wcDepth = z_w(N - 1, 0, 0) - z_w(0, 0, 0)

    #pragma unroll 1
    for (int k=1; k<N; n++)
    {
        auto dξz_rU = dξRtoU(z_r);
        auto dηz_rV = dηRtoV(z_r);
        vertW(0,0,k) = RtoW(UtoR(u(0,0,k)*dξz_rU(0,0,k)) +
                            VtoR(v(0,0,k)*dηz_rV(0,0,k)));
    }

    # Compute contribution due to time tendency of the free-surface,
    # d(zeta)/d(t), which is the vertical velocity at the free-surface
    # and it is expressed in terms of barotropic mass flux divergence.
    # Notice that it is divided by the total depth of the water column.
    # This is needed because this contribution is linearly distributed
    # throughout the water column by multiplying it by the distance from
    # the bottom to the depth at which the vertical velocity is computed.

    // Local depth (lc)
    double lcDepth = z_w(k,0,0) - z_w(0,0,0)

    double tmpW = W - UVdivR(DU_avg1)*lcDepth/wcDepth

    if (OMEGA_IMPLICIT) tmpW += Wi;

    wvel = pm*pn*tmpW + vertW;



    # Set lateral boundary conditions.
    bc_w3d(ng, N, wvel)

