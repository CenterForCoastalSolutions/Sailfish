import mod_param
import mod_grid
import mod_ncparam
import cupy as cp


def ana_grid(LBi, UBi, LBj, UBj, GRID):
    """
    This routine sets model grid using an analytical expressions.

    On Output:  stored in common blocks:

                             "grid"    (file grid.h)
                             "scalars" (file scalar.h)

       f        Coriolis parameter (1/seconds) at RHO-points.
       h        Bathymetry (meters; positive) at RHO-points.
       hmin     Minimum depth of bathymetry (m).
       hmax     Maximum depth of bathymetry (m).
       pm       Coordinate transformation metric "m" (1/meters)
                associated with the differential distances in XI
                at RHO-points.
       pn       Coordinate transformation metric "n" (1/meters)
                associated with the differential distances in ETA.
                at RHO-points.
       xp       XI-coordinates (m) at PSI-points.
       xr       XI-coordinates (m) at RHO-points.
       yp       ETA-coordinates (m) at PSI-points.
       yr       ETA-coordinates (m) at RHO-points.
    """



# Set grid parameters:
#    Xsize    Length (m) of domain box in the XI-direction.
#    Esize    Length (m) of domain box in the ETA-direction.
#    depth    Maximum depth of bathymetry (m).
#    f0       Coriolis parameter, f-plane constant (1/s).
#    beta     Coriolis parameter, beta-plane constant (1/s/m).



    if analiticalCase == 'BASIN':
          Xsize=3600.0E+03
          Esize=2800.0E+03
          depth=5000.0
          f0=1.0E-04
          beta=2.0E-11
    elif analiticalCase == 'BENCHMARK':
          Xsize=360.0              # degrees of longitude
          Esize=20.0               # degrees of latitude
          depth=4000.0
          f0=-1.0E-04
          beta=2.0E-11
    elif analiticalCase == 'BL_TEST':
          Xsize=100.0E+03
          Esize=5.0E+03
          depth=47.5
          f0=9.25E-04
          beta=0.0
    elif analiticalCase == 'CHANNEL':
          Xsize=600.0E+03
          Esize=360.0E+03
          depth=500.0
          f0=1.0E-04
          beta=0.0
    elif analiticalCase == 'CANYON':
          Xsize=128.0E+03
          Esize=96.0E+03
          depth=4000.0
          f0=1.0E-04
          beta=0.0
    elif analiticalCase == 'COUPLING_TEST':
          Xsize=6000.0*Lm
          Esize=6000.0*Mm
          depth=1500.0
          f0=5.0E-05
          beta=0.0
    elif analiticalCase == 'DOUBLE_GYRE':
          Xsize=1000.0E+03
          Esize=2000.0E+03
          depth=500.0
          f0=7.3E-05
          beta=2.0E-11
    elif analiticalCase == 'ESTUARY_TEST':
          Xsize=100000.0
          Esize=300.0
          depth=10.0
          f0=0.0
          beta=0.0
    elif analiticalCase == 'KELVIN':
          Xsize=20000.0*Lm
          Esize=20000.0*Mm
          depth=100.0
          f0=1.0E-04
          beta=0.0
    elif analiticalCase == 'FLT_TEST':
          Xsize=1.0E+03*Lm
          Esize=1.0E+03*Mm
          depth=10.0
          f0=0.0
          beta=0.0
    elif analiticalCase == 'GRAV_ADJ':
          Xsize=64.0E+03
          Esize=2.0E+03
          depth=20.0
          f0=0.0
          beta=0.0
    elif analiticalCase == 'LAB_CANYON':
          Xsize=0.55                  # width of annulus
          Esize=2.0*pi                # azimuthal length (radians)
          f0=4.0*pi/25.0
          beta=0.0
    elif analiticalCase == 'LAKE_SIGNELL':
          Xsize=50.0e3
          Esize=10.0e3
          depth=18.0
          f0=0.0E-04
          beta=0.0
    elif analiticalCase == 'LMD_TEST':
          Xsize=100.0E+03
          Esize=100.0E+03
          depth=50.0
          f0=1.09E-04
          beta=0.0
    elif analiticalCase == 'MIXED_LAYER':
          Xsize=500.0
          Esize=400.0
          depth=50.0
          f0=0.0
          beta=0.0
    elif analiticalCase == 'OVERFLOW':
          Xsize=4.0E+03
          Esize=200.0E+03
          depth=4000.0
          f0=0.0
          beta=0.0
    elif analiticalCase == 'RIVERPLUME1':
          Xsize=58.5E+03
          Esize=201.0E+03
          depth=150.0
          f0=1.0E-04
          beta=0.0
    elif analiticalCase == 'RIVERPLUME2':
          Xsize=100.0E+03
          Esize=210.0E+03
          depth=190.0
          f0=1.0E-04
          beta=0.0
    elif analiticalCase == 'SEAMOUNT':
          Xsize=320.0E+03
          Esize=320.0E+03
          depth=5000.0
          f0=1.0E-04
          beta=0.0
    elif analiticalCase == 'SOLITON':
          Xsize=48.0
          Esize=16.0
          depth=1.0
          f0=0.0
          beta=1.0
          g=1.0
    elif analiticalCase == 'SED_TEST1':
          Xsize=300.0
          Esize=36.0
          depth=10.0
          f0=0.0
          beta=0.0
    elif analiticalCase == 'SED_TOY':
          Xsize=40.0
          Esize=30.0
          depth=0.5
          f0=0.0
          beta=0.0
    elif analiticalCase == 'SHOREFACE':
          Xsize=1180.0
          Esize=140.0
          depth=15.0
          f0=0.0E-04
          beta=0.0
    elif analiticalCase == 'TEST_CHAN':
          Xsize=10000.0
          Esize=1000.0
          depth=10.0
          f0=0.0
          beta=0.0
    elif analiticalCase == 'UPWELLING':
          Xsize=1000.0*Lm
          Esize=1000.0*Mm
          depth=150.0
          f0=-8.26E-05
          beta=0.0
    elif analiticalCase == 'WEDDELL':
          Xsize=4000.0*Lm
          Esize=4000.0*Mm
          depth=4500.0
          f0=0.0
          beta=0.0
    elif analiticalCase == 'WINDBASIN':
          Xsize=2000.0*Lm
          Esize=1000.0*Mm
          depth=50.0
          f0=1.0E-04
          beta=0.0
    else:
        errorMsg('ana_grid.py: analytical case "%s" not found.' % analiticalCase)




    # Compute the (XI,ETA) coordinates at PSI- and RHO-points.
    # Set grid spacing (m).
    # -----------------------------------------------------------------------

    # Determine I- and J-ranges for computing grid data.  These ranges
    # are special in periodic boundary conditons since periodicity cannot
    # be imposed in the grid coordinates.

    Imin = Istr - 1
    Imax = Iend + 1
    Jmin = Jstr - 1
    Jmax = Jend + 1

    I, J = cp.meshgrid(cp.arange(Jmin, Jmax+1), cp.arange(Imin, Imax+1))

    dx = Xsize/Lm
    dy = Esize/Mm

    if analiticalCase == 'BENCHMARK':
        # Spherical coordinates set-up.

        spherical = True

        lonr = dx*(I - 0.5)
        latr = dy*(J - 0.5) - 70.0
        lonu = dx*I
        lonp = lonu
        latu = latr
        lonv = lonr
        latv = dy*J - 70.0
        latp = latv

    elif analiticalCase == 'LAB_CANYON':
        # Polar coordinates set-up.
        dth = 0.01                          # Azimulthal spacing
        cff = 4.0*pi/(dth*Mm) - 1.0         # F

        r = 0.35 + dx*(I - 1)
        theta = -pi + 0.5*dth*((cff + 1.0)*(J - 1      ) + (cff - 1.0)*(Mm/(2*pi))*cp.sin(2*pi*(J - 1      )/Mm))
        xp = r*cp.cos(theta)
        yp = r*cp.sin(theta)

        r = 0.35 + dx*(I - 1 + 0.5)
        theta = -pi + 0.5*dth*((cff + 1.0)*(J - 1 + 0.5) + (cff - 1.0)*(Mm/(2*pi))*cp.sin(2*pi*(J - 1 + 0.5)/Mm))
        xr = r*cp.cos(theta)
        yr = r*cp.sin(theta)

        xu = xp
        yu = yr

        xv = xr
        yv = yp

    else:

        if analiticalCase == 'BL_TEST':
            dx = 0.5*4000.0*(Lm + 1)*I + 675.0

        xp = dx*(I - 1)
        yp = dy*(J - 1)

        xr = xp + 0.5*dx
        yr = yp + 0.5*dy

        xu = xp
        yu = yr

        xv = xr
        yv = yp


    # Report statistics.
    if SPHERICAL:
        stats_2dfld(iNLM, p2dvar, GRID, (lonp, latp), ('longitude of PSI-points: lon_psi', 'latitude of PSI-points: lat_psi'))
        stats_2dfld(iNLM, r2dvar, GRID, (lonr, latr), ('longitude of RHO-points: lon_rho', 'latitude of RHO-points: lat_rho'))
        stats_2dfld(iNLM, u2dvar, GRID, (lonu, latu), ('longitude of   U-points: lon_u',   'latitude of   U-points: lat_u'))
        stats_2dfld(iNLM, v2dvar, GRID, (lonp, latp), ('longitude of   V-points: lon_v',   'latitude of   V-points: lat_v'))

    else:
        stats_2dfld(iNLM, p2dvar, GRID, (xp, yp), ('x-location of PSI-points: x_psi', 'y-location of PSI-points: y_psi'))
        stats_2dfld(iNLM, r2dvar, GRID, (xr, yr), ('x-location of RHO-points: x_rho', 'y-location of RHO-points: y_rho'))
        stats_2dfld(iNLM, u2dvar, GRID, (xu, yu), ('x-location of   U-points: x_u',   'y-location of   U-points: y_u'))
        stats_2dfld(iNLM, v2dvar, GRID, (xp, yp), ('x-location of   V-points: x_v',   'y-location of   V-points: y_v'))



    # Compute coordinate transformation metrics at RHO-points "pm" and
    # "pn"  (1/m) associated with the differential distances in XI and
    # ETA, respectively.
    # -----------------------------------------------------------------------


    if analiticalCase == 'BENCHMARK':

        # Spherical coordinates set-up.

        valX = Lm/(2.0*pi*Eradius)
        valY = 360*Mm/(2.0*pi*Eradius*Esize)

        pm[:,:] = valX/cp.cos(cp.deg2rad(-70.0 + dy*J - 0.5))
        pn[:,:] = valY

    elif analiticalCase == 'LAB_CANYON':

        # Polar coordinates set-up.

        r = 0.35 + dx*((I-1) + 0.5)
        theta = 0.5*dth*((cff + 1.0) +  (cff - 1.0)*cp.cos(2*pi*(J - 1)/Mm))
        pm[:,:] = 1.0/dx
        pn[:,:] = 1.0/(r*theta)

    else:
        if analiticalCase == 'BL_TEST':
              dx = 0.5*(4000.0/(Lm + 1))*I + 675.0

        pm[:,:] = 1.0/dx
        pn[:,:] = 1.0/dy


    # Report statistics.
    stats_2dfld(iNLM, r2dvar, Stats[9:10], GRID, [pm, pn], ['reciprocal XI-grid spacing: pm', 'reciprocal ETA-grid spacing: pn'])


    # Exchange boundary data.
    exchange_r2d(LBi, UBi, LBj, UBj, (pm, pn))


    XXXXif analiticalCase == 'CURVGRID' and analiticalCase == 'UV_ADV' :
    # Compute d(1/n)/d(xi) and d(1/m)/d(eta) at RHO-points.
    # -----------------------------------------------------------------------


    dndx(i,j) = Dξ(1.0/pn)
    dmde(i,j) = Dη(1.0/pm)


    # Report statistics.
    stats_2dfld(iNLM, r2dvar, GRID, (dmde, dndx), ('reciprocal XI-grid spacing: pm', 'reciprocal ETA-grid spacing: pn'))


    # Exchange boundary data.
    exchange_r2d(LBi, UBi, LBj, UBj, (dndx, dmde))


    # Angle (radians) between XI-axis and true EAST at RHO-points.
    # -----------------------------------------------------------------------

    IT, JT = cp.meshgrid(cp.arange(JstrT, JendT + 1), cp.arange(IstrT, IendT + 1))
    if analiticalCase == 'LAB_CANYON':

          theta = -pi + 0.5*dth*((cff+1.0)*((JT - 1) + 0.5) + (cff - 1.0)*(Mm/(2*pi))*cp.sin(2*pi*((JT - 1) + 0.5)/Mm))
          angler[:,:] = theta

    elif WEDDELL:
          angler[:,:] = 90.0*deg2rad

    else:
          angler[:,:] = 0.0

    # Report Statistics.
    stats_2dfld(iNLM, r2dvar, GRID, (angler), ('angle between XI-axis and EAST: angler'))

    # Exchange boundary data.
    exchange_r2d(LBi, UBi, LBj, UBj, angler)


    # Compute Coriolis parameter (1/s) at RHO-points.
    # -----------------------------------------------------------------------

    if analiticalCase == 'BENCHMARK':
        val1 = 2.0*(2.0*pi*366.25/365.25)/86400.0
        f[:,:] = val1*cp.sin(latr*deg2rad)

    elif analiticalCase == 'WEDDELL':
        val1 = 10.4/Lm
        f[:,:] = 2.0*7.2E-05*cp.sin((-79.0 + (IT - 1)*val1)*deg2rad)

    else:
        val1 = 0.5*Esize
        f[:,:] = f0 + beta*(yr - val1)


    # Report Statistics.
    stats_2dfld(iNLM, r2dvar, GRID, (f), ('Coriolis parameter at RHO-points: f'))

    # Exchange boundary data.
    exchange_r2d(LBi, UBi, LBj, UBj, f)



    # Set bathymetry (meters; positive) at RHO-points.
    # -----------------------------------------------------------------------


    if analiticalCase == 'BENCHMARK':
          h[:,:] = 500.0 + 1750.0*(1.0+cp.tanh((68.0 + latr)/dy))

    elif analiticalCase == 'BL_TEST'':
          val1 = (xr + 500.0)/15000.0
          h[:,:] = 14.0 + 25.0*(1.0 - cp.exp(-pi*xr*1.0E-05)) - 8.0*cp.exp(-val1*val1)

    elif analiticalCase == 'CANYON':
          val1 = 32000.0 - 16000.0*(cp.sin(pi*xr)/Xsize))**24
          h[:,:] = 20.0 + 0.5*(depth - 20.0)*(1.0 + cp.tanh((yr - val1)/10000.0))

    elif analiticalCase == 'ESTUARY_TEST':
          h[:,:] = 5.0 + (Xsize-xr)/Xsize*5.0

    elif analiticalCase == 'LAB_CANYON':
        r = 0.35 + dx*(i - 1 + 0.5)
        theta = -pi +  0.5*dth*((cff + 1.0)*(JT - 1 + 0.5) + (cff - 1.0)*(Mm/(2*pi))*cp.sin(dth*(JT - 1 + 0.5)/Mm))
        val1 = 0.55 # r_small
        val2 = 0.15 # lambda
        if cp.abs(theta) < 0.181818181818:
            val1 -= 0.15*(cp.cos(pi*theta*0.55/0.2)**2)
            val2 += 0.15*(cp.cos(pi*theta*0.55/0.2)**2)

        h[:,:] = 0.125 - 0.1*(cp.cos(0.5*pi*(r - val1)/val2)**2)
        h[r <= val1] = 0.025  # shelf
        h[r >=  0.7] = 0.125  # deep

    elif analiticalCase == 'LAKE_SIGNELL':
        h[:,:] = 18.0 - 16.0*(Mm - JT)/(Mm - 1)


    elif analiticalCase == 'MIXED_LAYER':
        h[:,:] = 50.0

    elif analiticalCase == 'OVERFLOW':
        val1 = 200.0
        h[:,:] = val1 + 0.5*(depth-val1)*(1.0 + cp.tanh((yr - 100000.0)/20000.0))


    elif analiticalCase == 'RIVERPLUME1 or RIVERPLUME2':
        h[:, :5] = 15.0
        h[:, 5:] = depth + (Lm - IT)*(15.0 - depth)/(Lm - 6)

    elif analiticalCase == 'SEAMOUNT':
        val1 = (xr - 0.5*Xsize)/40000.0
        val2 = (yr - 0.5*Esize)/40000.0
        h = depth - 4500.0*cp.exo(-(val1*val1 + val2*val2))

    elif analiticalCase == 'SED_TOY':
        h[:,:] = 20.0

    elif analiticalCase == 'SHOREFACE':
        h[:,:] = 11.75 - 0.0125*Xsize/(Lm+1)*IT

    elif analiticalCase == 'TEST_CHAN':
        h[:,:] = 10.0 + 0.4040*IT/(Lm+1)

    elif analiticalCase == 'UPWELLING':
        if NSperiodic:
            val1 = IT
            val1[IT > Lm/2] = Lm + 1 - IT
            h = cp.minimum(depth, 84.5 + 66.526*cp.tanh((val1 - 10.0)/7.0))
        elif EWperiodic:
            val1 = JT
            val1[JT > Mm/2] = Mm + 1 - JT
            h = cp.minimum(depth, 84.5 + 66.526*cp.tanh((val1 - 10.0) / 7.0))

    elif analiticalCase == 'WEDDELL':
        val1 = 98.80
        val2 = 0.8270
        for k in range(-1,27):
        # DO k=-1,26
            xwrk[k] = (k-1)*15.0*1000.0
            hwrk[k] = 375.0

        for k in range(27,233):
            zwrk = -2.0 + (k - 1)*0.020
            xwrk[k] = (520.0 + val1 + zwrk*val1 + val1*val2*LOG(COSH(zwrk)))*1000.0
            hwrk[k] = -75.0 + 2198.0*(1.0 + val2*TANH(zwrk))

        for k in range(233,236):
            xwrk[k] = (850.0 + (k - 228)*50.0)*1000.0
            hwrk[k] = 4000.0

        h[:,:] = 375.0
        for k in range(1,235):
            mWrk =  ((xwrk[k] <= xr) and (xr < xwrk[k+1]))
            cff = 1.0/(xwrk[k+1] - xwrk[k])
            h[mWrk] = cff*(xwrk[k+1]-xr(i,j))*hwrk[k] + cff*(xr(i,j)-xwrk[k])*hwrk[k+1]

    elif analiticalCase == 'WINDBASIN':
        ival = int(0.03*(Lm + 1))
        mIval1 = (IT < ival)
        mIval2 = ((Lm + 1 - IT) < ival)
        val1[:,:] = 1.0
        val1[mIval1] = 1.0 - ((IT + 1 - ival)/ival)**2
        val1[mIval2] = 1.0 - ((Lm + 1 - IT - ival)/ival)**2

        val2 = 2.0*(JT - (Mm + 1))/2/(Mm+1)
        h = depth*(0.08 + 0.92*val1*(1.0 - val2*val2))

    else:
        h[:.:] = depth


    # Report Statistics.
    stats_2dfld(iNLM, r2dvar, GRID, (h), ('bathymetry at RHO-points: h'))

    # Exchange boundary data.
    exchange_r2d(LBi, UBi, LBj, UBj, h)



