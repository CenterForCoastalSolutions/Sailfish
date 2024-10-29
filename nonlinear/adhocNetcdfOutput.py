from netCDF4 import Dataset

def create_roms_netcdf(filename, xi_rho, eta_rho, s_rho):
    # Dimension relationships
    xi_u = xi_rho - 1
    eta_u = eta_rho
    xi_v = xi_rho
    eta_v = eta_rho - 1
    s_w = s_rho + 1

    # Create a new NetCDF file
    ncfile = Dataset(filename, 'w', format='NETCDF4')

    # Define dimensions
    ncfile.createDimension('xi_rho', xi_rho)
    ncfile.createDimension('xi_u', xi_u)
    ncfile.createDimension('xi_v', xi_v)
    ncfile.createDimension('eta_rho', eta_rho)
    ncfile.createDimension('eta_u', eta_u)
    ncfile.createDimension('eta_v', eta_v)
    ncfile.createDimension('s_rho', s_rho)
    ncfile.createDimension('s_w', s_w)
    ncfile.createDimension('ocean_time', None)  # Unlimited dimension
    ncfile.createDimension('N', s_rho)  # Alias for s_rho

    # Define variables
    ntimes = ncfile.createVariable('ntimes', 'i4')
    ntimes.long_name = "number of long time-steps"

    xl = ncfile.createVariable('xl', 'f8')
    xl.long_name = "domain length in the XI-direction"
    xl.units = "meter"

    el = ncfile.createVariable('el', 'f8')
    el.long_name = "domain length in the ETA-direction"
    el.units = "meter"

    # Define s_rho, s_w, Cs_r, and Cs_w
    s_rho_var = ncfile.createVariable('s_rho', 'f8', ('s_rho',))
    s_rho_var.long_name = "S-coordinate at RHO-points"
    s_rho_var.valid_min = -1.0
    s_rho_var.valid_max = 0.0
    s_rho_var.positive = "up"
    s_rho_var.standard_name = "ocean_s_coordinate_g2"
    s_rho_var.formula_terms = "s: s_rho C: Cs_r eta: zeta depth: h depth_c: hc"
    s_rho_var.field = "s_rho, scalar"

    s_w_var = ncfile.createVariable('s_w', 'f8', ('s_w',))
    s_w_var.long_name = "S-coordinate at W-points"
    s_w_var.valid_min = -1.0
    s_w_var.valid_max = 0.0
    s_w_var.positive = "up"
    s_w_var.standard_name = "ocean_s_coordinate_g2"
    s_w_var.formula_terms = "s: s_w C: Cs_w eta: zeta depth: h depth_c: hc"
    s_w_var.field = "s_w, scalar"

    Cs_r = ncfile.createVariable('Cs_r', 'f8', ('s_rho',))
    Cs_r.long_name = "S-coordinate stretching curves at RHO-points"
    Cs_r.valid_min = -1.0
    Cs_r.valid_max = 0.0
    Cs_r.field = "Cs_r, scalar"

    Cs_w = ncfile.createVariable('Cs_w', 'f8', ('s_w',))
    Cs_w.long_name = "S-coordinate stretching curves at W-points"
    Cs_w.valid_min = -1.0
    Cs_w.valid_max = 0.0
    Cs_w.field = "Cs_w, scalar"

    # Bathymetry and coordinate variables
    h = ncfile.createVariable('h', 'f8', ('eta_rho', 'xi_rho'))
    h.long_name = "bathymetry at RHO-points"
    h.units = "meter"
    h.grid = "grid"
    h.location = "face"
    h.coordinates = "x_rho y_rho"
    h.field = "bath, scalar"

    x_rho = ncfile.createVariable('x_rho', 'f8', ('eta_rho', 'xi_rho'))
    x_rho.long_name = "x-locations of RHO-points"
    x_rho.units = "meter"
    x_rho.field = "x_rho, scalar"

    y_rho = ncfile.createVariable('y_rho', 'f8', ('eta_rho', 'xi_rho'))
    y_rho.long_name = "y-locations of RHO-points"
    y_rho.units = "meter"
    y_rho.field = "y_rho, scalar"

    x_u = ncfile.createVariable('x_u', 'f8', ('eta_u', 'xi_u'))
    x_u.long_name = "x-locations of U-points"
    x_u.units = "meter"
    x_u.field = "x_u, scalar"

    y_u = ncfile.createVariable('y_u', 'f8', ('eta_u', 'xi_u'))
    y_u.long_name = "y-locations of U-points"
    y_u.units = "meter"
    y_u.field = "y_u, scalar"

    x_v = ncfile.createVariable('x_v', 'f8', ('eta_v', 'xi_v'))
    x_v.long_name = "x-locations of V-points"
    x_v.units = "meter"
    x_v.field = "x_v, scalar"

    y_v = ncfile.createVariable('y_v', 'f8', ('eta_v', 'xi_v'))
    y_v.long_name = "y-locations of V-points"
    y_v.units = "meter"
    y_v.field = "y_v, scalar"

    # Time and ocean variables
    ocean_time = ncfile.createVariable('ocean_time', 'f8', ('ocean_time',))
    ocean_time.long_name = "time since initialization"
    ocean_time.units = "seconds since 2001-01-01 00:00:00"
    ocean_time.calendar = "proleptic_gregorian"
    ocean_time.field = "time, scalar, series"

    zeta = ncfile.createVariable('zeta', 'f4', ('ocean_time', 'eta_rho', 'xi_rho'))
    zeta.long_name = "free-surface"
    zeta.units = "meter"
    zeta.time = "ocean_time"
    zeta.grid = "grid"
    zeta.location = "face"
    zeta.coordinates = "x_rho y_rho ocean_time"
    zeta.field = "free-surface, scalar, series"

    ubar = ncfile.createVariable('ubar', 'f4', ('ocean_time', 'eta_u', 'xi_u'))
    ubar.long_name = "vertically integrated u-momentum component"
    ubar.units = "meter second-1"
    ubar.time = "ocean_time"
    ubar.grid = "grid"
    ubar.location = "edge1"
    ubar.coordinates = "x_u y_u ocean_time"
    ubar.field = "ubar-velocity, scalar, series"

    vbar = ncfile.createVariable('vbar', 'f4', ('ocean_time', 'eta_v', 'xi_v'))
    vbar.long_name = "vertically integrated v-momentum component"
    vbar.units = "meter second-1"
    vbar.time = "ocean_time"
    vbar.grid = "grid"
    vbar.location = "edge2"
    vbar.coordinates = "x_v y_v ocean_time"
    vbar.field = "vbar-velocity, scalar, series"

    u = ncfile.createVariable('u', 'f4', ('ocean_time', 's_rho', 'eta_u', 'xi_u'))
    u.long_name = "u-momentum component"
    u.units = "meter second-1"
    u.time = "ocean_time"
    u.grid = "grid"
    u.location = "edge1"
    u.coordinates = "x_u y_u s_rho ocean_time"
    u.field = "u-velocity, scalar, series"

    v = ncfile.createVariable('v', 'f4', ('ocean_time', 's_rho', 'eta_v', 'xi_v'))
    v.long_name = "v-momentum component"
    v.units = "meter second-1"
    v.time = "ocean_time"
    v.grid = "grid"
    v.location = "edge2"
    v.coordinates = "x_v y_v s_rho ocean"

    return ncfile