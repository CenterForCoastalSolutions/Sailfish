


# Adams-Moulton 3rd order coefficients
AM3_2 =  5.0 / 12.0
AM3_1 =  8.0 / 12.0
AM3_0 = -1.0 / 12.0


# HSIMT tracer advection coefficients for the TVD limiter (Wu and Zhu, 2010).
cc1 = 1.0/4.0
cc2 = 1.0/2.0
cc3 = 1.0/12.0



# Physical constants.
# -----------------------------------------------------------------------
Cp      = 3985.0         # Specific heat for seawater (Joules/Kg/degC).
Csolar  = 1353.0         # Solar irradiantion constant, 1360-1380 (W/m2).
Eradius = 6371315.0      # Earth equatorial radius (m).
StefBo  = 5.67E-8        # Stefan-Boltzmann constant (W/m2/K4).
emmiss  = 0.97           # Infrared emissivity (non dimensional)
rhow    = 1000.0         # fresh water density (kg/m3).
rho0    = 1025.0         # Mean density (Kg/m3) used when the Boussinesq approximation is inferred.
g       = 10.0           # Acceleration due to gravity (m/s2).
gorho0  = g/rho0         # gravity divided by mean density anomaly. m4/s2/kg
vonKar  = 0.41           # von Karman constant (non dimensional)


# Math constants.
# -----------------------------------------------------------------------
pi = 3.14159265358979323846
deg2rad = pi/180.0
rad2deg = 180.0/pi
day2sec = 86400.0


# Time constants.
# -----------------------------------------------------------------------
sec2day = 1.0/86400.0
jul_off = 2440000.0
