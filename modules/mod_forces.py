import mod_param

#  Surface momentum stresses.

# sustr        Surface momentum flux (wind stress) in the ξ-direction (m²/s²) at horizontal U-points.         !
# svstr        Surface momentum flux (wind stress) in the η-direction (m²/s²) at horizontal V-points.        !

# sustrG       Latest two-time snapshots of input "sustr" grided data used for interpolation.                         !
# svstrG       Latest two-time snapshots of input "svstr" grided data used for interpolation.

# Taux         Surface stress in the ξ-direction at rho points from atmospheric model.                                      !
# Tauy         Surface stress in the η-direction at rho points from atmospheric model.


#  Bottom momentum stresses.

# bustr        Bottom momentum flux (bottom stress) in the ξ-direction (m²/s²) at horizontal U-points.         !
# bvstr        Bottom momentum flux (bottom stress) in the η-direction (m²/s²) at horizontal V-points.         !


# Surface wind induced waves.

# Hwave        Surface wind induced wave height (m).
!  HwaveG       Latest two-time snapshots of input "Hwave" grided data used for interpolation.                         !
!  Dwave        Surface wind induced mean wave direction (radians).
!  DwaveG       Latest two-time snapshots of input "Dwave" grided data used for interpolation.                         !
!  Dwavep       Surface wind induced peak wave direction (radians).
!  DwavepG      Latest two-time snapshots of input "Dwavep" grided data used for interpolation.                         !
!  Lwave        Mean surface wavelength read in from swan output
!  LwaveG       Latest two-time snapshots of input "Lwave" grided data used for interpolation.                         !
!  Lwavep       Peak surface wavelength read in from swan output
!  LwavepG      Latest two-time snapshots of input "Lwavep" grideddata used for interpolation.                         !
!  Pwave_top    Wind induced surface wave period (s).
!  Pwave_topG   Latest two-time snapshots of input "Pwave_top" grided data used for interpolation.                         !
!  Pwave_bot    Wind induced bottom wave period (s).
!  Pwave_botG   Latest two-time snapshots of input "Pwave_bot" grided data used for interpolation.                         !
!  Uwave_rms    Bottom orbital velocity read in from swan output
!  Uwave_rmsG   Latest two-time snapshots of input "Uwave_rms" grided data used for interpolation.                         !
!  wave_dissip  Wave dissipation
!  wave_dissipG Latest two-time snapshots of input "wave_dissip" gridded data used for interpolation.                 !
!  Wave_break   Percent of wave breaking for use with roller model.
!  Wave_breakG  Latest two-time snapshots of input "wave_break" gridded data used for interpolation.                 !
!  Wave_ds      Wave directional spreading.
!  Wave_qp      Wave spectrum peakedness.
!
!  Solar shortwave radiation flux.
!
!  srflx        Surface shortwave solar radiation flux (degC m/s) at horizontal RHO-points                            !
!  srflxG       Latest two-time snapshots of input "srflx" grided data used for interpolation.                         !

!  Cloud fraction.                                                     !
!                                                                      !
!  cloud        Cloud fraction (percentage/100).                       !
!  cloudG       Latest two-time snapshots of input "cloud" grided      !
!                 data used for interpolation.                         !
!                                                                      !
!  Surface heat fluxes, Atmosphere-Ocean bulk parameterization.        !
!                                                                      !
!  lhflx        Latent heat flux (degC m/s).                           !
!  lrflx        Longwave radiation (degC m/s).                         !
!  shflx        Sensible heat flux (degC m/s).                         !
!                                                                      !
!  Surface air humidity.                                               !
!                                                                      !
!  Hair         Surface air specific (g/kg) or relative humidity       !
!                 (percentage).                                        !
!  HairG        Latest two-time snapshots of input "Hair" grided       !
!                 data used for interpolation.                         !
!                                                                      !
!  Surface air pressure.                                               !
!                                                                      !
!  Pair         Surface air pressure (mb).                             !
!  PairG        Latest two-time snapshots of input "Pair" grided       !
!                 data used for interpolation.                         !
!                                                                      !
!  Surface air temperature.                                            !
!                                                                      !
!  Tair         Surface air temperature (Celsius)                      !
!  TairG        Latest two-time snapshots of input "Tair" grided       !
!                 data used for interpolation.                         !
!  PotT         Surface air potential temperature (Kelvin)             !
!  Surface Winds.                                                      !
!                                                                      !
!  Uwind        Surface wind in the XI-direction (m/s) at              !
!                 horizontal RHO-points.                               !
!  UwindG       Latest two-time snapshots of input "Uwind" grided      !
!                 data used for interpolation.                         !
!  Vwind        Surface wind in the ETA-direction (m/s) at             !
!                 horizontal RHO-points.                               !
!  VwindG       Latest two-time snapshots of input "Vwind" grided      !
!                 data used for interpolation.                         !
!                                                                      !
!  Rain fall rate.                                                     !
!                                                                      !
!  evap         Evaporation rate (kg/m2/s).                            !
!  rain         Rain fall rate (kg/m2/s).                              !
!  rainG        Latest two-time snapshots of input "rain" grided       !
!                 data used for interpolation.                         !
!                                                                      !
!  Snow fall rate.                                                     !
!                                                                      !
!  snow         Snow fall rate (kg/m2/s).                              !
!  snowG        Latest two-time snapshots of input "snow" grided       !
!                 data used for interpolation.                         !
!                                                                      !
!  Surface tracer fluxes.                                              !
!                                                                      !
!  stflux       Forcing surface flux of tracer type variables from     !
!                 data, coupling, bulk flux parameterization, or       !
!                 analytical formulas.                                 !
!                                                                      !
!                 stflux(:,:,itemp)  surface net heat flux             !
!                 stflux(:,:,isalt)  surface net freshwater flux (E-P) !
!                                                                      !
!  stfluxG      Latest two-time snapshots of input "stflux" grided     !
!                 data used for interpolation.                         !
!                                                                      !
!  stflx        ROMS state surface flux of tracer type variables       !
!                 (TracerUnits m/s) at horizontal RHO-points, as used  !
!                 in the governing equations.                          !

                                                              !
!  Bottom tracer fluxes.                                               !
!                                                                      !
!  btflux       Forcing bottom flux of tracer type variables from      !
!                 data or analytical formulas. Usually, the bottom     !
!                 flux of tracer is zero.                              !
!                                                                      !
!                 btflux(:,:,itemp)  bottom heat flux                  !
!                 btflux(:,:,isalt)  bottom freshwater flux            !
!                                                                      !
!  btfluxG      Latest two-time snapshots of input "vtflux" grided     !
!                 data used for interpolation.                         !
!                                                                      !
!  btflx        ROMS state bottom flux of tracer type variables        !
!                 (TracerUnits m/s) at horizontal RHO-points, as used  !
!                 in the governing equations.                          !
!                                                                      !
!  Surface heat flux correction.                                       !
!                                                                      !
!  dqdt         Surface net heat flux sensitivity to SST,              !
!                 d(Q)/d(SST), (m/s).                                  !
!  dqdtG        Latest two-time snapshots of input "dqdt" grided       !
!                 data used for interpolation.                         !
!  sst          Sea surface temperature (Celsius).                     !
!  sstG         Latest two-time snapshots of input "sst" grided        !
!                 data used for interpolation.                         !
!                                                                      !
!  Surface freshwater flux correction.                                 !
!                                                                      !
!  sss          Sea surface salinity (PSU).                            !
!  sssG         Latest two-time snapshots of input "sss" grided        !
!                 data used for interpolation.                         !
!  sssflx       Sea surface salinity flux correction.                  !
!  sssflxG      Latest two-time snapshots of input "sssflx" grided     !
!                 data used for interpolation.                         !
!                                                                      !
!  Surface spectral downwelling irradiance.                            !
!                                                                      !
!  SpecIr       Spectral irradiance (NBands) from 400-700 nm at        !
!                 5 nm bandwidth.                                      !
!  avcos        Cosine of average zenith angle of downwelling          !
!                 spectral photons.                                    !
!                                                                      !


class Forces:

    # Nonlinear model state.

    real(r8), pointer :: sustr(:,:)
    real(r8), pointer :: svstr(:,:)

    real(r8), pointer :: bustr(:,:)
    real(r8), pointer :: bvstr(:,:)









forces = Forces()




def initialize_forces (LBi, UBi, LBj, UBj, model):
# This routine initialize all variables in the module using first
# touch distribution policy. In shared-memory configuration, this
# operation actually performs propagation of the  "shared arrays"
# across the cluster, unless another policy is specified to
# override the default.

    IniVal = 0.0

#include "set_bounds.h"
!
!  Set array initialization range.
!
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        Imin=BOUNDS(ng)%LBi(tile)
      ELSE
        Imin=Istr
      END IF
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        Imax=BOUNDS(ng)%UBi(tile)
      ELSE
        Imax=Iend
      END IF
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        Jmin=BOUNDS(ng)%LBj(tile)
      ELSE
        Jmin=Jstr
      END IF
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        Jmax=BOUNDS(ng)%UBj(tile)
      ELSE
        Jmax=Jend
      END IF

!  Initialize module variables.
!-----------------------------------------------------------------------
!
!  Nonlinear model state.
!
    if (model == 0) or (model == iNLM):
        FORCES.bustr[Jmin:Jmax, IMin:Imax] = IniVal
        FORCES.buvtr[Jmin:Jmax, IMin:Imax] = IniVal


      RETURN
      END SUBROUTINE initialize_forces

