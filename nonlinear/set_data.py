




import mod_param

def set_data():
    """
    This subroutine processes forcing, boundary, climatology, and other input data. It time-interpolates between snapshots.
    """
    import mod_param
    import mod_boundary
    import mod_forces
    import mod_grid
    import mod_ncparam
    import mod_ocean
    import mod_stepping
    import mod_scalars
    import ana_grid



    ana_grid.bulkFluxes()
    ana_grid.sources()
    ana_grid.kinematicSurfaceMomentumFlux()
    ana_grid.


    # implement later
    # do something with netCDF here

# Set point Sources/Sinks (river runoff).
# =======================================================================


# Point Source/Sink vertically integrated mass transport.

#
# Set forcing data.
# =======================================================================

# Set switch to process surface atmospheric fields.



#  if !defined BULK_FLUXES && !defined BULK_FLUXES2D

# Set kinematic surface momentum flux (m2/s2).

#


# Set surface air pressure (mb).

#  ifdef ANA_PAIR
        ana_pair(iNLM)

#  endif




# ======================================================================
#  Set open boundary conditions fields.
# ======================================================================
#  Set free-surface open boundary conditions.

# ifdef ANA_FSOBC
        ana_fsobc(iNLM)


# Set 2D momentum component open boundary conditions.

# ifdef ANA_M2OBC
        CALL ana_m2obc (ng, tile, iNLM)






