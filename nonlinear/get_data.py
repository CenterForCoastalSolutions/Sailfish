import mod_param
import mod_boundary
import mod_forces
import mod_grid
import mod_iounits
import mod_ncparam
import mod_scalars
import mod_stepping
from misc import *


def get_data(BOUNDS, BOUNDARY):
    """This routine reads in forcing, climatology and other data from NetCDF files.  If there is more than one
    time-record,  data is loaded into global  two-time  record arrays. The interpolation is carried elsewhere.
    """

    # At this point we are only using analytical functions, so this function is empty.
    msgInfo('**WARNING**, function get_data is empty because we are using only analytical functions')

    return


