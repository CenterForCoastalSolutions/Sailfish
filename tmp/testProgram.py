import sys
import cProfile

# import mod_operators

sys.path.insert(0, '../modules')
sys.path.insert(0, '../utility')
sys.path.insert(0, '../nonlinear')
sys.path.insert(0, '../functionals')

import io_utils
import mod_grid
import mod_boundary
import mod_io
import mod_ocean
import mod_comptimes
import mod_physical_params
import mod_operators

from main2d import main2d

from barotropicVelocityBC import barotropicVelocityBC
from zetabc import zetabc

import ana_grid



# Reads the file with meta information about the variables.
varInfoList = io_utils.VarInfoList(path = r'../modules')

# Reads the input file.
input = io_utils.Input(r'../utility/test.in', varInfoList)

GRID           = mod_grid           .Grid(input)

BOUNDARY       = mod_boundary       .Boundary(input, GRID)
io             = mod_io             .T_IO(input)
physicalParams = mod_physical_params.PhysicalParams(input, GRID)
compTimes      = mod_comptimes      .CompTimes(input)
OCEAN          = mod_ocean          .T_OCEAN(input, GRID)

mod_operators.initModule(GRID)


v = OCEAN.zeta[0,:,:]
# bc_2d.bc_r2d([v], BOUNDARY)


# zetabc(OCEAN.zeta_t2, compTimes, BOUNDARY)
# barotropicVelocityBC(OCEAN.ubar_t2, OCEAN.vbar_t2, compTimes, BOUNDARY)


# Print report of all input parameters read.
input.printReport()



# Generates the mesh.
ana_grid.ana_grid('Basin', GRID)
GRID.updateMetrics()

# main2d(compTimes,  GRID, OCEAN, BOUNDARY)

cProfile.run('main2d(compTimes,  GRID, OCEAN, BOUNDARY)')
pass