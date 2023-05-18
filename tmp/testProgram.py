import sys
sys.path.insert(0, '../modules')
sys.path.insert(0, '../utility')

import io_utils
import mod_grid
import mod_boundary
import mod_io
import mod_ocean
import mod_comptimes
import mod_physical_params
from zetabc import zetabc

import ana_grid



# Reads the file with meta information about the variables.
varInfoList = io_utils.VarInfoList(path = r'..\modules')

# Reads the input file.
input = io_utils.Input(r'..\utility\test.in', varInfoList)

GRID           = mod_grid     .Grid(input)
BOUNDARY       = mod_boundary .Boundary(input, GRID)
io             = mod_io       .T_IO(input)
physicalParams = read_phypar  .PhysicalParams(input, GRID)
compTimes      = mod_comptimes.CompTimes(input)
OCEAN          = mod_ocean    .T_OCEAN(input, GRID)


v = OCEAN.zeta[0,:,:]
# bc_2d.bc_r2d([v], BOUNDARY)


kout = 0
zetabc(OCEAN.zeta, kout, compTimes, BOUNDARY)
barotropicVelocityBC(OCEAN.ubar, OCEAN.vbar, kout, compTimes, BOUNDARY)


# Print report of all input parameters read.
input.printReport()



# Generates the mesh.
ana_grid.ana_grid('Basin', GRID)


pass