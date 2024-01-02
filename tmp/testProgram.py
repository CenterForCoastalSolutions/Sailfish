import sys
import cupy as cp

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
from set_coords import set_coords   # TODO I believe this can be inside mod_grid

from main2d import main2d
from main3d import main3d

import set_weights

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


# Generates the mesh.
ana_grid.ana_grid('Basin', GRID)
GRID.updateMetrics()

set_coords(GRID)



v = OCEAN.zeta[0,:,:]
# bc_2d.bc_r2d([v], BOUNDARY)


# zetabc(OCEAN.zeta_t2, compTimes, BOUNDARY)
# barotropicVelocityBC(OCEAN.ubar_t2, OCEAN.vbar_t2, compTimes, BOUNDARY)


# Print report of all input parameters read.
input.printReport()


set_weights.set_weights(compTimes)

mod_operators.initModule(GRID)
# mod_operators.initOperators((1,), (1,), (10, *(GRID.h.shape)))

from mod_operators import grsz, bksz, set_depth
set_depth(grsz, bksz, (GRID.Vtransform, OCEAN.Zt_avg1, GRID.z_w, GRID.z_r, GRID.h, GRID.hc, GRID.Hz,
                       GRID.sc_r,  GRID.sc_w, GRID.Cs_r, GRID.Cs_w))
OCEAN.AKv[:,:,:] =+0.00001
for i in range(GRID.N-2):
    OCEAN.AKv[i,:,:] = -0.3*((GRID.z_r[i,:,:])/(GRID.z_r[0])**2)+0.00001


# rhs3d(GRID, OCEAN, BOUNDARY)

main3d(compTimes, GRID, OCEAN, BOUNDARY)
main2d(compTimes, GRID, OCEAN, BOUNDARY)

# cProfile.run('main2d(compTimes,  GRID, OCEAN, BOUNDARY)')
pass