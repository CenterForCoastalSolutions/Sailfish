import datetime
import sys
import cupy as cp

cp.set_debug_mode(True)

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

OCEAN.AKv[:,:,:] =  2.0e-3 + 8.0e-3*cp.exp(GRID.z_w/150)


t0 = datetime.datetime.now()
main3d(compTimes, GRID, OCEAN, BOUNDARY)
t1 = datetime.datetime.now()
print('clock time = %.3f s' %(t1-t0).seconds)


# cProfile.run('main2d(compTimes,  GRID, OCEAN, BOUNDARY)')
pass