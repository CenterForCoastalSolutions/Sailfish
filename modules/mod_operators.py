import mod_grid


grid = None

def initModule(GRID _grid):
    global grid, pm...
    grid = _grid
    pm = grid.pm


def RtoU(r):
    return r*pm


def divUVtoR(u, v):
    global grid
    something like: grid.pm*GRID.pn*(GRID.on_u*(ux - ux) + GRID.om_v*(vy-vy))