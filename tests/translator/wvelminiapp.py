import cupy as cp
from cupyx.time import repeat

N = 11
Istr = 1
Jstr = 1
Iend = 1000
Jend = 1000
r8 = 8

LBi = Istr - 1
UBi = Iend + 1
LBj = Jstr - 1
UBj = Jend + 1

pm = cp.zeros((UBj - LBj + 1, UBi - LBi + 1), dtype=cp.float64)
pn = cp.zeros((UBj - LBj + 1, UBi - LBi + 1), dtype=cp.float64)
z_r = cp.zeros((N, UBj - LBj + 1, UBi - LBi + 1), dtype=cp.float64)
u = cp.zeros((2, N, UBj - LBj + 1, UBi - LBi + 1), dtype=cp.float64)
v = cp.zeros((2, N, UBj - LBj + 1, UBi - LBi + 1), dtype=cp.float64)

vert = cp.zeros((N, UBj - LBj, UBi - LBi ), dtype=cp.float64)
wrk = cp.zeros((N, UBj - LBj , UBi - LBi ), dtype=cp.float64)

Ninp = 1
res = 0

# Initialization
k, j, i = cp.indices(z_r.shape)
print(z_r.shape)
for l in range(2):
    u[l, :, :, :] = i - j - k * 2.0 - 2
    v[l, :, :, :] = i - j + k * 2.0 + 2
    pm[:, :] = i[0, :, :] + 1
    pn[:, :] = j[0, :, :] + 1
    z_r[:, :, :] = j - i - k - 1
print(z_r.sum(), u.sum(), v.sum(), pm.sum(), pn.sum())

def trivialCupy():
    res = 0.0
    for iter in range(10):

        wrk[:, :-1, :] = u[Ninp, :, 1:-1, 1:] * (z_r[:, 1:-1, 1:] - z_r[:, 1:-1, :-1]) * cp.tile(pm[1:-1, 0:-1] + pm[1:-1, 1:], [N,1,1])
        vert[:,:,1:] = 0.25*(wrk[:,:, 1:] + wrk[:,:,:-1])

        wrk[:, :, :-1] = v[Ninp, :, 1:, 1:-1] * (z_r[:, 1:, 1:-1] - z_r[:, :-1, 1:-1]) * cp.tile(pn[0:-1, 1:-1] + pn[1:, 1:-1], [N,1,1])
        vert[:, :-1, 1:] += 0.25 * (wrk[:, 1:, :-1] + wrk[:, :-1, :-1])

        res += vert[:,:-1,:].sum(axis=None)
        print(iter, res)

@cp.fuse()
def fuseRtoU(u1, u2):
    return 0.25*(u1 + u2)


@cp.fuse()
def fuseRtoV(u1, u2):
    return 0.25*(u1 + u2)

@cp.fuse()
def fuseDηRtoU(u1, u2):
    return (u1 - u2)

@cp.fuse()
def fuseDyRtoV(u1, u2):
    return (u1 - u2)



def fuseAndCupy():
    res = 0.0
    for iter in range(10):
        wrk[:, :-1, :] = u[Ninp, :, 1:-1, 1:] * fuseDηRtoU(z_r[:, 1:-1, 1:], z_r[:, 1:-1, :-1]) * (pm[None,1:-1, 0:-1] + pm[None,1:-1, 1:])
        vert[:,:,1:] = fuseRtoU(wrk[:,:, 1:], wrk[:,:,:-1])

        wrk[:, :, :-1] = v[Ninp, :, 1:, 1:-1] * fuseDyRtoV(z_r[:, 1:, 1:-1], z_r[:, :-1, 1:-1]) * (pn[None, 0:-1, 1:-1] + pn[None, 1:, 1:-1])
        vert[:, :-1, 1:] += fuseRtoV(wrk[:, 1:, :-1] , wrk[:, :-1, :-1])

        # res += vert[:,:-1,:].sum(axis=None)
    print(iter, res)

kernelfunc = cp.ElementwiseKernel(
    'float64 z_r2d, raw float64 z_r, raw float64 u, raw float64 v, raw float64 pn, raw float64 pm',
    'float64 vert',
'''
    double wrk;
    wrk = u[i]*(z_r[i] - z_r[i-1])*(pm[i-1]+pm[i])
//   printf("%i ", i)    ;
// _ind.size()    ;
//  vert=1.0
''',
    'kernelfunc', options=('--dopt=on', '--use_fast_math'))
# wrk[:, :-1, :] = u[Ninp, :, 1:-1, 1:] * (z_r[:, 1:-1, 1:] - z_r[:, 1:-1, :-1]) * cp.tile(pm[1:-1, 0:-1] + pm[1:-1, 1:], [N,1,1])
#         vert[:,:,1:] = 0.25*(wrk[:,:, 1:] + wrk[:,:,:-1])
#
#         wrk[:, :, :-1] = v[Ninp, :, 1:, 1:-1] * (z_r[:, 1:, 1:-1] - z_r[:, :-1, 1:-1]) * cp.tile(pn[0:-1, 1:-1] + pn[1:, 1:-1], [N,1,1])
#         vert[:, :-1, 1:] += 0.25 * (wrk[:, 1:, :-1] + wrk[:, :-1, :-1])
#
#         res += vert[:,:-1,:].sum(axis=None)


# def kernelCupy():
#     res = 0.0
#     for iter in range(100):
#
#         # vert = kernelfunc(z_r, z_r, u[Ninp,:,:,:], v[Ninp,:,:,:], pn, pm)
#         vert = fuseCupy(z_r)
#
# kernelCupy()
print(repeat(trivialCupy, (), n_repeat=10, n_warmup=10))
print(repeat(fuseAndCupy, (), n_repeat=10, n_warmup=10))

