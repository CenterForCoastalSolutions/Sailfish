import cupy as cp


def WtoR(varW):
    """ Interpolates W points into Rho points"""

    return 0.5*(varW[:-1,:,:] + varW[:-1,:,:])


def UtoR(varU):
    """ Interpolates U points into Rho points"""

    return 0.5*(varW[:,:,-1:] + varW[:,:,1:])


def dηUtoR(R, pm):
    """ Computes the centered η derivative of U points into Rho points
        Remember that pm are coordinate transformation metrics equal to the inverse of the distance between points"""

    return (R[1:, :, :] - R[:-1, :, :]) / (pm[1:, :, :] + pm[:-1, :, :])


def dξVtoR(R, pn):
    """ Computes the centered ξ derivative of V points into Rho points
        Remember that pm are coordinate transformation metrics equal to the inverse of the distance between points"""

    return (R[:, 1:, :] - R[:, :-1, :]) / (pn[:, 1:, :] + pn[:, :-1, :])


def ηξdivR(U, V):
    return U[1:,:] - U[:-1,:] + V[:,1:] - V[:,-1]


def RtoW(R):
    ''' Uses a 4 or 3 order Lagrange polynomial to interpolate W points into R ones.
    It uses the 4th order approximation in the inner points for computing the approximation at the center using a scheme like:

             R        R        R        R
                          W
           1/16     9/16     9/16     1/16

    Near the bottom or the free surface, there are not enough points to do this. For the points of index 1 and N-2 it does:

                 0        1        2     (example for W at idx 1)
                 R        R        R
                     W
                3/8     3/4      1/8


    At the bottom or the free surface nodes there are two problems. 1. We have to extrapolate, 2. The locations depend on
    time and the coefficients must be recomputed at every time step.

                 0        1        2     (example for W at idx 1)
                 R        R        R
            W
                 ?        ?        ?

     '''


    W = cp.zeros([R.shape[0] + 1] + R.shape[1:])

    # Interior
    W[2:-2,:,:] = (1/16)*(R[0:-4,:,:] + R[3:-1,:,:]) + (9/16)*(R[1:-3,:,:] + R[2:-2,:,:])


    # 1 and N-2
    W[ 1] = (3/8)*R[ 0,:,:] + (3/4)*R[ 1,:,:] + (1/8)*R[ 2,:,:]
    W[-2] = (3/8)*R[-1,:,:] + (3/4)*R[-2,:,:] + (1/8)*R[-3,:,:]


    # Free surface and Bottom
    slopeT = (z_r[:,:, 1] - z_w[:,:, 0])/(z_r[:,:, 1] - z_r[:,:, 0])  # extrapolation slope Top
    slopeB = (z_w[:,:,-1] - z_r[:,:,-1])/(z_r[:,:,-1] - z_r[:,:,-2])  # extrapolation slope Bottom
    W[ 0] = (3/8)*(R[ 0,:,:] + slopeT*(R[1,:,:] - R[0,:,:])) + (3/4)*R[ 1,:,:] + (1/8)*R[ 2,:,:]
    W[-1] = (3/8)*(R[-3,:,:] + slopeB*(R[1,:,:] - R[0,:,:])) + (3/4)*R[-2,:,:] + (1/8)*R[-1,:,:]

    return W




