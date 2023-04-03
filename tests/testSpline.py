import numpy as np
import matplotlib.pyplot as plt

N = 20

X = np.arange(N)/2
Y = np.sin(X)
Y[:] = 0.0
Y[10] = 1
H = np.ones(N)
CF = np.zeros(N)
dY = np.zeros(N)

for j,(x,y) in enumerate(zip(X[1:-1],Y[1:-1])):

    i = j + 1
    print(i)
    cff = 1.0 / (2.0 * H[i+1] + H[i] * (2.0 - CF[i-1]))
    CF[i] = cff * H[i+1]

    dY[i] = cff * (3.0 * 2*(Y[i+1] - Y[i]) - H[i] * dY[i-1])

for j,(x,y) in enumerate(zip(reversed(X[1:-1]),reversed(Y[1:-1]))):
    i = N - j - 2
    print(i)
    dY[i] -= CF[i] * dY[i+1]

pass


# dU(i,k)=cff*(3.0_r8*(u(i  ,j,k+1,nstp)-u(i,  j,k,nstp)+     &
#      &                           u(i+1,j,k+1,nstp)-u(i+1,j,k,nstp))-    &
#      &                   Hz(i,j,k)*dU(i,k-1))
