

kernel = cp.ElementwiseKernel(...)


CF = cp.zeros(Hz.shape, dtype=cp.float64)
kernel(Hz[0,:,:], dU[0,:,:], dV[0,:,:], u[nstp,0,:,:], v[nstp,0,:,:], CF[0,:,:], shear2[0,:,:])


preamble = '''
// PREAMBLE?
//-----------
template<class T>
class Stencil
{
    T *p;
    size_t strideI, strideJ, strideK;

public:
    Stencil(T *_p): p(_p)
    {
      INITIALIZE strideI, strideJ, strideK
    }

    T &operator()(size_t const k, size_t const j, size_t const i) const
    {
        return *(p + i*strideI + j*strideJ + k*strideK);
    }

    T operator=(T val) const
    {
        return (*p = val);
    }
}
'''

// KERNEL
//-----------
//void kernel(double *_Hz, double *_dU, double *_dV, double *_u, double *_v, double *_CF, double *_shear2)

    Stencil Hz(_Hz), dU(_dU), dV(_dV), u(_u), v(_v), CF(_CF), shear2(_shear2);

    dU(0,0,0) = 0.0;
    dV(0,0,0) = 0.0;
    dU(N,0,0) = 0.0;
    dU(N,0,0) = 0.0;

    #pragma unroll 1
    for (int k=1; k<N; k++)
    {
        const double cff = 1.0/(2.0*Hz(k+1,0,0) + Hz(k,0,0)*(2.0 - CF(k-1,0,0)));

        CF = cff*Hz(1,0,0);

        dU = cff*(3.0*(u(k+1,0,0) - u(k,0,0) + u(k+1,0,1) - u(k,0,1))) - Hz(k,0,0)*dU(k-1,0,0)
        dV = cff*(3.0*(v(k+1,0,0) - v(k,0,0) + v(k+1,1,0) - v(k,1,0))) - Hz(k,0,0)*dV(k-1,0,0)
    }

    #pragma unroll 1
    for (int k=N-1; k>=1; k--)
    {
        dU += dU - CF(k,0,0)*dU(k+1,0,0)   // Why not included in the previous line?
        dV += dV - CF(k,0,0)*dV(k+1,0,0)

        shear2 = dU(k,0,0)*dU(k,0,0) + dV(k,0,0)*dV(k,0,0)
    }




//ORIGINAL FORTRAN CODE


      DO j=Jstrm1,Jendp1
        DO i=Istrm1,Iendp1
          CF(i,0)=0.0_r8
          dU(i,0)=0.0_r8
          dV(i,0)=0.0_r8
        END DO

        DO k=1,N(ng)-1
          DO i=Istrm1,Iendp1
            cff=1.0_r8/(2.0_r8*Hz(i,j,k+1)+ Hz(i,j,k)*(2.0_r8-CF(i,k-1)))
            CF(i,k)=cff*Hz(i,j,k+1)
            dU(i,k)=cff*(3.0_r8*(u(i  ,j,k+1,nstp)-u(i,  j,k,nstp) + u(i+1,j,k+1,nstp)-u(i+1,j,k,nstp)) - Hz(i,j,k)*dU(i,k-1))
            dV(i,k)=cff*(3.0_r8*(v(i,j  ,k+1,nstp)-v(i,j  ,k,nstp) + v(i,j+1,k+1,nstp)-v(i,j+1,k,nstp)) - Hz(i,j,k)*dV(i,k-1))
          END DO
        END DO


        DO i=Istrm1,Iendp1
          dU(i,N(ng))=0.0_r8
          dV(i,N(ng))=0.0_r8
        END DO

        DO k=N(ng)-1,1,-1
         DO i=Istrm1,Iendp1
            dU(i,k)=dU(i,k)-CF(i,k)*dU(i,k+1)
            dV(i,k)=dV(i,k)-CF(i,k)*dV(i,k+1)
          END DO
        END DO

        DO k=1,N(ng)-1
          DO i=Istrm1,Iendp1
            shear2(i,j,k)=dU(i,k)*dU(i,k)+dV(i,k)*dV(i,k)
          END DO
        END DO

      END DO



