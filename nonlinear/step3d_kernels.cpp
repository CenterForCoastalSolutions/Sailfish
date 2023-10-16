#include "mod_cppkernels.h"

extern "C"  __global__
void createVertViscousOpMatrixU(const double *_DC,  const double *_BC, const double *_FC,
                                const double *_Akv, const double *_Hz, const double *_z_r,
                                const double *_u,   const double *_ru, const int N)
{
    //The Eq. is: Un+1 = Un - Δt*∂z(Akv*∂zU) = [1 - ∂z(Akv*∂z)]U = Ru
    //Integrating vertically between nodes k and k+1 we have (using Gauss central quadrature for the velocity):
    //
    //Hz*Un+1 = Hz*Un + Δt*Akv*(∂zU) = Hz*Un + Akv*(Un,k+1 - Un,k)/(Zk+1 - Zk)
    u  = u(nnew,:,:,:)
    ru = ru(nrhs,:,:,:)

    STENCIL3D(DC);
    STENCIL3D(BC);
    STENCIL3D(FC);

    if (isUnode(i))
    {
        auto AK  = RtoU(Akv);
        auto Hzk = RtoU(Hz)
        auto zU  = RtoU(z_r);

        DC = cff*RtoU(pm)*RtoU(pn);
        for (int k=1; k<N; k++)
        {
            u += DC*ru(k,0,0);

            double Δz = zU(k+1,0,0) - zU(k,0,0);

            FC(k) = (-lambda*Δt)*AK(i,k)/Δz;
        }

        FC(0,0,0) = 0.0;
        FC(N,0,0) = 0.0;

        for (int k=1; k<=N; k++)
        {
            // RHS
            DC[k] = u(i,j,k,nnew);

            // Diagonal includes Hz and two of the components of (Un,k+1 - Un,k)/(Zk+1 - Zk) coming from the elements up and down
            BC[k] = Hzk(i,j,k) - FC(k,0,0) - FC(k-1,0,0);
        }
    }

}




// Replace INTERIOR POINTS incorrect vertical mean with more accurate barotropic component, vbar=DV_avg1/(D*om_v).
// Recall that, D=CF(:,0).

       DO i=Istr,Iend
         CF(i,0)=Hzk(i,1)
         DC(i,0)=v(i,j,1,nnew)*Hzk(i,1)

       END DO
       DO k=2,N(ng)
         DO i=Istr,Iend
           CF(i,0)=CF(i,0) + Hzk(i,k)
           DC(i,0)=DC(i,0) + v(i,j,k,nnew)*Hzk(i,k)

         END DO
       END DO

       DO i=Istr,Iend
         DC(i,0) = (DC(i,0)*om_v(i,j) - DV_avg1(i,j))/(CF(i,0)*om_v(i,j))    ! recursive

       END DO

 # Couple and update new solution.

       DO k=1,N(ng)
         DO i=Istr,Iend
           v(i,j,k,nnew) = v(i,j,k,nnew) - DC(i,0)


         END DO
       END DO


     END IF
   END DO
