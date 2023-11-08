#include "mod_cppkernels.h"


XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
void applyPointSources()
{
//# # Apply momentum transport point sources (like river runoff), if any.
//    # # -----------------------------------------------------------------------
//    #   IF (LuvSrc(ng)) THEN
//    #     DO is=1,Nsrc(ng)
//    #       i=SOURCES(ng)%Isrc(is)
//    #       j=SOURCES(ng)%Jsrc(is)
//    #       IF (((IstrR.le.i).and.(i.le.IendR)).and.                      &
//    #  &        ((JstrR.le.j).and.(j.le.JendR))) THEN
//    #         IF (INT(SOURCES(ng)%Dsrc(is)).eq.0) THEN
//    #           DO k=1,N(ng)
//    #             cff1=1.0_r8/(on_u(i,j)*                                 &
//    #  &                       0.5_r8*(z_w(i-1,j,k)-z_w(i-1,j,k-1)+       &
//    #  &                               z_w(i  ,j,k)-z_w(i  ,j,k-1)))
//    #             u(i,j,k,nnew)=SOURCES(ng)%Qsrc(is,k)*cff1
//    #           END DO
//    #         ELSE IF (INT(SOURCES(ng)%Dsrc(is)).eq.1) THEN
//    #           DO k=1,N(ng)
//    #             cff1=1.0_r8/(om_v(i,j)*                                 &
//    #  &                       0.5_r8*(z_w(i,j-1,k)-z_w(i,j-1,k-1)+       &
//    #  &                               z_w(i,j  ,k)-z_w(i,j  ,k-1)))
//    #             v(i,j,k,nnew)=SOURCES(ng)%Qsrc(is,k)*cff1
//    #           END DO
//    #         END IF
//    #       END IF
//    #     END DO
//    #   END IF
}











extern "C"  __global__
void createVertViscousOpMatrixU(const double *_DC,  const double *_BC, const double *_FC,
                                const double *_Akv, const double *_Hz, const double *_z_r,
                                const double *_u,   const double *_ru, const int N)
{
    //The Eq. is: Un+1 = Un - Δt*∂z(Akv*∂zUn) = [1 - ∂z(Δt*Akv*∂z)]Un = Ru
    //Integrating vertically between nodes k and k+1 we have (using Gauss central quadrature for the velocity):
    //
    //Hz*Un+1 = Hz*Un + Δt*Akv*(∂zU) = Hz*Un + Akv*(Un,k+1 - Un,k)/(Zk+1 - Zk)
    u  = u(nnew,:,:,:)
    ru = ru(nrhs,:,:,:)

    int K;
    STENCIL3D(DC, K);
    STENCIL3D(BC, K);
    STENCIL3D(FC, K);

    if (isUnode(i))
    {
        auto AK  = RtoU(Akv);
        auto HzU = RtoU(Hz)
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
            BC[k] = HzU(i,j,k) - FC(k,0,0) - FC(k-1,0,0);
        }
    }

}



extern "C"  __global__
void adjustBarotropicVelocity(const double *_Hz, const double *_u, const double *_v,
                              const double *_DU_avg1, const double *_DV_avg1, const int N)
// Replace INTERIOR POINTS? incorrect vertical mean velocity with more accurate barotropic component, ubar=DU_avg1/(D*om_u)
// and vbar=DV_avg1/(D*om_v).
//
//      _HzU, _HzV: Hz interpolated at U and V points.
//      _u, _v    : velocity components.
{
    v := v(i,j,1,nnew);

    int K = 0;   // This variable is used by all stencils to define the level they act upon.

    STENCIL3D(Hz,      K);
    STENCIL3D(u,       K);
    STENCIL3D(v,       K);
    STENCIL3D(om_u,    K);
    STENCIL3D(on_v,    K);
    STENCIL3D(DU_avg1, K);
    STENCIL3D(DV_avg1, K);

    auto HzU = RtoU(Hz);
    auto HzV = RtoV(Hz);

    double HU  = 0.0, HV  = 0.0;
    double HuU = 0.0, HvV = 0.0;

    for (K = 0; k < N; k++)
    {
        HU  += (HzU  ).Eval(0,0,0);
        HV  += (HzV  ).Eval(0,0,0);

        HuU += (HzU*u).Eval(0,0,0);
        HvV += (HzV*v).Eval(0,0,0);
    }

    auto uUpdate = (HuU - DU_avg1/HU)*on_u;
    auto vUpdate = (HvV - DV_avg1/HV)*om_v;


    // Update new solution.
    for (K=0; K<N; K++)
    {
        u -= uUpdate;
        v -= vUpdate;
    }
}