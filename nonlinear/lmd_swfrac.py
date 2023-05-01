#if defined NONLINEAR && (defined LMD_SKPP || defined SOLAR_SOURCE) &&  defined SOLVE3D

def lmd_swfrac(ng, Zscale, Z, swdk)
# This subroutine computes the  fraction  of  solar shortwave flux    !
# penetrating to specified depth (times Zscale) due to exponential    !
# decay in Jerlov water type.                                         !
#                                                                     !
# On Input:                                                           !
#                                                                     !
#    Zscale   scale factor to apply to depth array.                   !
#    Z        vertical height (meters, negative) for                  !
#             desired solar short-wave fraction.                      !
#                                                                     !
# On Output:                                                          !
#                                                                     !
#    swdk     shortwave (radiation) fractional decay.                 !
#                                                                     !
# Reference:                                                          !
#                                                                     !
#   Use Paulson and Simpson (1977) two wavelength bands solar
#   absorption model.
#   Paulson, C.A., and J.J. Simpson, 1977: Irradiance meassurements     !
#   in the upper ocean, J. Phys. Oceanogr., 7, 952-956.              !


    mixing = MIXING[ng]

     # "Jerlov water type" of types I, IA, IB, II and III (see paper)
    Jindex = int(mixing.Jwtype(i,j));

'''
    // lmd_mu1,2: Reciprocals of the absorption coefficient for solar wavelength bands 1,2 as a function of the Jerlov water type.
    const double fac1 = Zscale/lmd_mu1[Jindex];
    const double fac2 = Zscale/lmd_mu2[Jindex];

    // Weight (R in paper)
    const double w = lmd_r1(Jindex);

    swdk(0,0) = w*exp(Z(0,0)*fac1) + (1-w)*exp(Z(0,0)*fac2);
'''