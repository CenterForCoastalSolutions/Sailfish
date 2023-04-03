/*
 * ** svn $Id: inlet_test.h 838 2008-11-17 04:22:18Z jcwarner $
 * *******************************************************************************
 * ** Copyright (c) 2002-2008 The ROMS/TOMS Group                               **
 * **   Licensed under a MIT/X style license                                    **
 * **   See License_ROMS.txt                                                    **
 * *******************************************************************************
 * **
 * ** Options for Inlet Test Case, waves-ocean (SWAN/ROMS) two-way coupling.
 * **
 * ** Application flag:   INLET_TEST
 * ** Input script:       ocean_inlet_test.in
 * **                     coupling_inlet_test.in
 * **                     sediment_inlet_test.in
 * */

#define ROMS_MODEL
#undef SWAN_MODEL

#ifdef SWAN_MODEL
#define MCT_LIB
#define WEC_VF
#define WDISS_WAVEMOD
#define UV_KIRBY
#endif
#undef NESTING

#define UV_VIS4
#define UV_C4VADVECTION
#define UV_U3HADVECTION
#define VAR_RHO_2D

#define DIFF_GRID
#define TS_DIF4
#define MIX_S_UV
#define MIX_S_TS
#define NONLIN_EOS
#define NONLINEAR
#define MASKING
#define UV_ADV
#define WET_DRY
#define UV_COR
#define CURVGRID
#define SPHERICAL
#define SALINITY
#define DJ_GRADPS
#define SSH_TIDES
#define UV_TIDES
#define RAMP_TIDES
#define PERFECT_RESTART
#define ADD_M2OBC
#define ADD_FSOBC
#define ANA_FSOBC
#define ANA_M2OBC
#define ANA_M3OBC
#undef ANA_TOBC

#define SOLVE3D
#define SPLINES_VVISC
#define SPLINES_VDIFF
#undef ANA_INITIAL

#define DEFLATE
#define DOUBLE_PRECISION

#undef UV_QDRAG
#define UV_LOGDRAG
#define GLS_MIXING
#ifdef GLS_MIXING
#define KANTHA_CLAYSON
#define N2S2_HORAVG
#define CHARNOK
#define CRAIG_BANNER
#endif

#define BULK_FLUXES
#define ATM_PRESS
#ifdef BULK_FLUXES
#ifdef SWAN_MODEL
#undef COARE_OOST
#define COARE_TAYLOR_YELLAND
#endif
#define EMINUSP
#define SOLAR_SOURCE
#undef ANA_SMFLUX
# else
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_SMFLUX
# endif

#define ANA_SPFLUX
#define ANA_BTFLUX
#define ANA_BSFLUX
#define ANA_BPFLUX
#undef ANA_SRFLUX
#define AKLIMIT

#undef SEDIMENT
#ifdef SEDIMENT
#define SUSPLOAD
#define BEDLOAD_SOULSBY
#define SED_SLUMP
#undef NONCOHESIVE_BED1
#define ANA_SEDIMENT
#define SED_MORPH
#endif



#define VERIFICATION
#define ENKF_RESTART