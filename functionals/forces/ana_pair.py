import mod_param
import mod_forces
import mod_ncparam
import mod_scalars


def ana_pair(model):
# Sets surface air pressure (mb) using an analytical expression.                                                         !

# Set analytical surface air pressure (mb).
# (1 mb = 100 Pa = 1 hPa, 1 bar = 1.0e+5 N/m2 = 1.0e+5 dynes/cm2).



#if defined BENCHMARK
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          Pair(i,j)=1025.0_r8
        END DO
      END DO
#elif defined BL_TEST
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          Pair(i,j)=1013.48_r8
        END DO
      END DO
#else
      ana_pair.h: no values provided for Pair.
#endif

