# This module sets vertical boundary conditons for momentum and tracers.                                                            !



def set_vbc():
"""Set kinematic barotropic bottom momentum stress (m2/s2)."""

#  if defined UV_LDRAG

# Set linear bottom stress.

    bustr[:,:] = 0.5*DηRtoU(rdrag)*ubar[krhs,:,:]
    bvstr[:,:] = 0.5*DxRtoV(rdrag)*vbar[krhs,:,:]


# #  elif defined UV_QDRAG
# !
# !  Set quadratic bottom stress.
# !
#       DO j=Jstr,Jend
#         DO i=IstrU,Iend
#           cff1=0.25_r8*(vbar(i  ,j  ,krhs)+                             &
#      &                  vbar(i  ,j+1,krhs)+                             &
#      &                  vbar(i-1,j  ,krhs)+                             &
#      &                  vbar(i-1,j+1,krhs))
#           cff2=SQRT(ubar(i,j,krhs)*ubar(i,j,krhs)+cff1*cff1)
#           bustr(i,j)=0.5_r8*(rdrag2(i-1,j)+rdrag2(i,j))*                &
#      &               ubar(i,j,krhs)*cff2
#         END DO
#       END DO
#       DO j=JstrV,Jend
#         DO i=Istr,Iend
#           cff1=0.25_r8*(ubar(i  ,j  ,krhs)+                             &
#      &                  ubar(i+1,j  ,krhs)+                             &
#      &                  ubar(i  ,j-1,krhs)+                             &
#      &                  ubar(i+1,j-1,krhs))
#           cff2=SQRT(cff1*cff1+vbar(i,j,krhs)*vbar(i,j,krhs))
#           bvstr(i,j)=0.5_r8*(rdrag2(i,j-1)+rdrag2(i,j))*                &
#      &               vbar(i,j,krhs)*cff2
#         END DO
#       END DO
# #  endif



    # Apply boundary conditions.
    # TODO: Is this necessary here?
    barotropicVelocityBC(bustr, bvstr)



