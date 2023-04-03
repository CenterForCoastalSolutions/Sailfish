program wvelminapp

      integer, parameter :: N = 11
      integer, parameter :: Istr = 1
      integer, parameter :: Jstr = 1
      integer, parameter :: Iend = 1000
      integer, parameter :: Jend = 1000
      integer, parameter :: r8 = 8

      integer, parameter :: LBi = Istr-1
      integer, parameter :: UBi = Iend+1
      integer, parameter :: LBj = Jstr-1
      integer, parameter :: UBj = Jend+1

      real(r8) :: pm(LBi:UBi,LBj:UBj)
      real(r8) :: pn(LBi:UBi,LBj:UBj)
      real(r8) :: z_r(LBi:UBi,LBj:UBj,N)
      real(r8) :: z_w(LBi:UBi,LBj:UBj,0:N)
      real(r8) :: u(LBi:UBi,LBj:UBj,N,2)
      real(r8) :: v(LBi:UBi,LBj:UBj,N,2)

      real(r8), dimension(Istr:Iend+1,Jstr:Jend+1,N) :: vert

      real(r8), dimension(Istr:Iend+1,Jstr:Jend+1) :: wrk

      integer :: i, j, k, NN
      integer :: Ninp, iter
      real(r8) :: res, t1, t2

      Ninp = 1


      DO NN = 1,2
        DO k=1,N
          DO j=LBj,UBj
            DO i=LBi,UBi
              u(i,j,k,NN) = i-j-k*2.0
              v(i,j,k,NN) = i-j+k*2.0
              pm(i,j) = i+1
              pn(i,j) = j+1
              z_r(i,j,k) = j-i-k
            END DO
          END DO
        END DO
      END DO

      print*, SHAPE(z_r)
          print*, SUM(z_r), SUM(u), SUM(v), SUM(pm), SUM(pn)


      call CPU_TIME(t1)
      !do iter2=1,10
          res = 0

          DO iter=1,1011

              DO k=1,N
                DO j=Jstr,Jend
                  DO i=Istr,Iend+1
                    wrk(i,j)=u(i,j,k,Ninp)*(z_r(i,j,k)-z_r(i-1,j,k))*(pm(i-1,j)+pm(i,j))

                  END DO
                  DO i=Istr,Iend
                    vert(i,j,k)=0.25_r8*(wrk(i,j)+wrk(i+1,j))
                    !print *, i, j, vert(i,j,k)

                  END DO
                END DO
                DO j=Jstr,Jend+1
                  DO i=Istr,Iend
                    wrk(i,j)=v(i,j,k,Ninp)*(z_r(i,j,k)-z_r(i,j-1,k))*(pn(i,j-1)+pn(i,j))

                  END DO
                END DO
                DO j=Jstr,Jend
                  DO i=Istr,Iend
                    vert(i,j,k)=vert(i,j,k)+0.25_r8*(wrk(i,j)+wrk(i,j+1))
                    !print *, i, j, 0.25_r8*(wrk(i,j)+wrk(i,j+1)), vert(i,j,k)
                  END DO
                END DO
              END DO
              !print *,shape(vert)
              res = res + sum(vert) ! (Istr+1:Iend,Jstr+1:Jend,:))
              !print *, res
          END DO
          print *,res
      !end do
      call CPU_TIME(t2)

      print *, "TIME: ", (t2-t1)/10

end program wvelminapp