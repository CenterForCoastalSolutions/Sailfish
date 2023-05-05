# This routine reads and reports physical model input parameters.
# =======================================================================

import mod_param
import mod_iounits
import mod_ncparam
import mod_netcdf
import mod_scalars
import mod_strings
import mod_boundary

import io_utils




    ifile    = 1             # multiple file counter
    ibcfile  = 1             # multiple BC file counter
    iclmfile = 1             # multiple CLM file counter
    itracer  = 0             # LBC tracer counter

    # Cdim = SIZE(Cval,1)
    # Clen = LEN(Cval(1))
    # Rdim = SIZE(Rval,1)
    # Nval=0



    # Read the input file into a dictionary.
    #--------------------------------------

    filenameInput = r'D:\projects\src\oceangpu\utility\test.in'

    input = Input(filenameInput)

    title    = input.getVal('TITLE')
    MyAppCPP = input.getVal('MyAppCPP')
    varname  = input.getVal('VARNAME')

    # CALL initialize_param    ! Continue allocating/initalizing
    # CALL allocate_scalars    ! variables since the application
    # CALL initialize_scalars  ! number of nested grids and
    # CALL allocate_ncparam    ! domain parameters are known
    # CALL initialize_ncparam

    BOUNDARY = Boundary(input)



    rho0    = getInputVal('RHO0',      dtype = float)
    bvf_bak = getInputVal('BVF_BAK',   dtype = float)


    time_ref = getInputVal('TIME_REF', dtype = float)
    CALL ref_clock (time_ref)




    # MIXING of MOMENTUM and  TRACERS

    nl_tnu2    = getInputVal('TNU2', count = NAT + NPT, dtype = float)
    nl_tnu4    = getInputVal('TNU4', count = NAT + NPT, dtype = float)


    nl_visc2   = getInputVal('VISC2',     dtype = float)
    nl_visc4   = getInputVal('VISC3',     dtype = float)
    LuvSponge  = getInputVal('LuvSponge', dtype = float)

    Akt_bak    = getInputVal('AKT_BAK',   count = NAT + NPT, dtype = float)
    Akt_limit  = getInputVal('AKT_LIMIT', count = NAT + NPT, dtype = float)

    Akv_bak    = getInputVal('AKV_BAK', dtype = float)
    Akv_limit  = getInputVal('AKV_LIMIT', dtype = float)
    ad_Akv_fac = getInputVal('ad_AKV_fac', dtype = float)
    Akk_bak    = getInputVal('AKK_BAK', dtype = float)
    Akp_bak    = getInputVal('AKP_BAK', dtype = float)
    tkenu2     = getInputVal('TKENU2', dtype = float)
    tkenu4     = getInputVal('TKENU4', dtype = float)

    # STATE EQUATION
    R0         = getInputVal('R0',    dtype = float, minVal = 0.0)
    T0         = getInputVal('T0',    dtype = float)
    S0         = getInputVal('S0',    dtype = float)
    Tcoef      = cp.abs(getInputVal('TCOEF', dtype = float))
    Scoef      = cp.abs(getInputVal('SCOEF', dtype = float))



    gamma2    = getInputVal('GAMMA2', dtype = float)



    LuvSrc    = getInputVal('LuvSrc', dtype = float)
    LwSrc     = getInputVal('LwSrc',  dtype = float)



    LsshCLM     = getInputVal('LsshCLM', dtype=float)
    Lm2CLM      = getInputVal('Lm2CLM', dtype=float)
    LnudgeM2CLM = getInputVal('LnudgeM2CLM', dtype=float)
    nl_tnu4 = getInputVal('TNU4', dtype=float)
    nl_tnu4 = getInputVal('TNU4', dtype=float)

    if (R0 < 100): R0 += 1000.0











            CASE ('NUDNAME')
              label='NUD - nudging coefficients'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Ngrids, Nfiles, NUD)
            CASE ('SSFNAME')
              label='SSF - Sources/Sinks forcing fields'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Ngrids, Nfiles, SSF)


            CASE ('NFFILES')
              Npts=load_i(Nval, Rval, Ngrids, nFfiles)
              DO ng=1,Ngrids
                IF (nFfiles(ng).le.0) THEN
                  IF (Master) WRITE (out,260) 'NFFILES', nFfiles(ng),   &
     &              'Must be equal or greater than one.'
                  exit_flag=4
                  RETURN
                END IF
              END DO
              max_Ffiles=MAXVAL(nFfiles)
              allocate ( FRC(max_Ffiles,Ngrids) )
              allocate ( FRCids(max_Ffiles,Ngrids) )
              allocate ( Ncount(max_Ffiles,Ngrids) )
              FRCids(1:max_Ffiles,1:Ngrids)=-1
              Ncount(1:max_Ffiles,1:Ngrids)=0
            CASE ('FRCNAME')
              label='FRC - forcing fields'
              DO ng=1,Ngrids
                IF (nFfiles(ng).lt.0) THEN
                  IF (Master) WRITE (out,290) 'nFfiles = ',             &
     &                                        nFfiles(ng),              &
     &              'KeyWord ''NFFILES'' unread or misssing from '//    &
     &              'input script ''roms.in''.'
                  exit_flag=4
                  RETURN
                END IF
              END DO
              Npts=load_s2d(Nval, Cval, Cdim, line, label, ifile,       &
     &                      igrid, Ngrids, nFfiles, Ncount, max_Ffiles,  FRC)







          END SELECT
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
      END DO
  10  IF (Master) WRITE (out,50) line
      exit_flag=4
      RETURN
  20  CLOSE (inp)






!
!  Check if both point sources methodologies are activated.  Only one
!  method is allowed for a particular grid.  Otherwise, the point
!  source will be applies twice.
!
      DO ng=1,Ngrids
        IF (LuvSrc(ng).and.LwSrc(ng)) THEN
          IF (Master) THEN
            WRITE (out,260) 'LuvSrc', LuvSrc(ng),                       &
     &            'Because LwSrc  is also ON; only one method is legal.'
            WRITE (out,260) 'LwSrc', LwSrc(ng),                         &
     &            'Because LuvSrc is also ON; only one method is legal.'
          END IF
          exit_flag=4
          RETURN
        END IF
      END DO

!
!  Make sure that both component switches are activated when processing
!  (Eastward,Northward) momentum components at RHO-points.
!
      DO ng=1,Ngrids
        IF (.not.Hout(idu2dE,ng).and.Hout(idv2dN,ng)) THEN
          Hout(idu2dE,ng)=.TRUE.
        END IF
        IF (.not.Hout(idv2dN,ng).and.Hout(idu2dE,ng)) THEN
          Hout(idv2dN,ng)=.TRUE.
        END IF

      END DO
!
!  Set various parameters.
!
      DO ng=1,Ngrids
!
!  Set switch to create history NetCDF file.
!
        IF ((nHIS(ng).gt.0).and.ANY(Hout(:,ng))) THEN
          LdefHIS(ng)=.TRUE.
        END IF
!
!  Set switch to create quicksave NetCDF file.
!
        IF ((nQCK(ng).gt.0).and.ANY(Qout(:,ng))) THEN
          LdefQCK(ng)=.TRUE.
        END IF


!  Set switch to process climatology file.
!

!
!  If appropriate, deactive outpur NetCDF files switches.
!
        IF (((nrrec(ng).eq.0).and.(nAVG(ng).gt.ntimes(ng))).or.         &
     &      (nAVG(ng).eq.0)) THEN
          LdefAVG(ng)=.FALSE.
        END IF
        IF (((nrrec(ng).eq.0).and.(nDIA(ng).gt.ntimes(ng))).or.         &
     &      (nDIA(ng).eq.0)) THEN
          LdefDIA(ng)=.FALSE.
        END IF
        IF (((nrrec(ng).eq.0).and.(nFLT(ng).gt.ntimes(ng))).or.         &
     &      (nFLT(ng).eq.0)) THEN
          LdefFLT(ng)=.FALSE.
        END IF
        IF (((nrrec(ng).eq.0).and.(nHIS(ng).gt.ntimes(ng))).or.         &
     &      (nHIS(ng).eq.0)) THEN
          LdefHIS(ng)=.FALSE.
        END IF
        IF (((nrrec(ng).eq.0).and.(nQCK(ng).gt.ntimes(ng))).or.         &
     &      (nQCK(ng).eq.0)) THEN
          LdefQCK(ng)=.FALSE.
        END IF
        IF (((nrrec(ng).eq.0).and.(nRST(ng).gt.ntimes(ng))).or.         &
     &      (nRST(ng).eq.0)) THEN
          LdefRST(ng)=.FALSE.
        END  IF
        IF (((nrrec(ng).eq.0).and.(nSTA(ng).gt.ntimes(ng))).or.         &
     &      (nSTA(ng).eq.0)) THEN
          LdefSTA(ng)=.FALSE.
        END IF
!
!  Determine switch to process boundary NetCDF file.
!
        ObcData(ng)=.FALSE.

      END DO
!
!  If multiple output files, edit derived type structure to store the
!  information about all multi-files.
!
      DO ng=1,Ngrids
        IF ((nHIS(ng).gt.0).and.(ndefHIS(ng).gt.0)) THEN
          OutFiles=ntimes(ng)/ndefHIS(ng)
          IF ((nHIS(ng).eq.ndefHIS(ng)).or.                             &
     &        (MOD(ntimes(ng),ndefHIS(ng)).ge.nHIS(ng))) THEN
            OutFiles=Outfiles+1
          END IF
          CALL edit_file_struct (ng, OutFiles, HIS)
        END IF
        IF ((nQCK(ng).gt.0).and.(ndefQCK(ng).gt.0)) THEN
          OutFiles=ntimes(ng)/ndefQCK(ng)
          IF ((nQCK(ng).eq.ndefQCK(ng)).or.                             &
     &        (MOD(ntimes(ng),ndefQCK(ng)).ge.nQCK(ng))) THEN
            OutFiles=Outfiles+1
          END IF
          CALL edit_file_struct (ng, OutFiles, QCK)
        END IF
      END DO


# !  Report input parameters.
# !-----------------------------------------------------------------------
# !
#       IF (Master.and.Lwrite) THEN
#         lstr=INDEX(my_fflags, 'free')-2
#         IF (lstr.le.0) lstr=LEN_TRIM(my_fflags)
#         WRITE (out,60) TRIM(title), TRIM(my_os), TRIM(my_cpu),          &
#      &                 TRIM(my_fort), TRIM(my_fc), my_fflags(1:lstr),   &
#      &                 TRIM(svn_url), TRIM(svn_rev),                    &
#      &                 TRIM(Rdir), TRIM(Hdir), TRIM(Hfile), TRIM(Adir)
# !
#         DO ng=1,Ngrids
# !
# !  Report grid size and domain decomposition.  Check for correct tile
# !  decomposition.
# !
#
#           WRITE (out,90) ng, Lm(ng), Mm(ng), N(ng), numthreads,         &
#      &                   NtileI(ng), NtileJ(ng)
#           IF (NtileI(ng)*NtileJ(ng).le.0) THEN
#             WRITE (out,100) ng
#             exit_flag=6
#             RETURN
#           END IF
#           IF (MOD(NtileI(ng)*NtileJ(ng),numthreads).ne.0) THEN
#             WRITE (out,100) ng
#             exit_flag=6
#             RETURN
#           END IF
#
# !
# !  Report physical parameters.





#      &          'Number of restart records to read from disk.'
#           WRITE (out,170) LcycleRST(ng), 'LcycleRST',                   &
#      &          'Switch to recycle time-records in restart file.'
#           WRITE (out,130) nRST(ng), 'nRST',                             &
#      &          'Number of timesteps between the writing of data',      &
#      &          'into restart fields.'
#           WRITE (out,130) ninfo(ng), 'ninfo',                           &
#      &          'Number of timesteps between print of information',     &
#      &          'to standard output.'
#           WRITE (out,170) ldefout(ng), 'ldefout',                       &
#      &          'Switch to create a new output NetCDF file(s).'
#           WRITE (out,130) nHIS(ng), 'nHIS',                             &
#      &          'Number of timesteps between the writing fields',       &
#      &          'into history file.'
#           IF (ndefHIS(ng).gt.0) THEN
#             WRITE (out,130) ndefHIS(ng), 'ndefHIS',                     &
#      &            'Number of timesteps between creation of new',        &
#      &            'history files.'
#           END IF
#           WRITE (out,130) nQCK(ng), 'nQCK',                             &
#      &          'Number of timesteps between the writing fields',       &
#      &          'into quicksave file.'
#           IF (ndefQCK(ng).gt.0) THEN
#             WRITE (out,130) ndefQCK(ng), 'ndefQCK',                     &
#      &            'Number of timesteps between creation of new',        &
#      &            'brief snpashots files.'
#           END IF
#
#           WRITE (out,200) rdrg(ng), 'rdrg',                             &
#      &          'Linear bottom drag coefficient (m/s).'
#           WRITE (out,200) rdrg2(ng), 'rdrg2',                           &
#      &          'Quadratic bottom drag coefficient.'
#           WRITE (out,200) Zob(ng), 'Zob',                               &
#      &          'Bottom roughness (m).'
#
#
#           WRITE (out,140) rho0, 'rho0',                                 &
#      &          'Mean density (kg/m3) for Boussinesq approximation.'
#           WRITE (out,140) dstart, 'dstart',                             &
#      &          'Time-stamp assigned to model initialization (days).'
#           WRITE (out,150) time_ref, 'time_ref',                         &
#      &          'Reference time for units attribute (yyyymmdd.dd)'
#
#
#           WRITE (out,210) Znudg(ng), 'Znudg',                           &
#      &          'Nudging/relaxation time scale (days)',                 &
#      &          'for free-surface.'
#           WRITE (out,210) M2nudg(ng), 'M2nudg',                         &
#      &          'Nudging/relaxation time scale (days)',                 &
#      &          'for 2D momentum.'
#
#           WRITE (out,210) obcfac(ng), 'obcfac',                         &
#      &          'Factor between passive and active',                    &
#      &          'open boundary conditions.'
#           WRITE (out,170) VolCons(1,ng), 'VolCons(1)',                  &
#      &          'NLM western  edge boundary volume conservation.'
#           WRITE (out,170) VolCons(2,ng), 'VolCons(2)',                  &
#      &          'NLM southern edge boundary volume conservation.'
#           WRITE (out,170) VolCons(3,ng), 'VolCons(3)',                  &
#      &          'NLM eastern  edge boundary volume conservation.'
#           WRITE (out,170) VolCons(4,ng), 'VolCons(4)',                  &
#      &          'NLM northern edge boundary volume conservation.'
#           WRITE (out,160) gamma2(ng), 'gamma2',                         &
#      &          'Slipperiness variable: free-slip (1.0) or ',           &
#      &          '                     no-slip (-1.0).'
#           IF (LuvSrc(ng)) THEN
#             WRITE (out,170) LuvSrc(ng), 'LuvSrc',                       &
#      &          'Turning ON  momentum point Sources/Sinks.'
#           ELSE
#             WRITE (out,170) LuvSrc(ng), 'LuvSrc',                       &
#      &          'Turning OFF momentum point Sources/Sinks.'
#           END IF
#           IF (LwSrc(ng)) THEN
#             WRITE (out,170) LwSrc(ng), 'LwSrc',                         &
#      &          'Turning ON  volume influx point Sources/Sinks.'
#           ELSE
#             WRITE (out,170) LwSrc(ng), 'LwSrc',                         &
#      &          'Turning OFF volume influx point Sources/Sinks.'
#           END IF
#           IF (LsshCLM(ng)) THEN
#             WRITE (out,170) LsshCLM(ng), 'LsshCLM',                     &
#      &          'Turning ON  processing of SSH climatology.'
#           ELSE
#             WRITE (out,170) LsshCLM(ng), 'LsshCLM',                     &
#      &          'Turning OFF processing of SSH climatology.'
#           END IF
#           IF (Lm2CLM(ng)) THEN
#             WRITE (out,170) Lm2CLM(ng), 'Lm2CLM',                       &
#      &          'Turning ON  processing of 2D momentum climatology.'
#           ELSE
#             WRITE (out,170) Lm2CLM(ng), 'Lm2CLM',                       &
#      &          'Turning OFF processing of 2D momentum climatology.'
#           END IF
#           IF (LnudgeM2CLM(ng)) THEN
#             WRITE (out,170) LnudgeM2CLM(ng), 'LnudgeM2CLM',             &
#      &          'Turning ON  nudging of 2D momentum climatology.'
#           ELSE
#             WRITE (out,170) LnudgeM2CLM(ng), 'LnudgeM2CLM',             &
#      &          'Turning OFF nudging of 2D momentum climatology.'
#           END IF
#           IF ((nHIS(ng).gt.0).and.ANY(Hout(:,ng))) THEN
#             WRITE (out,'(1x)')
#             IF (Hout(idFsur,ng)) WRITE (out,170) Hout(idFsur,ng),       &
#      &         'Hout(idFsur)',                                          &
#      &         'Write out free-surface.'
#             IF (Hout(idUbar,ng)) WRITE (out,170) Hout(idUbar,ng),       &
#      &         'Hout(idUbar)',                                          &
#      &         'Write out 2D U-momentum component.'
#             IF (Hout(idVbar,ng)) WRITE (out,170) Hout(idVbar,ng),       &
#      &         'Hout(idVbar)',                                          &
#      &         'Write out 2D V-momentum component.'
#             IF (Hout(idu2dE,ng)) WRITE (out,170) Hout(idu2dE,ng),       &
#      &         'Hout(idu2dE)',                                          &
#      &         'Write out 2D U-eastward  at RHO-points.'
#             IF (Hout(idv2dN,ng)) WRITE (out,170) Hout(idv2dN,ng),       &
#      &         'Hout(idv2dN)',                                          &
#      &         'Write out 2D V-northward at RHO-points.'
#             IF (Hout(idUsms,ng)) WRITE (out,170) Hout(idUsms,ng),       &
#      &         'Hout(idUsms)',                                          &
#      &         'Write out surface U-momentum stress.'
#             IF (Hout(idVsms,ng)) WRITE (out,170) Hout(idVsms,ng),       &
#      &         'Hout(idVsms)',                                          &
#      &         'Write out surface V-momentum stress.'
#             IF (Hout(idUbms,ng)) WRITE (out,170) Hout(idUbms,ng),       &
#      &         'Hout(idUbms)',                                          &
#      &         'Write out bottom U-momentum stress.'
#             IF (Hout(idVbms,ng)) WRITE (out,170) Hout(idVbms,ng),       &
#      &         'Hout(idVbms)',                                          &
#      &         'Write out bottom V-momentum stress.'
#           END IF
#
#           IF ((nQCK(ng).gt.0).and.ANY(Qout(:,ng))) THEN
#             WRITE (out,'(1x)')
#             IF (Qout(idFsur,ng)) WRITE (out,170) Qout(idFsur,ng),       &
#      &         'Qout(idFsur)',                                          &
#      &         'Write out free-surface.'
#             IF (Qout(idUbar,ng)) WRITE (out,170) Qout(idUbar,ng),       &
#      &         'Qout(idUbar)',                                          &
#      &         'Write out 2D U-momentum component.'
#             IF (Qout(idVbar,ng)) WRITE (out,170) Qout(idVbar,ng),       &
#      &         'Qout(idVbar)',                                          &
#      &         'Write out 2D V-momentum component.'
#             IF (Qout(idu2dE,ng)) WRITE (out,170) Qout(idu2dE,ng),       &
#      &         'Qout(idu2dE)',                                          &
#      &         'Write out 2D U-eastward  at RHO-points.'
#             IF (Qout(idv2dN,ng)) WRITE (out,170) Qout(idv2dN,ng),       &
#      &         'Qout(idv2dN)',                                          &
#      &         'Write out 2D V-northward at RHO-points.'
#             IF (Qout(idUsms,ng)) WRITE (out,170) Qout(idUsms,ng),       &
#      &         'Qout(idUsms)',                                          &
#      &         'Write out surface U-momentum stress.'
#             IF (Qout(idVsms,ng)) WRITE (out,170) Qout(idVsms,ng),       &
#      &         'Qout(idVsms)',                                          &
#      &         'Write out surface V-momentum stress.'
#             IF (Qout(idUbms,ng)) WRITE (out,170) Qout(idUbms,ng),       &
#      &         'Qout(idUbms)',                                          &
#      &         'Write out bottom U-momentum stress.'
#             IF (Qout(idVbms,ng)) WRITE (out,170) Qout(idVbms,ng),       &
#      &         'Qout(idVbms)',                                          &
#      &         'Write out bottom V-momentum stress.'
#
#           END IF
# !
# !-----------------------------------------------------------------------
# !  Report output/input files and check availability of input files.
# !-----------------------------------------------------------------------
# !
#           WRITE (out,220)
#
#
#           fname=varname
#           IF (.not.find_file(ng, fname, 'VARNAME')) GO TO 30
#           WRITE (out,230) 'ROMS I/O variables Metadata File:  ',        &
#      &                    TRIM(fname)
#           GO TO 40
#   30      IF (Master) THEN
#             IF (LEN_TRIM(fname).eq.0) THEN
#               WRITE (out,270) ng, 'Oops unassigned file name. '//       &
#      &                            'Check standard input script...'
#             ELSE
#               WRITE (out,270) ng, TRIM(fname)
#             END IF
#           END IF
#           exit_flag=4
#           IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   40      CONTINUE
#         END DO
#         IF (Nuser.gt.0) THEN
#           WRITE (out,230) '          Input/Output USER File:  ',        &
#      &                    TRIM(USRname)
#         END IF
# !
# !-----------------------------------------------------------------------
# !  Report generic USER parameters.
# !-----------------------------------------------------------------------
# !
#         IF (Nuser.gt.0) THEN
#           WRITE (out,240)
#           DO i=1,Nuser
#             WRITE (out,250) user(i), i, i
#           END DO
#         END IF
#       END IF
#
# !
#
#   50  FORMAT (/,' READ_PhyPar - Error while processing line: ',/,a)
#   60  FORMAT (/,1x,a,/,                                                 &
#      &        /,1x,'Operating system  : ',a,                            &
#      &        /,1x,'CPU/hardware      : ',a,                            &
#      &        /,1x,'Compiler system   : ',a,                            &
#      &        /,1x,'Compiler command  : ',a,                            &
#      &        /,1x,'Compiler flags    : ',a,                            &
#      &        /,1x,'SVN Root URL     : ',a,                             &
#      &        /,1x,'SVN Revision     : ',a,/,                           &
#      &        /,1x,'Local Root       : ',a,                             &
#      &        /,1x,'Header Dir       : ',a,                             &
#      &        /,1x,'Header file      : ',a,                             &
#      &        /,1x,'Analytical Dir   : ',a)
#   70  FORMAT (/,' Resolution, Grid ',i2.2,': ',i0,'x',i0,'x',i0,        &
#      &        ',',2x,'Parallel Nodes: ',i0,',',2x,'Tiling: ',i0,        &
#      &        'x',i0)
#   80  FORMAT (/,' ROMS/TOMS: Wrong choice of grid ',i2.2,1x,            &
#      &        'partition or number of parallel nodes.',                 &
#      &        /,12x,a,1x,i0,/,12x,                                      &
#      &        'must be equal to the number of parallel processes = ',   &
#      &        i0,/,12x,'Change -np value to mpirun or',                 &
#      &        /,12x,'change domain partition in input script.')
#   90  FORMAT (/,' Resolution, Grid ',i2.2,': ',i0,'x',i0,'x',i0,        &
#      &        ',',2x,'Parallel Threads: ',i0,',',2x,'Tiling: ',i0,      &
#      &        'x',i0)
#  100  FORMAT (/,' ROMS/TOMS: Wrong choice of grid ',i2.2,1x,            &
#      &        'partition or number of parallel threads.',               &
#      &        /,12x,'NtileI*NtileJ must be a positive multiple of the', &
#      &        ' number of threads.',                                    &
#      &        /,12x,'Change number of threads (environment variable) ', &
#      &        'or',/,12x,'change domain partition in input script.')
#  110  FORMAT (/,/,' Physical Parameters, Grid: ',i2.2,                  &
#      &        /,  ' =============================',/)
#  120  FORMAT (1x,i10,2x,a,t32,a)
#  130  FORMAT (1x,i10,2x,a,t32,a,/,t34,a)
#  140  FORMAT (f11.3,2x,a,t32,a)
#  150  FORMAT (f11.2,2x,a,t32,a)
#  160  FORMAT (f11.3,2x,a,t32,a,/,t34,a)
#  170  FORMAT (10x,l1,2x,a,t32,a)
#  180  FORMAT (10x,l1,2x,a,t32,a,i2.2,':',1x,a)
#  185  FORMAT (10x,l1,2x,a,'(',i2.2,')',t32,a,i2.2,':',1x,a)
#  190  FORMAT (1p,e11.4,2x,a,'(',i2.2,')',t32,a,/,t34,a,i2.2,':',1x,a)
#  195  FORMAT (1p,e11.4,2x,a,t32,a,i2.2,':',1x,a)
#  200  FORMAT (1p,e11.4,2x,a,t32,a)
#  210  FORMAT (1p,e11.4,2x,a,t32,a,/,t34,a)
#  220  FORMAT (/,' Output/Input Files:',/)
#  230  FORMAT (2x,a,a)
#  240  FORMAT (/,' Generic User Parameters:',/)
#  250  FORMAT (1p,e11.4,2x,'user(',i2.2,')',t32,                         &
#      &        'User parameter ',i2.2,'.')
#  260  FORMAT (/,' READ_PHYPAR - Invalid input parameter, ',a,           &
#      &        i4,/,15x,a)
#  265  FORMAT (/,' READ_PHYPAR - Invalid input parameter, ',a,           &
#      &        1p,e11.4,/,15x,a)
#  270  FORMAT (/,' READ_PHYPAR - Grid ',i2.2,                            &
#      &        ', could not find input file:  ',a)
#  280  FORMAT (/,' READ_PHYPAR - Variable index not yet loaded, ', a)
#  290  FORMAT (/,' READ_PHYPAR - Invalid dimension parameter, ',a,i0,    &
#      &        /,15x,a)
#  300  FORMAT (/,' READ_PHYPAR - Invalid dimension parameter, ',a,'(',   &
#      &        i2.2,')',/,15x,a)
#  310  FORMAT (2x,a,i2.2,a,a)
#  320  FORMAT (/,' READ_PHYPAR - Could not find input parameter: ', a,   &
#      &        /,15x,'in ROMS standard input script.',/,15x,a)
#  330  FORMAT (/,' READ_PHYPAR - Invalid input parameter, ',a,i4,/,15x,a)
#  340  FORMAT (/,' READ_PHYPAR - Inconsistent time-stepping period:',    &
#      &        /,15x,'Grid ',i2.2,':',f14.1,' (sec)',2x,f14.2,' (days)', &
#      &        /,15x,'Grid ',i2.2,':',f14.1,' (sec)',2x,f14.2,' (days)', &
#      &        /,15x,'Adjust standard input parameter NTIMES in ',       &
#      &              '''roms.in''.'/)
#  350  FORMAT (/,' READ_PHYPAR - Invalid input parameter, ',a,i0,        &
#      &        ', for grid ',i2.2,/,15x,a,i0,', ',a,i0,/,15x,a,/,15x,a)
#       RETURN
#       END SUBROUTINE read_PhyPar
