
import misc




class T_IO:

    def __init__(self):


        nRST  = getInputVal('NRST',   minVal = -1, dtype = int)
        nSTA  = getInputVal('NSTA',   minVal = -1, dtype = int)
        nFLT  = getInputVal('NFLT',   minVal = -1, dtype = int)
        ninfo = getInputVal('NINFO',  minVal = 0,  dtype = int)

        ldefout = getInputVal('LDEFOUT', minVal = 0, dtype = bool)


        nHIS    = getInputVal('NHIS',      minVal = -1, dtype = int)
        ndefHIS = getInputVal('NDEFHIS',   minVal = -1, dtype = int)
        nQCK    = getInputVal('NQCK',      minVal = -1, dtype = int)
        ndefQCK = getInputVal('NDEFQCK',   minVal = -1, dtype = int)




        # HISTORY output file.
        # =============================================================================
        idVars = ('idUvel', 'idVvel', 'idWvel', 'idOvel', 'idUbar', 'idVbar', 'idFsur', 'idu2dE', 'idv2dN', 'idu3dE', 'idv3dN', 'idTvar',
                  'idpthR', 'idpthU', 'idpthV', 'idpthW', 'idUsms', 'idVsms', 'idUbms', 'idVbms', 'idHsbl', 'idHbbl', 'idMtke', 'idMtls')


        self.Hout = {}
        for idVar in idVars:
            self.Hout[id] = getInputVal('Hout(%s)' % idVar, dtype = bool)




        # QUICK output file.
        # =============================================================================
        idVars = ('idUvel', 'idVvel', 'idWvel', 'idOvel', 'idUbar', 'idVbar', 'idFsur', 'idu2dE', 'idv2dN', 'idu3dE', 'idv3dN', 'idTvar',
                  'idUsur', 'idVsur', 'idUsuE', 'idVsuN', 'idsurT',
                  'idpthR', 'idpthU', 'idpthV', 'idpthW', 'idUsms', 'idVsms', 'idUbms', 'idVbms', 'idW2xx', 'idW2xy',
                  'idW2yy', 'idU2rs', 'idV2rs', 'idU2Sd',
                  'idV2Sd', 'idW3xx', 'idW3xy', 'idW3yy', 'idW3zx', 'idW3zy', 'idU3rs', 'idV3rs', 'idU3Sd', 'idV3Sd',
                  'idWamp', 'idWlen', 'idWdir', 'idWdip',
                  'idLhea', 'idShea', 'idLrad', 'idSrad', 'idEmPf', 'idevap', 'idrain', 'idDano', 'idVvis', 'idTdif',
                  'idHsbl', 'idHbbl', 'idMtke', 'idMtls'
                  )


        self.Qout = {}
        for idVar in idVars:
            self.Qout[id] = getInputVal('Qout(%s)' % idVar, dtype=bool)





        # data for NETCDF
        # =============================================================================
        self.shuffle = getInputVal('NC_SHUFFLE', dtype = 'int')
        self.shuffle = getInputVal('NC_DEFLATE', dtype = 'int')
        self.shuffle = getInputVal('NC_DLEVEL',  dtype = 'int')

