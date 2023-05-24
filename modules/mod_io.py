import cupy as cp
from misc import *




class T_IO:

    def __init__(self, input):


        self.nRST    = input.getVal('NRST',   minVal = -1, dtype = int)
        self.nSTA    = input.getVal('NSTA',   minVal = -1, dtype = int)
        self.nFLT    = input.getVal('NFLT',   minVal = -1, dtype = int)
        self.ninfo   = input.getVal('NINFO',  minVal = 0,  dtype = int)

        self.ldefout = input.getVal('LDEFOUT', minVal = 0, dtype = bool)


        self.nHIS    = input.getVal('NHIS',      minVal = -1, dtype = int)
        self.ndefHIS = input.getVal('NDEFHIS',   minVal = -1, dtype = int)
        self.nQCK    = input.getVal('NQCK',      minVal = -1, dtype = int)
        self.ndefQCK = input.getVal('NDEFQCK',   minVal = -1, dtype = int)




        # HISTORY output file.
        # =============================================================================
        idVars = ('idUvel', 'idVvel', 'idWvel', 'idOvel', 'idUbar', 'idVbar', 'idFsur', 'idu2dE', 'idv2dN', 'idu3dE', 'idv3dN', 'idTvar',
                  'idpthR', 'idpthU', 'idpthV', 'idpthW', 'idUsms', 'idVsms', 'idUbms', 'idVbms', 'idHsbl', 'idHbbl', 'idMtke', 'idMtls')


        self.Hout = {}
        for idVar in idVars:
            try:
                self.Hout[idVar] = input.getVal('Hout(%s)' % idVar, dtype = bool)
            except:
                self.Hout[idVar] = False




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
            try:
                self.Qout[idVar] = input.getVal('Qout(%s)' % idVar, dtype=bool)
            except:
                self.Qout[idVar] = False





        # data for NETCDF
        # =============================================================================
        self.shuffle = input.getVal('NC_SHUFFLE', dtype = int)
        self.shuffle = input.getVal('NC_DEFLATE', dtype = int)
        self.shuffle = input.getVal('NC_DLEVEL',  dtype = int)



        # Make sure that both component switches are activated when processing (Eastward, Northward) momentum components at RHO-points.
        try:
            if self.Hout['idu2dE'] or self.Hout['idv2dN']:
                self.Hout['idu2dE'] = True
                self.Hout['idv2dN'] = True
        except:
            pass


        # Set switches to create NetCDF files.
        self.LdefHIS = (self.nHIS > 0) and cp.any(list(self.Hout.values()))
        self.LdefQCK = (self.nQCK > 0) and cp.any(list(self.Hout.values()))




