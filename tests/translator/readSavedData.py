import os
import cupy as np


def readVar(file, isInput, dictIn):
    vardef = file.readline().strip()
    if vardef == 'endofdatafile':
        return False
    if (vardef != 'vardef'):
        print('ERROR: "vardef" expected, found: "%s"' % vardef)
    var = file.readline().strip()
    origVar = var
    ioFlag = int(file.readline().strip())
    strDataType = file.readline().strip()


    # if 'DOMAIN'

    if 'array' in strDataType:
        dims = int(file.readline().strip())
        strShape1 = np.array(file.readline().strip().split()).astype(int)
        strShape2 = np.array(file.readline().strip().split()).astype(int)


        print(strShape1)
        if ((strShape1!=1).any()):
            print(strShape1)

        if (strDataType == 'float64_array'):
            dataType = np.float64
        elif (strDataType == 'float32_array'):
            dataType = np.float32
        elif (strDataType == 'int_array'):
            dataType = np.int
        elif (strDataType == 'logical_array'):
            dataType = np.bool

        varLong = 'globals()["%s"]' % var
        try:
            while var.find('(')>-1:
                idx = var.find('(')
                var = var[:idx] + '[' + var[idx+1:]
                idx = var.find(')')
                var = var[:idx] + ']' + var[idx+1:]
            while var.find('%') > -1:
                idx = var.find('%')
                var = var[:idx] + '.' + var[idx+1:]

            class Object(object):
                pass

            idx1 = var.find('[')
            idx2 = var.find(']')
            if idx1>-1 and idx2>-2:
                varLong = 'globals()["%s"][globals()["%s"]]' % (var[:idx1], var[idx1 + 1:idx2])
                if var[:idx1] not in globals():
                    exec('globals()["%s"] = [Object(), Object()]' % var[:idx1])
                    # exec('globals()["%s"] = 0' % (var[idx1 + 1:idx2]))
                    # exec('%s = 0' % varLong)
                varLong = '%s.%s' % (varLong ,var[idx2 + 2:])
                # exec('tempVar = globals()["%s"][%s]' % (var[:idx1],var[idx1+1:idx2]))



            try:
                exec('%s = np.zeros(np.flip(strShape2), dtype = dataType)' % (varLong))
            except:
                exec('tempVar.%s = %s' % (var[idx2 + 2:]))

        except:
            print (var, strDataType)

        print(var)

        line = np.array(file.readline().strip().split()).astype(dataType)
        numValsPerLine = line.size
        exec('%s.flat[:%i] = line' % (varLong, numValsPerLine))  # (var, idx, idx + line.size))

        count = strShape2.prod()-numValsPerLine
        fullRows = int(count / numValsPerLine)
        exec('%s.flat[:] = np.loadtxt(file, max_rows=%i)' % (varLong, fullRows))

        # idx = 0
        # while (count>0):
        if (count> fullRows*numValsPerLine):
            try:
                line = np.array(file.readline().strip().split()).astype(dataType)
                exec('%s.flat[%i:%i] = line' %  (varLong, fullRows*3, count))   #(var, idx, idx + line.size))
            except:
                print(222222)
            # count -= line.size
            # idx += line.size

            exec('dictIn["%s"] = %s' % (origVar,varLong))

    else:
        if (strDataType == 'float64'):
            dataType = np.float64
        elif (strDataType == 'float32'):
            dataType = np.float32
        elif (strDataType == 'int'):
            dataType = int
        elif (strDataType == 'logical'):
            dataType = np.bool

        strVal = file.readline().strip()
        val = dataType(strVal)

        try:
            exec('globals()["%s"] = val' % var)
            exec('dictIn["%s"] = val' % origVar)
        except:
            exec('globals()["temp_%s"] = val' % var)
            exec('dictIn["temp_%s"] = val' % origVar)






        # shape = int()
    # print(tmpArray)

    return True

inputFilename = 'step2d          10 .dat'
dictIn = {}
with open(inputFilename, 'r') as file:
    while readVar(file, True, dictIn):
        pass


pass