from os.path import join, exists
from misc import *
# This class contains general information about the input variables. This information is stored in the text dile input_vars.dat

class VarInfo:
    def __init__(self, varName, desc, obsolete):
        self.varName         = varName
        self.desc            = desc
        self.obsolete        = obsolete
        self.additionalInfo  = 'Not used.'    # This information is updated when (and if) the variable is parsed using getVal()


# This class contains information about the variables as given by the text file 'input_var.dat'. All variables used in the input file MUST be included
# in 'input_var.dat'.
# The intention here is to simplify the code and make it more flexible by removing this information from the code itself. Think for example in how
# ROMS code prints out the values read from the input file.

class VarInfoList:
    def __init__(self, path = '.'):

        self.obsolete = []
        varInfoFilename = join(path, 'input_params_info.dat')
        self.info = {}
        lineNum = 0
        for line in open(varInfoFilename):

            lineNum += 1

            line = line.strip()
            if (len(line) == 0) or (line[0] == '!'):
                continue

            # At this point we are in a line with information. These lines are of the form: VAR         Description@more description@more lines...
            idxSpace = line.find(' ')
            var  = line[:idxSpace].strip()
            desc = line[idxSpace:].strip()

            # if a variable name starts with a * it means it is obsolete (we keep this info to produce meaningful error messages).
            obsolete = False
            if var[0] == '*':
                var = var[1:]
                obsolete = True
                self.obsolete += [var]

            # For the model, we only keep the short description (up the the first "@")
            idxSpace = line.find('@')
            if idxSpace>0:
                desc = line[:idxSpace].strip()


            # Stores the information in the dictionary.
            if var in self.info:
                msgError('Variable %s in file %s already defined' % (var, varInfoFilename))

            self.info[var] = VarInfo(varName=var, desc=desc, obsolete=obsolete)





# This class represents the input file as a dictionary ('inputDict'). It reads the input file literally (removing comments and continuation symbols '/')
# The actual parsing is done in method 'getVal()'.

class Input:
    def __init__(self, fileName, varInfoList):

        self.inputDict = {}
        lineNum = 0
        li = ''
        for line in open(fileName):


            lineNum += 1

            li += line.strip()

            # Removes comments.
            while '!' in li:
                li = li[:li.find('!')].strip()

            if li != '' and li[-1] != '\\':
                # If it is a non empty line that doesn't continue in the next, adds it to the dictionary.

                # We are considering that = and == are quivalent
                i = li.find('==')
                i2 = i + 2
                if i < 0:
                    i = li.find('=')
                    i2 = i + 1
                if i < 0:
                    msgError('No equal sign found reading input file "%s", line %i: %s', (fileName, lineNum, li))

                self.inputDict[li[:i].strip()] = li[i2:].strip()
                li = ''
            else:
                li = li[:-1]   # Removes the continuation symbol '/'



        for var in varInfoList.obsolete:
            if var in self.inputDict.keys():
                msgInfo('Variable %s is now obsolete and its value is not taken into account. It can be safely removed from the input file' % var)

        self.varInfoList = varInfoList


    def printReport(self, filename = None):

        print('Input file contents report')
        print('--------------------------')
        print()

        for var in self.inputDict:

            val = self.inputDict[var]

            if (var not in self.varInfoList.info.keys()):
                msgInfo('Variable %s is not found in "input_params_info", where all variables must be described.' % var)
                print('%s = %s*' % (var, val))
            else:
                info = self.varInfoList.info[var]
                print('%s = %s [%s]: %s' % (var, val, info.additionalInfo, info.desc))



            if (var == 'VARNAME') and not exists(val):
                msgError('Netcdf meta information file %s not found' % val, 4)


    def getVal(self, varName, minVal = None, maxVal = None, dtype = None, count = None):
        # This function parses the value

        val = self.inputDict[varName]

        if dtype is None:
            if minVal is not None or maxVal is not None:
                msgError('minVal and maxVal are only valid when the variable has defined dtype. See input_params_info.dat (var = %s)' % varName)

            res = val


        else:
            # For typed variables.

            if count is None:
                count = 1
                values = [val]
            else:
                values = val.split()

            res = []

            for i in range(count):
                val = values[i]

                if (dtype == float) and ('d' in val.lower()):
                    # Fortran uses "d" like in 1.57d-2 for double floats.
                    # Substitutes the 'd' by 'e'
                    idx = val.lower().find('d')
                    val = val[:idx] + 'e' + val[idx+1:]

                if (dtype == bool):
                    if   val == 'T':
                        val = '1'       # This is converted to True in dtype(val)
                    elif val == 'F':
                        val = ''        # This is converted to False in dtype(val)

                try:
                    val = dtype(val)
                except:
                    msgError('Parsing the text "%s" into a variable of type %s. See input_params_info.dat (var = %s)' % (val, repr(dtype), varName))

                if minVal is not None and minVal >= val:
                    msgError('Value %s is smaller or equal than the minimum. See input_params_info.dat (var = %s)' % (val, varName))

                if maxVal is not None and maxVal <= val:
                    msgError('Value %s is larger or equal than the maximum. See input_params_info.dat (var = %s)' % (val, varName))

                res += [val]

            if count == 1:
                res = res[0]

        return res

