from os.path import join
from misc import *
# This class contains general information about the input variables. This information is stored in the text dile input_vars.dat

class VarInfo:
    def __init__(self, varName, desc, obsolete):
        self.varName  = varName
        self.desc     = desc
        self.obsolete = obsolete


# This class contains information about the variables as given by the text file 'input_var.dat'. All variables used in the input file MUST be included
# in 'input_var.dat'.
# The intention here is to simplify the code and make it more flexible by removing this information from the code itself. Think for example in how
# ROMS code prints out the values read from the input file.

class VarInfoList:
    def __init__(self, path = '.'):

        self.obsolete = []
        varInfoFilename = join(path, 'input_vars.dat')
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
    def __init__(self, fileName):

        varInfoList = VarInfoList()

        self.inputDict = {}
        lineNum = 0
        for line in open(fileName):

            li  = ''
            lineNum += 1

            if len(li) > 0 and li[-1] == '\\':
                li += line.strip()
            else:
                li = line.strip()

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
                    print('No equal sign found reading input file "%s", line %i: %s', (filenameInput, lineNum, li))

                self.inputDict[li[:i].strip()] = li[i2:].strip()
            else:
                li = li[:-1]

            for var in varInfoList.obsolete:
                if var in self.inputDict.keys():
                    print('Variable %s is now obsolete and its value is not taken into account', var)





    def getVal(self):
        pass