import sys
import os


# exePath = r'D:\projects\src\oceangpu'
# exePath = r'/home/jo.gonzalez/src/Sailfish'

# compilationOptions = ('-default-device', '--restrict', '--std=c++17', )
compilationOptions = ('-default-device', '--restrict', '--std=c++17', )
GPUMUL = 1

LwrtInfo = True

filePath = os.path.dirname(os.path.abspath(__file__))
exePath, _ = os.path.split(filePath)

blockSize = 1024

# Execution termination errors.
exit_flag = {0: 'No error', 1: 'Blows up', 2: 'Input error', 3: 'Output error', 4: 'IO error',
             5: 'Configuration error', 6: 'Partition error', 7: 'Illegal input parameter',
             8: 'Fatal algorith result', 9: 'Coupling error'}



def exitProgram(str = None):
    if str is not None:
        print(str)

    print('PROGRAM ENDED SUCCESFULY')
    sys.exit(0)


def msgError(str, error = 1000):
    print('ERROR: ', str)
    sys.exit(error)


def msgInfo(str, level = 0):
    print('INFO: ', str )