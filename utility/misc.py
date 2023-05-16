import sys



# Execution termination errors.
exit_flag = {0: 'No error', 1: 'Blows up', 2: 'Input error', 3: 'Output error', 4: 'IO error',
             5: 'Configuration error', 6: 'Partition error', 7: 'Illegal input parameter',
             8: 'Fatal algorith result', 9: 'Coupling error'}


def msgError(str, error = 1000):
    print('ERROR: ', str)
    sys.exit(error)


def msgInfo(str, level = 0):
    print('INFO: ', str )