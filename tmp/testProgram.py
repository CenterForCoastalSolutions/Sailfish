import sys
sys.path.insert(0, '../modules')
sys.path.insert(0, '../utility')


import mod_io
import io_utils



# Reads the file with meta information about the variables.
varInfoList = io_utils.VarInfoList('../modules')

# Reads the input file.
input = io_utils.Input('../utility/test.in')

pass