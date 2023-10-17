import os
import json

with open('parameters.json', 'r') as archivo:
    param = json.load(archivo)

numrows = param['numrow']

numprocessers = param['numprocesses']

os.system('python proteins-generator.py ' + str(numrows))

os.system('python serial-proteins.py')

os.system('mpiexec -n ' + str(numprocessers) + ' python mpi-proteins.py')