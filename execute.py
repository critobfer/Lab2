import os

numrows = 50000

numprocessers = 5

os.system('python proteins-generator.py ' + str(numrows))

os.system('python serial-proteins.py')

# os.system('mpiexec -n ' + str(numprocessers) + ' python mpi-proteins.py')