import os

numrows = 50000
os.system('python proteins-generator.py ' + str(numrows))

os.system('python serial-proteins.py')

os.system('python mpi-proteins.py')