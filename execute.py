import os
import json

numrows = 50000

# First we update the number of rows
# Step 1: Read the existing JSON file
with open('parameters.json', 'r') as file:
    data = json.load(file)

# Step 2: Modify the data structure
data['numrow'] = numrows

# Step 3: Write the updated JSON file
with open('parameters.json', 'w') as file:
    json.dump(data, file, indent=4)

numprocessers = 5

os.system('python proteins-generator.py ' + str(numrows))

os.system('python serial-proteins.py')

# os.system('mpiexec -n ' + str(numprocessers) + ' python mpi-proteins.py')