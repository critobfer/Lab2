import time
import csv
import matplotlib.pyplot as plt
import numpy as np
from mpi4py import MPI

# Initialize MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

def read_data():
    if rank == 0:
        with open('proteins.csv', newline='') as csvfile:
            csvreader = csv.DictReader(csvfile)
            data = [row for row in csvreader]
    else:
        data = None

    # Broadcast the data to all processes
    data = comm.bcast(data, root=0)
    return data

def search_pattern(data, pattern):
    ocurrences_dict = {}
    for (id, sequence) in data:
        ocurrences = sequence.count(pattern)
        if ocurrences != 0:
            ocurrences_dict[id] = ocurrences
    return ocurrences_dict

def plot_histogram(ocurrences_dict, n):
    ocurrences_dict_sorted = dict(sorted(ocurrences_dict.items(), key=lambda item: item[1], reverse=True))

    # Extract x and y values into separate lists
    x_values = list(ocurrences_dict_sorted.keys())[:n]
    y_values = list(ocurrences_dict_sorted.values())[:n]

    x_positions = np.arange(len(x_values))
    plt.bar(x_positions, y_values, align='center')
    plt.xticks(x_positions, x_values, rotation=45)
    plt.xlabel('Protein id')
    plt.ylabel('Occurrences')
    plt.title(f"Histogram of the {n} proteins with more matches")
    plt.show()

    return ocurrences_dict_sorted

def max_occurrences_protein(ocurrences_dict_sorted, patter_in_upper_case):
    (id, n_occurrences) = next(iter(ocurrences_dict_sorted.items()))

    if rank == 0:
        print('The protein with max occurrences has de id', id, 'and it has', n_occurrences, 'occurrence of the pattern', patter_in_upper_case, '.')

def main():
    # Reading and changing to uppercase the pattern
    if rank == 0:
        pattern = input("Please insert the pattern you want to search for: ")
        patter_in_upper_case = pattern.upper()
    else:
        patter_in_upper_case = None

    # Broadcast the pattern to all processes
    patter_in_upper_case = comm.bcast(patter_in_upper_case, root=0)

    # Measure time
    t0 = time.time()
    proteins = read_data()
    local_ocurrences = search_pattern(proteins, patter_in_upper_case)

    # Gather local occurrences to the root process
    all_ocurrences = comm.gather(local_ocurrences, root=0)
    t1 = time.time()

    if rank == 0:
        # Merge the dictionaries from all processes
        ocurrences_dict = {}
        for local_dict in all_ocurrences:
            ocurrences_dict.update(local_dict)

        print(f"Time elapsed: {t1-t0} seconds")
        ocurrences_dict_sorted = plot_histogram(ocurrences_dict, 10)
        max_occurrences_protein(ocurrences_dict_sorted, patter_in_upper_case)

if __name__ == "__main__":
    main()
