import time
import csv
import matplotlib.pyplot as plt
import numpy as np
from mpi4py import MPI
from itertools import islice


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

def read_data(numrows):
    subset_size = numrows // size

    with open('proteins.csv', newline='') as csvfile:
        csvreader = csv.DictReader(csvfile)

        for _ in range(rank * subset_size):
            next(csvreader)
        
        data_subset = []
        for _ in range(subset_size):
            try:
                row = next(csvreader)
                data_subset.append((int(row['structureId']), row['sequence']))
            except StopIteration:
                break
    return data_subset
    
def search_pattern(data_chunk, pattern):
    ocurrences_dict = {}
    for (id, sequence) in data_chunk:
        ocurrences = sequence.count(pattern)
        if ocurrences != 0:
            ocurrences_dict[id] = ocurrences
    return ocurrences_dict

def plot_histogram(ocurrences_dict, n):
    ocurrences_dict_sorted = dict(sorted(ocurrences_dict.items(), key=lambda item: item[1], reverse=True))
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

def max_occurrences_protein(ocurrences_dict_sorted, pattern_in_upper_case):
    (id, n_occurrences) = next(iter(ocurrences_dict_sorted.items()))
    print(f'The protein with max occurrences has the id {id} and it has {n_occurrences} occurrences of the pattern {pattern_in_upper_case}.')

def main():
    if rank == 0:
        # Only the root process collects user input
        pattern = input("Please insert the pattern you want to search for: ")
    else:
        pattern = None

    # Broadcast the pattern to all processes
    pattern = comm.bcast(pattern, root=0)
    pattern_in_upper_case = pattern.upper()

    t0 = time.time()
    local_data = read_data(500000)

    # Search for the pattern in each process's data chunk
    local_ocurrences_dict = search_pattern(local_data, pattern_in_upper_case)

    # Gather results from all processes
    ocurrences_dicts = comm.gather(local_ocurrences_dict, root=0)

    t1 = time.time()

    if rank == 0:
        # Combine results from all processes
        ocurrences_dict = {}
        for d in ocurrences_dicts:
            ocurrences_dict.update(d)
        
        print(f"Time elapsed: {t1 - t0} seconds")
        ocurrences_dict_sorted = plot_histogram(ocurrences_dict, 10)
        max_occurrences_protein(ocurrences_dict_sorted, pattern_in_upper_case)

if __name__ == "__main__":
    main()




