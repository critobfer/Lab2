import time
import csv
import matplotlib.pyplot as plt
import numpy as np
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

def read_data():
    with open('proteins.csv', newline='') as csvfile:
        csvreader = csv.DictReader(csvfile)
        data = []
        for row in csvreader:
            id = int(row['structureId'])
            sequence = row['sequence']
            data.append((id, sequence))
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
    print('The protein with max ocurrences has de id', id, 'and it has', n_occurrences, 'ocurrence of the pattern', patter_in_upper_case, '.')

def main():
    pattern = input("Please insert the pattern you want to search for: ")
    patter_in_upper_case = pattern.upper()

    t0 = time.time()
    proteins = read_data()

    # Split the data among processes
    chunk_size = len(proteins) // size
    start = rank * chunk_size
    end = (rank + 1) * chunk_size if rank != size - 1 else len(proteins)

    # Each process works on its portion of the data
    local_proteins = proteins[start:end]
    local_ocurrences_dict = search_pattern(local_proteins, patter_in_upper_case)

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
        max_occurrences_protein(ocurrences_dict_sorted, patter_in_upper_case)

if __name__ == "__main__":
    main()

