import time 
import csv
import matplotlib.pyplot as plt
import numpy as np

def read_data():
    with open('proteins.csv', newline='') as csvfile:
        # Read the CSV file using DictReader
        csvreader = csv.DictReader(csvfile)
        
        # Initialize a list to store the tuples
        data = []
        
        # Iterate through the rows of the CSV file
        for row in csvreader:
            # Get the values of the 'structureId' and 'sequence' columns
            id = int(row['structureId'])
            sequence = row['sequence']
        
            # Create a tuple with the values and add it to the data list
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
    # Sort the dictionary to take only the 10 proteins with more matches, so we order by value
    ocurrences_dict_sorted =  dict(sorted(ocurrences_dict.items(), key=lambda item: item[1], reverse=True))
    # Extract x and y values into separate lists
    x_values = list(ocurrences_dict_sorted.keys())[:n]
    y_values = list(ocurrences_dict_sorted.values())[:n]

    # We set up the positions of the bar
    x_positions = np.arange(len(x_values))
    # Create a scatter plot using matplotlib.pyplot
    plt.bar(x_positions, y_values, align='center')
    #In the bar position we put the index
    plt.xticks(x_positions, x_values, rotation=45)
    # Set labels, title
    plt.xlabel('Protein id')
    plt.ylabel('Occurrences')
    plt.title(f"Histogram of the {n} proteins with more matches")
    plt.show()

    return ocurrences_dict_sorted

def max_occurrences_protein(ocurrences_dict_sorted, patter_in_upper_case):
    # We take the first element of the dictionary
    (id, n_occurrences) = next(iter(ocurrences_dict_sorted.items()))

    print('The protein with max ocurrences has de id', id, 'and it has', n_occurrences, 'ocurrence of the pattern', patter_in_upper_case, '.')

def main():
    # Reading and changing to uppercase the pattern 
    pattern = input("Please insert the pattern you want to search for: ")
    patter_in_upper_case = pattern.upper()

    # Measure time
    t0 = time.time()
    proteins = read_data()
    ocurrences_dict = search_pattern(proteins, patter_in_upper_case)
    t1 = time.time()
    # End measure time
    print(f"Time elapsed: {t1-t0} seconds")
    ocurrences_dict_sorted = plot_histogram(ocurrences_dict, 10)
    max_occurrences_protein(ocurrences_dict_sorted, patter_in_upper_case)


if __name__ == "__main__":
    main()