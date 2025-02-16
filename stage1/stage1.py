import numpy as np
import pandas as pd

# Genetic Code Dictionary
genetic_code = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
    'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W'
}

# DNA to Protein Function
def dna_to_protein(dna_sequence):
    """
    Translates a given DNA sequence into a protein sequence.
    Stops translation at a stop codon (*).

    Parameters:
        dna_sequence (str): DNA sequence in 5' to 3' direction.

    Returns:
        str: Translated protein sequence.
    """
    dna_sequence = dna_sequence.upper()
    protein_sequence = []

    for i in range(0, len(dna_sequence) - 2, 3):
        codon = dna_sequence[i:i + 3]
        amino_acid = genetic_code.get(codon, 'X')  # 'X' for unknown codons
        if amino_acid == '*':
            break
        protein_sequence.append(amino_acid)

    return ''.join(protein_sequence)
# Test cases
test_sequences = 'GCCATGGAGGCCCATCAGTTTATTAAGGCTCCTGGCATCACTACTGCTATTGAGCAGGCTGCTCTAGCAGCGGCCAACTCTGCCCTTGCGAATGCTGTGGTAGTTAGGCCTTTTCTCTCTCACCAGCAGATTGAGATCCTTATTAACCTAATGCAACCTCGCCAGCTTGTTTTCCGCCCCGAGGTTTTCTGGAACCATCCCATCCAGCGTGTTATCCATAATGAGCTGGAGCTTTACTGTCGCGCCCGTTCCGGCCGCTGCCTTGAAATTGGTGCCCACCCCCGCTCAATAAATGATAACCCTAATGTGGTCCACCGCTGCTTCCTCCGCCCTGCCGGGCGTGATGTTCAGCGTTGGTATACTGCCCCTACCCGCGGGCCGGCTGCTAATTGCCGGCGTTCCGCACTGCGCGGGCTCCCCGCTGCTGACCGCACTTACTGCTTCGACGGGTTTTCTGGCTGTAACTTTCCCGCCGAGACTGGCGTCGCCCTCTATTCTCTCCATGATATGTCACCATCTGATGTCGCCGAGGCTATGTTCCGCCATGGTATGACGCGGCTTTACGCTGCCCTCCACCTCCCGCCTGAGGTCCTGTTGCCCCCTGGCACATACCGCACCGCGTCGTACTTGCTGATCCATGACGGCAGGCGCGTTGTGGTGACGTATGAGGGTGACACTAGTGCTGGTTATAACCACGATGTTTCCAACCTGCGCTCCTGGATTAGAACCACTAAGGTTACCGGAGACCATCCTCTCGTCATTGAGCGGGTTAGGGCCATTGGCTGCCACTTTGTCCTCTTACTCACGGCAGCCCCGGAGCCATCACCTACGCCCTATGTTCCTTACCCCCGGTCTACCGAGGTCTATGTCCGATCGATCTTCGGCCCGGGTGGTACCCCCTCCCTATTTTCCAACCTCATGCTCCACTAAGTCGACCTTCCATGCTGTCCCTGCCCATATCTGGGACCGTCTCATGTTGTTCGGGGCCACCCTAGATGACCAAGCCTTTTGCTGCTCCCGCCTAATGACTTACCTCCGTGGCATTAGCTACAAGGTTACTGTGGGCACCCTTGTTGCCAATGAAGGCTGGAACGCCTCTGAGGTCGCTCTTACAGCTGTCATCACTGCCGCCTACCTTACCATCTGCCACCAGCGGTACCTCCGCACTCAGGCTATATCTAAGGGGATGCGCCGTCTGGAGCGGGAGCATGCTCAGAAGTTTATAACACGCCTCTACAGTTGGCTCTTTGAGAAGTCCGGCCGTGATTATATCCCCGGCCGTCAGTTGGAGTTCTACGCTCAGTGTAGGCGCTGGCTCTCGGCCGGCTTTTCATCTTGACCCACGGGTGTTGGTTTTTGATGAGTCGGCCCCCTGCCACTGTAGGACTGCGATTCGTAAGGCGGTCTCAAAGTTTTGCTGTTTTATGAAGTGGCTGGGCCAGGAGTGCACCTGTTTCCTACAACCTGCAGAAGGCGCCGTCGGCGACCAGGGCCATGACAACGAGGCCTATGAGGGGTCTGATGTCGACCCCGCTGAATCCGCTATTAGTGACATATCTGGGTCCTACGTCGTCCCTGGCACTGCCCTCCAACCGCTTTACCAAGCCCTTGACCTCCCCGCTGAGATTGTGGCTCGTGCAGGCCGGCTGACCGCCACAGTAAAGGTCTCCCAGGTCGACGGGCGGATCGATTGTGAGACCCTTCTCGGTAATAAAACCTTCCGCACGTCGTTTGTTGACGGGGCGGTTTTAGAGACTAATGGCCCAGAGCGCCACAATCTCTCTTTTGATGCCAGTCAGAGCACTATGGCCGCCGGCCCTTTCAGTCTCACCTATGCCGCCTCTGCTGCTGGGCTGGAGGTGCGCTATGTCGCTGCCGGGCTTGACCACCGGGCGGTTTTTGCCCCCGGCGTTTCACCCCGGTCAGCCCCTGGCGAGGTCACCGCCTTTTGTTCTGCCCTATACAGGTTTAATCGCGAGGCCCAGCGCCTTTCGCTCACCGGTAATTTTTGGTTCCATCCTGAGGGGCTCCTTGGCCCCTTTGCCCCGTTTTCCCCCGGGCATGTTTGGGAGTCGGCTAATCCATTCTGTGGAGAGAGCACACTTTACACCCGCACTTGGTCGGAGGTTGATGCTGTTTCTAGTCCAGCCCAGCCCGACTTAGGTTTTATATCTGAGCCTTCTATACCTAGTAGGGCCGCCACACTTACCCCGGCGGCCCCTCTACCCCCCCCTGCACGGAATCCTTCCCCTACTCCCTCTGCTCCGGCGCGTGGTGAGCCGGCTCCTGGCGCTACCGCCCGGGCCCCGGCCATAACCCACCAGGCGGCCCGGCATCGCCGCCTGCTCTTTACCTACCCGGATGGCTCTAAGGTATTCGCCGGCTCGCTGTTTGAGTCGACATGTACCTGGCTCGTTAACGCGTCTAATGTTG'
print(dna_to_protein(test_sequences))
# Logistic Growth Curve Function
def logistic_growth_curve(time_points, carrying_capacity, max_growth_rate, lag_duration, exp_duration):
    """
    Simulates a logistic population growth curve with randomized lag and exponential phases.

    Parameters:
        time_points (array): Array of time points.
        carrying_capacity (int): Maximum population size.
        max_growth_rate (float): Maximum growth rate.
        lag_duration (int): Duration of the lag phase.
        exp_duration (int): Duration of the exponential phase.

    Returns:
        array: Population values at each time point.
    """
    population_values = np.zeros_like(time_points, dtype=float)

    # Lag Phase
    lag_end = lag_duration
    exp_end = lag_end + exp_duration

    for i in range(len(time_points)):
        time = time_points[i]

        if time < lag_end:
            population_values[i] = 1  # Initial cell count in lag phase

        elif lag_end <= time < exp_end:
            population_values[i] = population_values[i - 1] * np.exp(max_growth_rate * (time - lag_end))

        else:  # Logistic Growth Phase
            N = population_values[i - 1]
            dNdt = max_growth_rate * ((carrying_capacity - N) / carrying_capacity) * N
            population_values[i] = N + dNdt

            # Ensure realistic population values
            population_values[i] = max(0, min(population_values[i], carrying_capacity))

    return population_values

# Generate 100 Growth Curves
num_curves = 100
time_points = np.linspace(0, 50, 100)  # Time points from 0 to 50
carrying_capacity = 1000  # Carrying Capacity
max_growth_rate = 0.2  # Maximum Growth Rate

# Create DataFrame
growth_df = pd.DataFrame({'Time': time_points})

for i in range(num_curves):
    lag_phase_duration = np.random.randint(5, 15)
    exp_phase_duration = np.random.randint(10, 20)
    population_curve = logistic_growth_curve(time_points, carrying_capacity, max_growth_rate, lag_phase_duration, exp_phase_duration)
    growth_df[f'Curve_{i+1}'] = population_curve

# Save to CSV
growth_df.to_csv("growth_curves.csv", index=False)
print(growth_df.head())

# Function to Find Time to Reach 80% of Carrying Capacity
def time_to_80_percent_capacity(time_points, population_values):
    """
    Determines the time at which population reaches 80% of carrying capacity.

    Parameters:
        time_points (array): Array of time values.
        population_values (array): Population values.

    Returns:
        float or None: Time when 80% of carrying capacity is reached, or None if not reached.
    """
    max_population = np.max(population_values)
    threshold_population = 0.8 * max_population

    indices = np.where(population_values >= threshold_population)[0]

    if indices.size > 0:
        return time_points[indices[0]]
    else:
        return None  # If 80% isn't reached

# Test 80% Capacity Function
time_test = np.linspace(0, 50, 100)  # Example time points
population_test = np.array([0, 10, 50, 100, 200, 300, 400, 500, 600, 700, 750, 780,
                            795, 800, 805, 810, 800, 790, 750, 700])  # Example population data

time_80 = time_to_80_percent_capacity(time_test, population_test)

if time_80 is not None:
    print(f"Time to reach 80% of carrying capacity: {time_80:.2f}")
else:
    print("Population did not reach 80% of carrying capacity within the given timeframe.")

# Hamming Distance Function
def hamming_distance(string1, string2):
    """
    Computes the Hamming distance between two strings.

    Parameters:
        string1 (str): First string.
        string2 (str): Second string.

    Returns:
        int: Number of differing characters between the two strings.
    """
    max_length = max(len(string1), len(string2))
    string1, string2 = string1.ljust(max_length), string2.ljust(max_length)  # Pad shorter string
    return sum(char1 != char2 for char1, char2 in zip(string1, string2))

# Test Hamming Distance
slack_username = "alvinnnnn"
x_handle = "alvinoooo"
print(f"Hamming Distance: {hamming_distance(slack_username, x_handle)}")
