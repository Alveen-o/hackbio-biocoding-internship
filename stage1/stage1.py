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


def dna_to_protein(seq):
    """
    Translates a given DNA sequence into a protein sequence.
    Stops translation at a stop codon (*).

    Parameters:
        seq (str): DNA sequence (assumed to be in 5' to 3' direction)

    Returns:
        str: Translated protein sequence
    """
    seq = seq.upper()  # Convert sequence to uppercase
    protein_seq = []  # Use list for efficiency

    for i in range(0, len(seq) - 2, 3):  # Ensure valid triplets
        codon = seq[i:i + 3]
        amino_acid = genetic_code.get(codon, 'X')  # 'X' for unknown codons
        if amino_acid == '*':  # Stop translation at stop codon
            break
        protein_seq.append(amino_acid)

    return ''.join(protein_seq)  # Join list into string


# Test cases
test_sequences = [
    'TCGCGCACGCTGATCGTGGGGTGA',
    'agtaaaactttaattgttggttaa',
    'GCCATGGAGGCCCATCAGTTTATTAAGGCTCCTGGCATCACTACTGCTATTGAGCAGGCTGCTCTAGCAGCGGCCAACTCTGCCCTTGCGAATGCTGTGGTAGTTAGGCCTTTTCTCTCTCACCAGCAGATTGAGATCCTTATTAACCTAATGCAACCTCGCCAGCTTGTTTTCCGCCCCGAGGTTTTCTGGAACCATCCCATCCAGCGTGTTATCCATAATGAGCTGGAGCTTTACTGTCGCGCCCGTTCCGGCCGCTGCCTTGAAATTGGTGCCCACCCCCGCTCAATAAATGATAACCCTAATGTGGTCCACCGCTGCTTCCTCCGCCCTGCCGGGCGTGATGTTCAGCGTTGGTATACTGCCCCTACCCGCGGGCCGGCTGCTAATTGCCGGCGTTCCGCACTGCGCGGGCTCCCCGCTGCTGACCGCACTTACTGCTTCGACGGGTTTTCTGGCTGTAACTTTCCCGCCGAGACTGGCGTCGCCCTCTATTCTCTCCATGATATGTCACCATCTGATGTCGCCGAGGCTATGTTCCGCCATGGTATGACGCGGCTTTACGCTGCCCTCCACCTCCCGCCTGAGGTCCTGTTGCCCCCTGGCACATACCGCACCGCGTCGTACTTGCTGATCCATGACGGCAGGCGCGTTGTGGTGACGTATGAGGGTGACACTAGTGCTGGTTATAACCACGATGTTTCCAACCTGCGCTCCTGGATTAGAACCACTAAGGTTACCGGAGACCATCCTCTCGTCATTGAGCGGGTTAGGGCCATTGGCTGCCACTTTGTCCTCTTACTCACGGCAGCCCCGGAGCCATCACCTACGCCCTATGTTCCTTACCCCCGGTCTACCGAGGTCTATGTCCGATCGATCTTCGGCCCGGGTGGTACCCCCTCCCTATTTTCCAACCTCATGCTCCACTAAGTCGACCTTCCATGCTGTCCCTGCCCATATCTGGGACCGTCTCATGTTGTTCGGGGCCACCCTAGATGACCAAGCCTTTTGCTGCTCCCGCCTAATGACTTACCTCCGTGGCATTAGCTACAAGGTTACTGTGGGCACCCTTGTTGCCAATGAAGGCTGGAACGCCTCTGAGGTCGCTCTTACAGCTGTCATCACTGCCGCCTACCTTACCATCTGCCACCAGCGGTACCTCCGCACTCAGGCTATATCTAAGGGGATGCGCCGTCTGGAGCGGGAGCATGCTCAGAAGTTTATAACACGCCTCTACAGTTGGCTCTTTGAGAAGTCCGGCCGTGATTATATCCCCGGCCGTCAGTTGGAGTTCTACGCTCAGTGTAGGCGCTGGCTCTCGGCCGGCTTTTCATCTTGACCCACGGGTGTTGGTTTTTGATGAGTCGGCCCCCTGCCACTGTAGGACTGCGATTCGTAAGGCGGTCTCAAAGTTTTGCTGTTTTATGAAGTGGCTGGGCCAGGAGTGCACCTGTTTCCTACAACCTGCAGAAGGCGCCGTCGGCGACCAGGGCCATGACAACGAGGCCTATGAGGGGTCTGATGTCGACCCCGCTGAATCCGCTATTAGTGACATATCTGGGTCCTACGTCGTCCCTGGCACTGCCCTCCAACCGCTTTACCAAGCCCTTGACCTCCCCGCTGAGATTGTGGCTCGTGCAGGCCGGCTGACCGCCACAGTAAAGGTCTCCCAGGTCGACGGGCGGATCGATTGTGAGACCCTTCTCGGTAATAAAACCTTCCGCACGTCGTTTGTTGACGGGGCGGTTTTAGAGACTAATGGCCCAGAGCGCCACAATCTCTCTTTTGATGCCAGTCAGAGCACTATGGCCGCCGGCCCTTTCAGTCTCACCTATGCCGCCTCTGCTGCTGGGCTGGAGGTGCGCTATGTCGCTGCCGGGCTTGACCACCGGGCGGTTTTTGCCCCCGGCGTTTCACCCCGGTCAGCCCCTGGCGAGGTCACCGCCTTTTGTTCTGCCCTATACAGGTTTAATCGCGAGGCCCAGCGCCTTTCGCTCACCGGTAATTTTTGGTTCCATCCTGAGGGGCTCCTTGGCCCCTTTGCCCCGTTTTCCCCCGGGCATGTTTGGGAGTCGGCTAATCCATTCTGTGGAGAGAGCACACTTTACACCCGCACTTGGTCGGAGGTTGATGCTGTTTCTAGTCCAGCCCAGCCCGACTTAGGTTTTATATCTGAGCCTTCTATACCTAGTAGGGCCGCCACACTTACCCCGGCGGCCCCTCTACCCCCCCCTGCACGGAATCCTTCCCCTACTCCCTCTGCTCCGGCGCGTGGTGAGCCGGCTCCTGGCGCTACCGCCCGGGCCCCGGCCATAACCCACCAGGCGGCCCGGCATCGCCGCCTGCTCTTTACCTACCCGGATGGCTCTAAGGTATTCGCCGGCTCGCTGTTTGAGTCGACATGTACCTGGCTCGTTAACGCGTCTAATGTTG'
]
import numpy as np
import pandas as pd

def growth_curve(time, max_pop_size, r_max, lag_dur, exp_dur):
    population = np.zeros_like(time, dtype=float)

    # Lag Phase
    lag_end = lag_dur
    exp_end = lag_end + exp_dur

    for i in range(len(time)):
        t = time[i]

        if t < lag_end:
            population[i] = 1  # Initial cell count in lag phase

        elif lag_end <= t < exp_end:
            population[i] = population[i - 1] * np.exp(r_max * (t - lag_end))

        else:  # Logistic Growth Phase
            N = population[i - 1]
            dNdt = r_max * ((max_pop_size - N) / max_pop_size) * N
            population[i] = N + dNdt

            # Ensure realistic population values
            if population[i] < 0:
                population[i] = 0
            elif population[i] > max_pop_size:
                population[i] = max_pop_size

    return population

# Parameters
num_curves = 100
time = np.linspace(0, 50, 100)  # Time points from 0 to 50
K = 1000  # Carrying capacity
r_max = 0.2  # Growth rate

# Create DataFrame
df = pd.DataFrame({'Time': time})

# Generate Growth Curves
for i in range(num_curves):
    lag_phase_add = np.random.randint(5, 15)
    exp_phase_add = np.random.randint(10, 20)
    population_curve = growth_curve(time, K, r_max, lag_phase_add, exp_phase_add)
    df[f'Curve_{i+1}'] = population_curve

# Save CSV
df.to_csv("growth_curves.csv", index=False)
print(df.head())

def growth_curve(time, max_pop_size, r_max, lag_dur, exp_dur):
    population = np.zeros_like(time, dtype=float)

    # Lag phase
    lag_end = np.random.randint(5, 15) + lag_dur
    for i in range(len(time)):
        if time[i] < lag_end:
            population[i] = 1
        else:
            break

    # Log phase
    exp_end = np.random.randint(lag_end + 10, lag_end + 20) + exp_dur
    for i in range(len(time)):
        if time[i] >= lag_end and time[i] < exp_end:
            population[i] = 1 * np.exp(r_max * (time[i] - lag_end))  # start log growth from 1 cell
        elif time[i] >= exp_end:
            break

    # Logistic growth
    for i in range(len(time)):
        if time[i] >= exp_end:
            N = population[i - 1]
            dNdt = r_max * ((max_pop_size - N) / max_pop_size) * N
            population[i] = N + dNdt

            if population[i] < 0:
                population[i] = 0  # prevent negative cell numbers
            elif population[i] > max_pop_size:
                population[i] = max_pop_size  # prevent overgrowth value

    return population

# Generate 100 growth curves
num_curves = 100
time = np.linspace(0, 50, 100)  # Time points from 0 to 50
K = 1000  # Example maximum population (K)
r_max = 0.2  # Example maximum growth rate

df = pd.DataFrame({'Time': time})

for i in range(num_curves):
    lag_phase_add = np.random.randint(0, 5)
    exp_phase_add = np.random.randint(0, 5)
    population_curve = growth_curve(time, K, r_max, lag_phase_add, exp_phase_add)
    df[f'Curve_{i+1}'] = population_curve

print(df.head())
df.to_csv("growth_curves.csv", index=False)  # Save to CSV


# Function to find the time at which population reaches 80% of max
def time_80percent(time, population):
    max_pop = np.max(population)
    target_pop = 0.8 * max_pop

    indices = np.where(population >= target_pop)[0]

    if indices.size > 0:
        return time[indices[0]]
    else:
        return None  # If 80% isn't reached

# Test for function 3:
time = np.linspace(0, 50, 100)  # Random time points
population = np.array([0, 10, 50, 100, 200, 300, 400, 500, 600, 700, 750, 780,
                       795, 800, 805, 810, 800, 790, 750, 700])  # Random population data

time_80 = time_80percent(time, population)

if time_80 is not None:
    print(f"Time to reach 80% of carrying capacity: {time_80:.2f}")
else:
    print("Population did not reach 80% of carrying capacity within the given timeframe.")


# Function for calculating the Hamming distance
def hamm_dist(s1, s2):
    max_len = max(len(s1), len(s2))
    s1, s2 = s1.ljust(max_len), s2.ljust(max_len)  # Pad shorter string with spaces
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

# Test Hamming Distance
slackname = "alvinnnnn"
xname = "alvinoooo"
print(f"Hamming Distance: {hamm_dist(slackname, xname)}")