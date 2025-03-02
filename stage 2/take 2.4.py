import pandas as pd
import matplotlib.pyplot as plt

sift_data = 'https://raw.githubusercontent.com/HackBio-Internship/public_datasets/main/R/datasets/sift.tsv'
fold_data = 'https://raw.githubusercontent.com/HackBio-Internship/public_datasets/main/R/datasets/foldX.tsv'

# Read file with proper whitespace handling, when \t and ' ' was used, the dataframe was not properly print
sift_df = pd.read_csv(sift_data, sep=r'\s+')
fold_df = pd.read_csv(fold_data, sep=r'\s+')
# print(sift_df) # --> check the printing of dataframe

sift_df["specific_Protein_aa"] = sift_df["Protein"] + '_' + sift_df["Amino_Acid"]
fold_df["specific_Protein_aa"] = fold_df["Protein"] + '_' + fold_df["Amino_Acid"]
# print(sift_df, '\n', fold_df) # --> check if the columns were added correctly

# merge datasets with specific_Protein_aa
merged_df = pd.merge(fold_df, sift_df, on="specific_Protein_aa")

# Drop duplicate columns
merged_df = merged_df.drop(columns=["Protein_y", "Amino_Acid_y"])
# Rename remaining columns for clarity
merged_df = merged_df.rename(columns={"Protein_x": "Protein", "Amino_Acid_x": "Amino_Acid"})
# print(merged_df) # --> check

# A SIFT Score below 0.05 is deleterious
delet_sift = merged_df[merged_df["sift_Score"] < 0.05]
# A FoldX score greater than 2 kCal/mol is deleterious
delet_fold = merged_df[merged_df["foldX_Score"] > 2]
# print(delet_fold, delet_sift) --> check

mutants = []
for i in delet_sift["specific_Protein_aa"]:
    if i in delet_fold["specific_Protein_aa"].values:
        mutants.append(i)
print(mutants)
# print(len(mutants)) # --> check

""" Study the amino_acid_substitution_nomenclature
(https://atlasgeneticsoncology.org/teaching/30067/nomenclature-for-the-description-of-mutations-and-other-sequence-variations#section-3)"""

# Extract the original (wild-type) amino acid
merged_df["Wild_Type_AA"] = merged_df["Amino_Acid"].str[0]
# print(merged_df) # --> check

# Generate frequency table:
freq_aa = merged_df["Wild_Type_AA"].value_counts()
print(freq_aa)

# Generate a barplot:
plt.bar(x=freq_aa.index, height=freq_aa.values)
plt.xlabel("Wild Type AA")
plt.ylabel("Frequency")
plt.show()

# Generate a pie-chart:
plt.pie(freq_aa.values, labels=freq_aa.index, textprops={'fontsize': 6})
plt.show()

