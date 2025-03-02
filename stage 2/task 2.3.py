import pandas as pd 
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# Load the pesticide treatment data
data_url = 'https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/Pesticide_treatment_data.txt'
treatment_df = pd.read_csv(data_url, sep='\t')

# Select only numeric metabolite columns for subtraction
numeric_metabolites = treatment_df.select_dtypes(include=['int64', 'float64']).columns

# Initialize an empty DataFrame for storing delta responses
delta_response_df = pd.DataFrame()

# Calculate Wild-Type metabolic response: subtract row 3 from row 0 for all numeric metabolites
delta_response_df["Delta_WT"] = treatment_df.loc[0, numeric_metabolites] - treatment_df.loc[3, numeric_metabolites] 

# Calculate Mutant metabolic response: subtract row 6 from row 4 for all numeric metabolites
delta_response_df["Delta_MT"] = treatment_df.loc[4, numeric_metabolites] - treatment_df.loc[6, numeric_metabolites]

# Calculate the distance of each metabolite from the y=x line (diagonal distance)
delta_response_df["Diagonal_Distance"] = np.abs(delta_response_df["Delta_MT"] - delta_response_df["Delta_WT"]) / np.sqrt(2)
print(delta_response_df["Diagonal_Distance"])

# Classify metabolites based on a defined residual cutoff:
cutoff_classification = []
for distance in delta_response_df["Diagonal_Distance"]:
    if -0.3 <= distance <= 0.3:
        classification = 'inside'
    else:
        classification = 'outside'
    cutoff_classification.append(classification)

delta_response_df["Cutoff_Status"] = cutoff_classification
print(f"These are metabolites outside the residual cut-off:\n{delta_response_df[delta_response_df['Cutoff_Status'] == 'outside'].index}")
print(f"These are metabolites inside the residual cut-off:\n{delta_response_df[delta_response_df['Cutoff_Status'] == 'inside'].index}")

# Create a scatter plot to visualize the metabolic responses
fig, scatter_ax = plt.subplots()

sns.scatterplot(
    data=delta_response_df, 
    x="Delta_MT", 
    y="Delta_WT", 
    hue="Cutoff_Status", 
    palette={'inside': 'grey', 'outside': 'salmon'}, 
    ax=scatter_ax
)

# Add a reference line (y = x)
ref_x = []
ref_y = []
def generate_reference_line():
    for i in range(-3, 3):
        ref_x.append(i)
        ref_y.append(i)
generate_reference_line()
scatter_ax.plot(ref_x, ref_y)

scatter_ax.set_title('Wild-Type and Mutant Metabolic Responses')
scatter_ax.set_xlabel('Delta MT')
scatter_ax.set_ylabel('Delta WT')
plt.show()

# Comment: Grey metabolites (inside) with low residual show similar response between WT and MT plants,
# while Salmon metabolites (outside) with high residual show different responses.
# There is a significant difference in metabolic response between wild-type (WT) and mutant (MT) plants.

# Select metabolites outside the residual cutoff (first 6 for visualization)
metabolites_outside_cutoff = delta_response_df[delta_response_df["Cutoff_Status"] == 'outside'].head(6)

# Define time points for the x-axis (in hours)
time_points = [0, 8, 24]

# Create a new figure for the time-course plot
fig, time_ax = plt.subplots()

# Loop over each selected metabolite and plot its time-course metabolic response (rows 1 to 3 correspond to WT)
for metabolite in metabolites_outside_cutoff.index:
    wt_response = treatment_df.loc[1:3, metabolite].values
    time_ax.plot(time_points, wt_response, marker='o', linestyle='-', label=f"WT {metabolite}")

time_ax.set_title('Metabolic Response Over Time for Selected Metabolites')
time_ax.set_xlabel('Time (hours)')
time_ax.set_ylabel('Metabolic Response')
time_ax.legend()

plt.show()
