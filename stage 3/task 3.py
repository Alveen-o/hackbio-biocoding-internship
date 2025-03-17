
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score
import seaborn as sns

# Load the data
data = pd.read_csv('https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/drug_class_struct.txt', delim_whitespace=True)

# Separate features and docking scores
features = data.drop(columns=['score'])
docking_scores = data['score']

# Filter out non-numeric columns
features = features.select_dtypes(include=[np.number])

# Standardize the data
scaler = StandardScaler()
features_scaled = scaler.fit_transform(features)

# Perform PCA
pca = PCA(n_components=2)
principal_components = pca.fit_transform(features_scaled)

# Perform K-means clustering
kmeans = KMeans(n_clusters=5, random_state=42, n_init=10)
clusters = kmeans.fit_predict(principal_components)

# Create a DataFrame for visualization
pca_df = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'])
pca_df['Cluster'] = clusters
pca_df['Docking Score'] = docking_scores
print(pca_df)
# Plot chemical space with clusters
plt.figure(figsize=(10, 6))
sns.scatterplot(x='PC1', y='PC2', hue='Cluster', size='Docking Score', sizes=(20, 200), data=pca_df, palette='viridis')
plt.title('Chemical Space Clustering by Docking Score')
plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.legend(title='Cluster')
plt.show()

# Identify sub-cluster with the lowest docking scores
cluster_docking_scores = pca_df.groupby('Cluster')['Docking Score'].mean()
lowest_score_cluster = cluster_docking_scores.idxmin()
print(f'The cluster with the lowest average docking score is Cluster {lowest_score_cluster}')

# Regression model to predict docking score
regressor = LinearRegression()
regressor.fit(features_scaled, docking_scores)

# Predicting the docking score
predicted_scores = regressor.predict(features_scaled)

# Calculate performance metrics
mse = mean_squared_error(docking_scores, predicted_scores)
r2 = r2_score(docking_scores, predicted_scores)
print(f'Mean Squared Error: {mse}')
print(f'R-squared: {r2}')
