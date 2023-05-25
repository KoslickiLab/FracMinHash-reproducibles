import seaborn as sns
import matplotlib.pyplot as plt

# Create the figure and subplots
fig, axes = plt.subplots(1, 3, figsize=(12, 4))

# Plot on the first subplot using Seaborn
sns.lineplot([1, 2, 3], [4, 5, 6], ax=axes[0])
axes[0].set_title('Plot 1')

# Plot on the second subplot using Seaborn
sns.scatterplot([1, 2, 3], [4, 5, 6], ax=axes[1])
axes[1].set_title('Plot 2')

# Plot on the third subplot using Seaborn
sns.barplot([1, 2, 3], [4, 5, 6], ax=axes[2])
axes[2].set_title('Plot 3')

# Add a common title for the figure
fig.suptitle('Three Subplots in One Row')

# Adjust the spacing between subplots
plt.tight_layout()

# Display the figure
plt.show()
