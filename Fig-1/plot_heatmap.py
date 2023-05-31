import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl


mpl.rcParams['font.family'] = 'Arial'

def generate_containment(p: float, k: int) -> float:
    return (1 - p) ** k

# Define the range of p and k values
#p_values = np.array([0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10])
p_values = np.array([0.001, 0.01, 0.02, 0.04, 0.06, 0.08, 0.10])
#p_values = np.array([0.001]) + np.linspace(0.01, 0.10, num=9)
#k_values = np.array([21, 31, 41, 51, 61, 71, 81, 91, 100])
k_values = np.array([20, 40, 60, 80, 100])

# Generate a 2D array of containment probabilities for all combinations of p and k
probabilities = np.zeros((len(k_values), len(p_values)))
for i, k in enumerate(k_values):
    for j, p in enumerate(p_values):
        probabilities[i, j] = generate_containment(p, k)

# Create a heatmap plot of the containment probabilities
fig, ax = plt.subplots(figsize=(6.5, 4))

im = ax.imshow(probabilities, cmap='GnBu_r')

# Add x and y axis labels and a colorbar legend
ax.set_xticks(np.arange(len(p_values)))
ax.set_yticks(np.arange(len(k_values)))
ax.set_xticklabels([f'{p_values[0]:.3f}'] + [f'{p:.2f}' for p in p_values[1:]])
ax.set_yticklabels(k_values)
ax.set_xlabel('Probability of point mutation ($p$)')
ax.set_ylabel('$k$-mer size ($k$)')

plt.title('Expected containment\nacross $k$-mer sizes and mutation rates')

plt.subplots_adjust(left=0)


# Add a colorbar and specify its position using a separate axes object
cax = fig.add_axes([0.82, 0.11, 0.03, 0.77])  # [left, bottom, width, height]

cbar = ax.figure.colorbar(im, ax=ax, cax=cax)
cbar.ax.set_ylabel('Containment of $k$-mers of \nmutated sequence in the original', rotation=90)

# Add the probability values to the heatmap
for i in range(len(k_values)):
    for j in range(len(p_values)):
        ax.text(j, i, f'{probabilities[i, j]:.3f}', ha='center', va='center')

# Set the title of the plot
#

# Display the plot
plt.savefig('containment_vs_k_p_heatmap.pdf')
#plt.show()
