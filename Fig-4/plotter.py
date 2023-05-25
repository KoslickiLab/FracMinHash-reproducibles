import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import linregress
import matplotlib as mpl
from scipy.stats.stats import pearsonr
import numpy as np
from sklearn.metrics import mean_squared_error

mpl.rcParams['font.family'] = 'Arial'

fig, axes = plt.subplots(2, 2, figsize=(6.5, 5.5))

# read the CSV file
df = pd.read_csv("combinedANI.csv")

def rmse(predictions, targets):
    predictions = np.array(predictions)
    targets = np.array(targets)
    return np.sqrt(((predictions - targets) ** 2).mean())

colors = {
'pyani':'green',
'MashANI':'red',
'avg_cANI_k21_sc10':'blue',
'fastani':'orange'
}

axes[0][0].grid(alpha=0.2)
axes[0][1].grid(alpha=0.2)
axes[1][0].grid(alpha=0.2)
axes[1][1].grid(alpha=0.2)

# create the scatter plot on first subplot
sns.scatterplot(data=df, x="orthoANI", y="avg_cANI_k21_sc10", alpha=0.5, label="FracMinHash", color=colors['avg_cANI_k21_sc10'], ax=axes[0][0])
axes[0][0].set_title('FMH vs fastani')
axes[0][0].set_xlabel('OrthoANI')
axes[0][0].set_ylabel('FMH ANI')

sns.scatterplot(data=df, x="orthoANI", y="fastani", alpha=0.5, label="fastani", color=colors['fastani'], ax=axes[0][1])
#axes[0].set_title('FMH vs fastani')
axes[0][1].set_xlabel('OrthoANI')
axes[0][1].set_ylabel('fastani')

# Plot on the second subplot using Seaborn
sns.scatterplot(data=df, x="orthoANI", y="pyani", alpha=0.5, label="pyani", color=colors['pyani'], ax=axes[1][0])
#axes[1].set_title('FMH vs pyani')
axes[1][0].set_xlabel('OrthoANI')
axes[1][0].set_ylabel('pyani')

# Plot on the third subplot using Seaborn
# create the scatter plot on first subplot
sns.scatterplot(data=df, x="orthoANI", y="MashANI", alpha=0.5, label="Mash", color=colors['MashANI'], ax=axes[1][1])
#axes[2].set_title('FMH vs MashANI')
axes[1][1].set_xlabel('OrthoANI')
axes[1][1].set_ylabel('MashANI')

# Add a common title for the figure
fig.suptitle('OrthoANI compared to different ANI tools')

# TODO
# add title and labels
#plt.title("OrthoANI vs ANI")
#plt.xlabel("OrthoANI")
#plt.ylabel("ANI using different tools")

# add linear regression line
x_indices = {
    "pyani":1, "avg_cANI_k21_sc10":0
}
y_indices = {
    "pyani":0, "avg_cANI_k21_sc10":0
}

for col in ["pyani", "avg_cANI_k21_sc10"]:
    slope, intercept, r_value, p_value, std_err = linregress(df["orthoANI"], df[col])
    #if col == 'avg_cANI_k21_sc10':
        #slope -= 0.01
        #intercept += 0.01
    axes[x_indices[col]][y_indices[col]].plot(df["orthoANI"], intercept + slope * df["orthoANI"], color=colors[col], alpha=0.7)

min_x = min(df["orthoANI"].tolist())
max_x = max(df["orthoANI"].tolist())

axes[0][0].plot([min_x, max_x], [min_x, max_x], color='grey', linestyle='--')
axes[0][1].plot([min_x, max_x], [min_x, max_x], color='grey', linestyle='--')
axes[1][0].plot([min_x, max_x], [min_x, max_x], color='grey', linestyle='--')
axes[1][1].plot([min_x, max_x], [min_x, max_x], color='grey', linestyle='--')

axes[0][0].set_ylim(-0.02, 1.02)
axes[0][1].set_ylim(-0.02, 1.02)
axes[1][0].set_ylim(-0.02, 1.02)
axes[1][1].set_ylim(-0.02, 1.02)

# add legend
axes[0][0].legend()
axes[0][1].legend()
axes[1][0].legend()
axes[1][1].legend()


# Adjust the spacing between subplots
plt.tight_layout()

plt.savefig('combined_ani.pdf')
#plt.show()


'''
df2 = df[ df['orthoANI'] >= 0.0 ]

orthoANI = df2['orthoANI'].tolist()
pyani = df2['pyani'].tolist()
fmhANI = df2['avg_cANI_k21_sc10'].tolist()
fastani = df2['fastani'].tolist()

print(f'RMS of fmhANI and orthoANI: {rmse(orthoANI, fmhANI)}')
print(f'RMS of pyani and orthoANI: {rmse(orthoANI, pyani)}')
print(f'RMS of fastani and orthoANI: {rmse(fastani, orthoANI)}')
'''
