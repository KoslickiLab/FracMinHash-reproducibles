import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import linregress
import matplotlib as mpl
from scipy.stats.stats import pearsonr
import numpy as np
from sklearn.metrics import mean_squared_error

mpl.rcParams['font.family'] = 'Arial'

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

plt.grid(alpha=0.2)

# create the scatter plot
sns.scatterplot(data=df, x="orthoANI", y="fastani", alpha=0.5, label="fastani", color=colors['fastani'])
sns.scatterplot(data=df, x="orthoANI", y="pyani", alpha=0.5, label="pyani", color=colors['pyani'])
sns.scatterplot(data=df, x="orthoANI", y="avg_cANI_k21_sc10", alpha=0.5, label="FracMinHash", color=colors['avg_cANI_k21_sc10'])
sns.scatterplot(data=df, x="orthoANI", y="MashANI", alpha=0.5, label="Mash", color=colors['MashANI'])


# add title and labels
plt.title("OrthoANI vs ANI")
plt.xlabel("OrthoANI")
plt.ylabel("ANI using different tools")

# add linear regression line
for col in ["pyani", "avg_cANI_k21_sc10"]:
    slope, intercept, r_value, p_value, std_err = linregress(df["orthoANI"], df[col])
    plt.plot(df["orthoANI"], intercept + slope * df["orthoANI"], color=colors[col], alpha=0.7)

min_x = min(df["orthoANI"].tolist())
max_x = max(df["orthoANI"].tolist())

plt.plot([min_x, max_x], [min_x, max_x], color='grey', linestyle='--')

# add legend
plt.legend()
plt.savefig('combined_ani.pdf')

df2 = df[ df['orthoANI'] >= 0.0 ]

orthoANI = df2['orthoANI'].tolist()
pyani = df2['pyani'].tolist()
fmhANI = df2['avg_cANI_k21_sc10'].tolist()
fastani = df2['fastani'].tolist()

print(f'RMS of fmhANI and orthoANI: {rmse(orthoANI, fmhANI)}')
print(f'RMS of pyani and orthoANI: {rmse(orthoANI, pyani)}')
print(f'RMS of fastani and orthoANI: {rmse(fastani, orthoANI)}')
