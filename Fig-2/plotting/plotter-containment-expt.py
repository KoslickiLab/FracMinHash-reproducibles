import matplotlib as mpl
#mpl.rcParams['font.family'] = 'Lato'
#mpl.rcParams['font.sans-serif'] = 'Lato'
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.pyplot import figure

figure(figsize=(5, 4), dpi=80)
plt.style.use('seaborn-bright')

df = pd.read_csv('containment-results.txt', delimiter=',', header=None)
df.columns=['c', 'true_c', 'mash_c', 'mash_c_var', 'scaled_c', 'scaled_c_var']
df_old = df

n = len(df['true_c'].tolist())
s = [10 for i in range(n)]

delta_x = 0.005
linewidth = 0.8
opacity=0.8

plt.plot( df['true_c'].tolist(), df['true_c'].tolist(), alpha=0.8, linestyle='--', linewidth=1, color='grey' )

plt.scatter( df['true_c'].tolist(), df['mash_c'].tolist(), marker='o', s=s, alpha=0.8, label="MinHash Containment", color='red')
for (c_t, c, c_var) in zip(df['true_c'].tolist(), df['mash_c'].tolist(), df['mash_c_var'].tolist()):
    c_stddev = c_var ** 0.5
    plt.plot( [c_t, c_t], [c-c_stddev, c+c_stddev], color='red', linewidth=linewidth, alpha=opacity)
    plt.plot( [c_t-delta_x, c_t+delta_x], [c-c_stddev, c-c_stddev], color='red', linewidth=linewidth, alpha=opacity)
    plt.plot( [c_t-delta_x, c_t+delta_x], [c+c_stddev, c+c_stddev], color='red', linewidth=linewidth, alpha=opacity)


plt.scatter( df['true_c'].tolist(), df['scaled_c'].tolist(), marker='o', s=s, alpha=0.8, label="FMH Containment", color='blue')
for (c_t, c, c_var) in zip(df['true_c'].tolist(), df['scaled_c'].tolist(), df['scaled_c_var'].tolist()):
    c_stddev = c_var ** 0.5
    plt.plot( [c_t, c_t], [c-c_stddev, c+c_stddev], color='blue', linewidth=linewidth, alpha=opacity)
    plt.plot( [c_t-delta_x, c_t+delta_x], [c-c_stddev, c-c_stddev], color='blue', linewidth=linewidth, alpha=opacity)
    plt.plot( [c_t-delta_x, c_t+delta_x], [c+c_stddev, c+c_stddev], color='blue', linewidth=linewidth, alpha=opacity)


plt.grid(alpha=0.2)
plt.legend()
plt.xlabel("True containment index")
plt.ylabel("Computed containment index")
plt.savefig('compare-containment-all.pdf')
plt.close()
