import matplotlib as mpl
#mpl.rcParams['font.family'] = 'Lato'
#mpl.rcParams['font.sans-serif'] = 'Lato'
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.pyplot import figure

figure(figsize=(8,8), dpi=80)
plt.style.use('seaborn-bright')

df = pd.read_csv('results', delimiter=' ', header=None)
df.columns=['mut_rate', 'true_c', 'mash_c', 'mash_c_var', 'scaled_c', 'scaled_c_var', 'mash_dist', 'p_est', 'p_low', 'p_high']

df.sort_values('mut_rate')


n = len(df['true_c'].tolist())
s = [12 for i in range(n)]


figure(figsize=(8,8), dpi=80)
plt.style.use('seaborn-bright')

plt.scatter( df['mut_rate'].tolist(), df['mash_dist'].tolist(), marker='o', s=s, alpha=0.6, label="Mash Distance")
plt.scatter( df['mut_rate'].tolist(), df['p_est'].tolist(), marker='o', s=s, alpha=0.6, label="Theorem 9 Distance", color='red')
plt.plot( [min(df['mut_rate'].tolist()), max(df['mut_rate'].tolist())] , [min(df['mut_rate'].tolist()), max(df['mut_rate'].tolist())], alpha=0.8, linestyle='--', linewidth=1, color='grey' )

delta_x = 0.005
linewidth = 0.9
opacity=0.5
for (p, ci_high, ci_low) in list( zip( df['mut_rate'].tolist(), df['p_low'].tolist(), df['p_high'].tolist() ) ):
    plt.plot( [p, p], [ci_low, ci_high], color='red', linewidth=linewidth, alpha=opacity)
    plt.plot( [p-delta_x, p+delta_x], [ci_low, ci_low], color='red', linewidth=linewidth, alpha=opacity)
    plt.plot( [p-delta_x, p+delta_x], [ci_high, ci_high], color='red', linewidth=linewidth, alpha=opacity)
    
    

plt.grid(alpha=0.2)
plt.legend()
plt.xlabel("True distance ($1.0 -$ANI)")
plt.ylabel("Predicted distance")
plt.savefig('compare-ani-200k-genomes.pdf')