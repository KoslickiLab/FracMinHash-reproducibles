import matplotlib as mpl
#mpl.rcParams['font.family'] = 'Lato'
#mpl.rcParams['font.sans-serif'] = 'Lato'
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.pyplot import figure

df = pd.read_csv('results', delimiter=' ', header=None)
df.columns=['mut_rate', 'true_c', 'mash_c', 'mash_c_var', 'scaled_c', 'scaled_c_var', 'mash_dist', 'p_est', 'p_low', 'p_high', 'p_est_j']
df.sort_values('mut_rate')

n = len(df['true_c'].tolist())
s = [24 for i in range(n)]


figure(figsize=(5, 5), dpi=80)
plt.style.use('seaborn-bright')

plt.plot( df['mut_rate'].tolist(), df['mut_rate'].tolist(), alpha=0.8, linestyle='-', linewidth=1, color='grey' )
plt.scatter( df['mut_rate'].tolist(), df['mash_dist'].tolist(), marker='.', s=s, alpha=0.6, label="Mash Distance", color='red')
plt.scatter( df['mut_rate'].tolist(), df['p_est'].tolist(), marker='.', s=s, alpha=0.6, label="FracMinHash Distance (Containment)", color='blue')
plt.scatter( df['mut_rate'].tolist(), df['p_est'].tolist(), marker='<', s=s, alpha=0.6, label="FracMinHash Distance (Jaccard)", color='purple')

delta_x = 0.003
linewidth = 0.9
opacity=0.6
for (p, ci_high, ci_low) in list( zip( df['mut_rate'].tolist(), df['p_low'].tolist(), df['p_high'].tolist() ) ):
    plt.plot( [p, p], [ci_low, ci_high], color='blue', linewidth=linewidth, alpha=opacity)
    plt.plot( [p-delta_x, p+delta_x], [ci_low, ci_low], color='blue', linewidth=linewidth, alpha=opacity)
    plt.plot( [p-delta_x, p+delta_x], [ci_high, ci_high], color='blue', linewidth=linewidth, alpha=opacity)

plt.grid(alpha=0.2)
plt.legend()
plt.xlabel("True distance ($1.0 -$ANI)")
plt.ylabel("Predicted distance")
plt.savefig('compare-ani-real-genomes.pdf')
