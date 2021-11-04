import matplotlib as mpl
#mpl.rcParams['font.family'] = 'Lato'
#mpl.rcParams['font.sans-serif'] = 'Lato'
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.pyplot import figure

figure(figsize=(5,5), dpi=80)
plt.style.use('seaborn-bright')

df = pd.read_csv('stats', delimiter=' ', header=None)
df.columns=['mut_rate', 'true_c', 'mash_c', 'mash_c_var', 'scaled_c', 'scaled_c_var', 'mash_dist', 'p_est', 'p_low', 'p_high']

n = len(df['true_c'].tolist())
s = [12 for i in range(n)]


plt.scatter( df['true_c'].tolist(), df['mash_c'].tolist(), marker='^', s=s, alpha=0.5, label="Mash Containment")
plt.scatter( df['true_c'].tolist(), df['scaled_c'].tolist(), marker='o', s=s, alpha=0.5, label="Scaled Containment")
plt.plot( [ 1.0*i/20 for i in range(21) ], [ 1.0*i/20 for i in range(21) ], alpha=0.4, linestyle='--', linewidth=1, color='grey' )
plt.grid(alpha=0.2)
plt.legend()
plt.xlabel("True containment index")
plt.ylabel("Predicted containment index")
plt.savefig('compare-containemnt-simulated.pdf')
plt.close()


figure(figsize=(5,5), dpi=80)
plt.style.use('seaborn-bright')

plt.scatter( df['mut_rate'].tolist(), df['mash_dist'].tolist(), marker='o', s=s, alpha=0.8, label="Mash Distance")
plt.scatter( df['mut_rate'].tolist(), df['p_est'].tolist(), marker='o', s=s, alpha=0.8, label="Theorem 9 Distance", color='red')
plt.plot( [ 1.0*i/20 for i in range(11) ], [ 1.0*i/20 for i in range(11) ], alpha=0.8, linestyle='--', linewidth=1, color='grey' )

delta_x = 0.01
linewidth = 0.7
opacity=0.8
for (p, ci_high, ci_low) in list( zip( df['mut_rate'].tolist(), df['p_low'].tolist(), df['p_high'].tolist() ) ):
    plt.plot( [p, p], [ci_low, ci_high], color='red', linewidth=linewidth, alpha=opacity)
    plt.plot( [p-delta_x, p+delta_x], [ci_low, ci_low], color='red', linewidth=linewidth, alpha=opacity)
    plt.plot( [p-delta_x, p+delta_x], [ci_high, ci_high], color='red', linewidth=linewidth, alpha=opacity)
    
    

plt.grid(alpha=0.2)
plt.legend()
plt.xlabel("True distance ($1.0 -$ANI)")
plt.ylabel("Predicted distance")
plt.savefig('compare-ani-staphylococcus.pdf')