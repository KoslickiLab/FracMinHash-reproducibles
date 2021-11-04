import subprocess
from mutate_genome import mutate_file

mut_rates = [0.001, 0.01, 0.025, 0.05, 0.07, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6]
seed = 2
stats_filename = 'stats'
reference_genome = 'ecoli.fasta'
mutated_filename_prefix = 'ecoli_mutated_'
k = 21
scale_factor = 0.01
num_runs = 10

for mutation_rate in mut_rates:
    mutate_file(reference_genome, mutated_filename_prefix+str(mutation_rate*100)+'.fasta', mutation_rate, seed)

f = open(stats_filename, 'w')
f.close()

for mutation_rate in mut_rates:
    mutated_filename = mutated_filename_prefix+str((mutation_rate*100))+'.fasta'
    cmd = "python test_code.py ecoli.fasta " + mutated_filename + " -k " + str(k) + " -s " + str(scale_factor) + " --seed " + str(seed) + " -c 0.95 -N " + str(num_runs) + " -p " + str(mutation_rate) + " --fout " + stats_filename
    args = cmd.split(' ')
    subprocess.call(args)