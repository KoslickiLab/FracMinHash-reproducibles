import screed
import random

bases = ['A', 'C', 'G', 'T']

def get_names_and_sequences_in_file(filename):
	with screed.open(filename) as f:
		for record in f:
			yield record.name, record.sequence
            
def mutate_sequence(sequence, mutation_rate, seed):
    random.seed(seed)
    sequence = sequence.upper()
    mutated_seq = []
    for i in range(len(sequence)):
        if sequence[i] != 'N':
            r = random.random()
            if r <= mutation_rate:
                base = bases[ random.randint(0,3) ]
                while base == sequence[i]:
                    base = bases[ random.randint(0,3) ]
                mutated_seq.append(base)
            else:
                mutated_seq.append(sequence[i])
    return "".join(mutated_seq)

def write_fasta(filename, names, sequences):
    f = open(filename, 'w')
    for (name, seq) in list( zip(names, sequences) ):
        f.write('> ' + name + '\n')
        f.write(seq + '\n')
    f.close()
    
def mutate_file(filename, output_filename, mutation_rate, seed):
    names = []
    sequences = []
    for name, seq in get_names_and_sequences_in_file(filename):
        names.append(name)
        mutated_sequence = mutate_sequence(seq, mutation_rate, seed)
        sequences.append(mutated_sequence)
    write_fasta(output_filename, names, sequences)
    
mut_rates = [0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]
seed = 1
for mutation_rate in mut_rates:
    mutate_file('ecoli.fasta', 'ecoli_mutated_'+str(int(mutation_rate*100))+'.fasta', mutation_rate, seed)