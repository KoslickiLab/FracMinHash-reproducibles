import screed
import subprocess
import mmh3
import numpy as np

__complementTranslation = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N", "R": "N"}

class ScaledMinHash:
    def __init__(self, scale_factor, max_hash_value):
        self.hash_set = set()
        self.H = max_hash_value
        self.scale_factor = scale_factor
        self.raw_elements = set()
        
    def add_value(self, hash_value):
        if hash_value <= self.H * self.scale_factor:
            self.hash_set.add(hash_value)
        self.raw_elements.add(hash_value)
            
    def add_values(self, hash_values):
        for hash_value in hash_values:
            self.add_value(hash_value)
            
    def remove(self, hash_value):
        self.hash_set -= hash_value
        
    def print_hash_set(self):
        print(self.H, self.scale_factor, self.hash_set)
        
    def get_containment(self, smh):
        return 1.0 * len(self.hash_set.intersection(smh.hash_set)) / len(self.hash_set)
    
    def get_scaled_containment(self, smh):
        bf = 1 - (1 - self.scale_factor) ** len(self.raw_elements)
        return 1.0 * len(self.hash_set.intersection(smh.hash_set)) / ( len(self.hash_set) * bf )

    def get_sketch_size(self):
        return len( self.hash_set )

def reverse_complement(s):
    """
    Return reverse complement of 's'.
    """
    c = "".join(reversed([__complementTranslation[n] for n in s]))
    return c

def canonical_kmers(seq, k):
    for start in range(len(seq) - k + 1):
        kmer = seq[start: start + k].upper()
        rev_kmer = reverse_complement(kmer)
        if rev_kmer < kmer:
            kmer = rev_kmer

        if any(c not in "ACGT" for c in kmer):
            continue
        yield kmer

def get_kmers_in_file(filename, k):
	with screed.open(filename) as f:
		for record in f:
			for kmer in canonical_kmers(record.sequence, k):
				yield kmer

def get_hash_from_kmer(kmer, seed=0):
	hash_value = mmh3.hash64(kmer, seed=seed)[0]
	if hash_value < 0:
		hash_value += 2**64
	return hash_value

# c is float, 0 < c < 1
def extract_part_of_genome(c, genome_filename, out_filename):
    subprocess.call(['rm', out_filename])
    length = 0
    with screed.open(genome_filename) as f:
        for record in f:
            length += len(record.sequence)
    required_length = int(length * c)
    with screed.open(genome_filename) as f:
        for record in f:
            if len(record.sequence) > required_length:
                small_str = record.sequence[:required_length]
                break
    f2 = open(out_filename, 'w')
    f2.write('>small_seq\n')
    f2.write(small_str)
    f2.write('\n')
    f2.close()
    
def create_super_metagenome(metagenome_filename, small_genome_filename, super_mg_filename):
    args = ['cat', metagenome_filename, small_genome_filename]
    f = open(super_mg_filename, 'w')
    subprocess.call(args, stdout=f)
    f.close()

def count_num_kmers_in_file(filename, k):
    kmer_hashes = set()
    for kmer in get_kmers_in_file(filename, k):
        kmer_hashes.add(get_hash_from_kmer(kmer))
    return len(kmer_hashes)

def get_mash_containments(f1, f2, sketch_size, size_union, size_1):
    f = open("script.sh", 'w')
    for i in range(num_runs):
        seed_for_mash = i + 1
        command = "mash dist " + f1 + " " + f2 + " -s " +str(sketch_size)+ " -S " + str(seed_for_mash)
        f.write(command)
        f.write("\n")
    f.close()
    
    f = open('mash_output', 'w')
    cmd = "bash script.sh"
    cmd_args = cmd.split(' ')
    subprocess.call(cmd_args, stdout=f)
    f.close()
    
    f = open('mash_jaccards', 'w')
    cmd = 'cut -f5 mash_output'
    cmd_args = cmd.split(' ')
    subprocess.call(cmd_args, stdout=f)
    f.close()
    
    mash_jaccards = []
    f = open('mash_jaccards', 'r')
    lines = f.readlines()
    for line in lines:
        v1 = float(line.split('/')[0])
        v2 = float(line.split('/')[1])
        mash_jaccards.append( 1.0 * v1 / v2 )
    #print(mash_jaccards)
    f.close()
    
    mash_containments = []
    for j in mash_jaccards:
        c = j * 1.0 * size_union / size_1
        mash_containments.append(c)
    return mash_containments
    

def compare_two_files_to_get_multiple_containments(filename_1, filename_2, k, scale_facor, num_runs):
    seeds = [i+1 for i in range(num_runs)]
    H = int(2**64)

    sketch_sizes = []
    scaled_containments = []
    for seed in seeds:
        kmer_hashes_1 = set()
        kmer_hashes_2 = set()
        for kmer in get_kmers_in_file(filename_1, k):
            kmer_hashes_1.add(get_hash_from_kmer(kmer, seed=seed))
        for kmer in get_kmers_in_file(filename_2, k):
            kmer_hashes_2.add(get_hash_from_kmer(kmer, seed=seed))
            
        size_1 = len(kmer_hashes_1)
        size_2 = len(kmer_hashes_2)
        size_union = len( kmer_hashes_1.union( kmer_hashes_2 ) )
        size_intersection = len( kmer_hashes_1.intersection( kmer_hashes_2 ) )            
        
        smh1 = ScaledMinHash(scale_facor, H)
        smh1.add_values(kmer_hashes_1)
        smh2 = ScaledMinHash(scale_facor, H)
        smh2.add_values(kmer_hashes_2)
        
        scaled_containment = smh1.get_containment(smh2)
        sketch_size = smh1.get_sketch_size()
        
        sketch_sizes.append(sketch_size)
        scaled_containments.append(scaled_containment)
    
    true_containment = 1.0*size_intersection/size_1
    return size_1, size_2, size_union, size_intersection, true_containment, scaled_containments, sketch_sizes
    

g_filename = 'ecoli.fasta'
smallg_filename = 'temp.fasta'
smg_filename = 'supermg.fasta'
mg_filename = 'SRR492190.contigs.fa'

scale_factor = 0.0005
k = 21
C = 0.1
num_runs = 20
seed = 1
num_hashes = 2000

num_kmers = count_num_kmers_in_file(g_filename, k)
scale_factor = 1.0*num_hashes/num_kmers

for C in [0.99, 0.95, 0.90]:
    extract_part_of_genome(C, g_filename, smallg_filename)
    create_super_metagenome(mg_filename, smallg_filename, smg_filename)
    expected_sketch_size = int(num_kmers * scale_factor)
    size_1, size_2, size_union, size_intersection, true_containment, scaled_containments, sketch_sizes = compare_two_files_to_get_multiple_containments(g_filename, smg_filename, k, scale_factor, num_runs)
    mash_containments = get_mash_containments(g_filename, smg_filename, expected_sketch_size, size_union, size_1)
    #print(size_1, size_2, size_union, size_intersection)
    #print(C, scale_factor)
    sc_c_avg = np.average(scaled_containments)
    sc_c_var = np.var(scaled_containments)
    mash_c_avg = np.average(mash_containments)
    mash_c_var = np.var(mash_containments)
    print(C, scale_factor, true_containment, sc_c_avg, sc_c_var, mash_c_avg, mash_c_var)

