import sys
import screed
import mmh3
import argparse
import string

__complementTranslation = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N", "R": "N"}
for char in string.ascii_uppercase:
    if char not in __complementTranslation.keys():
        __complementTranslation[char] = 'N'
"""
Simple implementation of scaled minhash
"""
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

def get_hash_from_kmer(kmer, seed):
	hash_value = mmh3.hash64(kmer, seed=seed)[0]
	if hash_value < 0:
		hash_value += 2**64
	return hash_value

def compare_two_files(filename_1, filename_2, k, scale_facor, seed):
	H = int(2**64)
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

	print("True containment: ")
	print(1.0*size_intersection/size_1, 1.0*size_intersection/size_2)
	print("Scaled containment: ")
	print(smh1.get_containment(smh2), smh2.get_containment(smh1))
	true_containment = 1.0*size_intersection/size_1
	scaled_containment = smh1.get_containment(smh2)
	sketch_size = smh1.get_sketch_size()
	return true_containment, scaled_containment, sketch_size


"""
will return: size_1, size_2, size_union, size_intersection, true_containment, list_of_scaled_containments
the size of the list = num_runs
"""
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

def parse_arguments(sys_args):
	parser = argparse.ArgumentParser()
	parser.add_argument("-k", "--ksize", type=int, default=21)
	parser.add_argument("-s", "--scale-factor", type=float, default=0.001)
	parser.add_argument("--f1", default=None)
	parser.add_argument("--f2", default=None)
	parser.add_argument("--seed", type=int, default=1)
	args = parser.parse_args()
	return args

if __name__=="__main__":
	args = parse_arguments(sys.argv)
	filename_1 = args.f1
	filename_2 = args.f2
	k = args.ksize
	scale_factor = args.scale_factor
	seed = args.seed
	
	values = compare_two_files(filename_1, filename_2, k, scale_factor, seed)
	print(values)
