import screed
import mmh3
import subprocess
import string

__complementTranslation = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N", "R": "N", 'K': 'N'}
for char in string.ascii_uppercase:
    if char not in __complementTranslation.keys():
        __complementTranslation[char] = 'N'

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

def count_num_kmers_in_file(filename, k):
    kmer_hashes = set()
    for kmer in get_kmers_in_file(filename, k):
        kmer_hashes.add(get_hash_from_kmer(kmer))
    return len(kmer_hashes)

def get_true_mut_rate(filename1, filename2):
    cmd = "java -jar OAU.jar -u ./usearch11.0.667_i86linux32 --f1 " + filename1 + " -f2 " + filename2 + " ."
    args = cmd.split(' ')
    f = open('temp', 'w')
    subprocess.call(args, stdout=f)
    f.close()
    f = open('temp', 'r')
    true_ani = float(f.readlines()[-1].split('\t')[1])
    f.close()
    return (100.0-true_ani)/100.0
    

seed = 1
stats_filename = 'results'
k = 21
scale_factor = 0.01
num_runs = 20

f = open(stats_filename, 'w')
f.close()

#pairs = {('data/taxid_1159055_genomic.fna', 'data/taxid_210_10_genomic.fna'), ('data/taxid_197_350_genomic.fna', 'data/taxid_1246_2_genomic.fna'), ('data/taxid_210_345_genomic.fna', 'data/taxid_2044939_18_genomic.fna'), ('data/taxid_1926614_genomic.fna', 'data/taxid_1951748_genomic.fna'), ('data/taxid_197_129_genomic.fna', 'data/taxid_210_675_genomic.fna'), ('data/taxid_1385722_genomic.fna', 'data/taxid_1926614_genomic.fna'), ('data/taxid_197_350_genomic.fna', 'data/taxid_210_517_genomic.fna'), ('data/taxid_1303688_genomic.fna', 'data/taxid_197_341_genomic.fna'), ('data/taxid_2026734_61_genomic.fna', 'data/taxid_197_595_genomic.fna'), ('data/taxid_210_173_genomic.fna', 'data/taxid_1946093_genomic.fna'), ('data/taxid_195_531_genomic.fna', 'data/taxid_478547_genomic.fna'), ('data/taxid_1303688_genomic.fna', 'data/taxid_889213_genomic.fna'), ('data/taxid_887320_genomic.fna', 'data/taxid_197_501_genomic.fna'), ('data/taxid_197_736_genomic.fna', 'data/taxid_197_596_genomic.fna'), ('data/taxid_197_15_genomic.fna', 'data/taxid_195_543_genomic.fna'), ('data/taxid_1660076_genomic.fna', 'data/taxid_1899355_0_genomic.fna'), ('data/taxid_1946015_genomic.fna', 'data/taxid_2024894_186_genomic.fna'), ('data/taxid_197_830_genomic.fna', 'data/taxid_32022_51_genomic.fna'), ('data/taxid_197_170_genomic.fna', 'data/taxid_210_672_genomic.fna'), ('data/taxid_102617_0_genomic.fna', 'data/taxid_1442601_genomic.fna'), ('data/taxid_1385721_genomic.fna', 'data/taxid_197_408_genomic.fna'), ('data/taxid_2161947_genomic.fna', 'data/taxid_1903292_genomic.fna'), ('data/taxid_197_208_genomic.fna', 'data/taxid_2133953_genomic.fna'), ('data/taxid_197_871_genomic.fna', 'data/taxid_888826_genomic.fna'), ('data/taxid_210_266_genomic.fna', 'data/taxid_992077_genomic.fna'), ('data/taxid_866346_genomic.fna', 'data/taxid_1946153_genomic.fna')}
pairs = {('data/taxid_210_429_genomic.fna', 'data/taxid_992018_genomic.fna'), ('data/taxid_710393_3_genomic.fna', 'data/taxid_195_543_genomic.fna'), ('data/taxid_1300157_genomic.fna', 'data/taxid_2133953_genomic.fna'), ('data/taxid_210_283_genomic.fna', 'data/taxid_210_197_genomic.fna'), ('data/taxid_1952326_genomic.fna', 'data/taxid_210_449_genomic.fna'), ('data/taxid_210_238_genomic.fna', 'data/taxid_210_43_genomic.fna'), ('data/taxid_1156915_genomic.fna', 'data/taxid_2093811_4_genomic.fna'), ('data/taxid_197_573_genomic.fna', 'data/taxid_691338_genomic.fna'), ('data/taxid_197_316_genomic.fna', 'data/taxid_197_623_genomic.fna'), ('data/taxid_1951807_genomic.fna', 'data/taxid_1781242_genomic.fna'), ('data/taxid_148814_8_genomic.fna', 'data/taxid_123841_0_genomic.fna'), ('data/taxid_1345580_genomic.fna', 'data/taxid_210_634_genomic.fna'), ('data/taxid_2026734_61_genomic.fna', 'data/taxid_197_731_genomic.fna'), ('data/taxid_1840518_genomic.fna', 'data/taxid_210_149_genomic.fna'), ('data/taxid_195_684_genomic.fna', 'data/taxid_197_845_genomic.fna'), ('data/taxid_1262857_genomic.fna', 'data/taxid_730_17_genomic.fna'), ('data/taxid_1345524_genomic.fna', 'data/taxid_1856286_genomic.fna'), ('data/taxid_197_181_genomic.fna', 'data/taxid_887320_genomic.fna'), ('data/taxid_1262857_genomic.fna', 'data/taxid_1947718_genomic.fna'), ('data/taxid_197_826_genomic.fna', 'data/taxid_887320_genomic.fna'), ('data/taxid_210_680_genomic.fna', 'data/taxid_210_672_genomic.fna')}
for (filename1, filename2) in pairs:
    print('Analyzing pair:')
    print(filename1, filename2)
    mutation_rate = get_true_mut_rate(filename1, filename2)
    print('Mut rate: ' + str(mutation_rate))
    cmd = "python test_code.py " + filename1 + " " + filename2 + " -k " + str(k) + " -s " + str(scale_factor) + " --seed " + str(seed) + " -c 0.95 -N " + str(num_runs) + " -p " + str(mutation_rate) + " --fout " + stats_filename
    args = cmd.split(' ')
    subprocess.call(args)
    f = open(stats_filename, 'a')
    f.write(filename1 + ' ' + filename2 + '\n')
    f.close()
    

