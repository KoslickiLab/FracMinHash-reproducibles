import subprocess
import numpy as np
from compare_two_genomes import compare_two_files_to_get_multiple_containments
from p_from_scaled_containment import compute_confidence_interval_one_step
import argparse
import sys

def parse_arguments(sys_args):
    parser = argparse.ArgumentParser()
    parser.add_argument("-k", "--ksize", type=int, default=21)
    parser.add_argument("-s", "--scale-factor", type=float, default=0.001)
    parser.add_argument("f1")
    parser.add_argument("f2")
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("-c", "--confidence", type=float, default=0.95)
    parser.add_argument("-N", "--num-runs", type=int, default=10)
    parser.add_argument("-p", "--mutation-rate", type=float, default=-1.0)
    parser.add_argument("--fout", help="will append to end of this file", default=None)
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_arguments(sys.argv)
    f1 = args.f1
    f2 = args.f2
    num_runs = int(args.num_runs)
    s = float(args.scale_factor)
    k = int(args.ksize)
    confidence = float(args.confidence)
    known_mutation_rate = float(args.mutation_rate)
    fout = args.fout
    
    size_1, size_2, size_union, size_intersection, true_containment, scaled_containments, sketch_sizes = compare_two_files_to_get_multiple_containments(f1, f2, k, s, num_runs)
    
    #print(scaled_containments)
    
    f = open("script.sh", 'w')
    for i in range(num_runs):
        seed_for_mash = i + 1
        command = "mash dist " + f1 + " " + f2 + " -s " +str(sketch_sizes[i])+ " -S " + str(seed_for_mash)
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
    #print(mash_containments)
    #print(scaled_containments)
    
    mash_c_avg = np.average(mash_containments)
    mash_c_var = np.var(mash_containments)
    scaled_c_avg = np.average(scaled_containments)
    scaled_c_var = np.var(scaled_containments)
    
    # get p from mash
    f = open('mash_distances', 'w')
    cmd = 'cut -f3 mash_output'
    cmd_args = cmd.split(' ')
    subprocess.call(cmd_args, stdout=f)
    f.close()
    
    mash_distances = []
    f = open('mash_distances', 'r')
    lines = f.readlines()
    for line in lines:
        d = float( line.strip() )
        mash_distances.append(d)
    f.close()
    
    #print(mash_distances)
    
    mash_distance = np.average(mash_distances)
    # get p from scaled containment
    scaled_containment = np.average(scaled_containments)
    L = int((size_1 + size_2)/2)
    conf_interval = compute_confidence_interval_one_step([scaled_containment], L, k, confidence, s)[0]
    # get true p from fast ani if == -1
    
    if fout is not None:
        sys.stdout = open(fout, 'a')
    print(known_mutation_rate, true_containment, mash_c_avg, mash_c_var, scaled_c_avg, scaled_c_var, mash_distance, conf_interval[6], conf_interval[4], conf_interval[5])
    if fout is not None:
        sys.stdout.close()
    