## python api comparisons

import os
import sys
import argparse
import glob
import pprint

import pandas as pd

import sourmash
from collections import namedtuple
from sourmash.sourmash_args import load_file_as_signatures
from sourmash.sketchcomparison import FracMinHashComparison
from sourmash.picklist import SignaturePicklist
from sourmash.logging import notify



def main(args):
    ksize=args.ksize
    scaled=args.scaled
    ani_thresh=args.ani_threshold

    CompareResult = namedtuple('CompareResult', ['comparison_name', 'identA', 'identB',
                                                 'ksize', 'scaled', 'avg_cANI', 'jaccard',
                                                 'idA_in_idB', 'idA_in_idB_cANI',
                                                 'idA_in_idB_cANI_low', 'idA_in_idB_cANI_high',
                                                 'idB_in_idA', 'idB_in_idA_cANI',
                                                 'idB_in_idA_cANI_low', 'idB_in_idA_cANI_high',
                                                 'idA_hashes', 'idB_hashes', 'num_common'])

    ANIResult = namedtuple('ANIResult', ['comparison_name', 'identA', 'identB',
                                                 f'avg_cANI_k{ksize}_sc{scaled}'])

    db = sourmash.load_file_as_index(args.db)

    compareInfo = pd.read_csv(args.comparison_csv)
    # make list of comparisons
    comparisons = list(zip(compareInfo.identA, compareInfo.identB))

    # loop through comparisons
    results,ani_results = [], []
    for n, (idA, idB) in enumerate(comparisons):
        if idA ==idB:
            raise ValueError(f"Cannot compare same ident, {idA}")
        if n !=0 and n % 50 == 0:
            notify(f"... assessing {n}th comparison\n")

        comparison_name = f"{idA}_x_{idB}"
        picklist = SignaturePicklist('ident')
        picklist.init([idA, idB])
        sigs = db.select(picklist=picklist, ksize=args.ksize)
        # select and load sigA
        siglist = list(sigs.signatures())
        ss1 = siglist[0]
        # select and load sigB
        ss2 = siglist[1]
        print(ss1.name, ss2.name)

        cmp = FracMinHashComparison(ss1.minhash, ss2.minhash, cmp_scaled=args.scaled, estimate_ani_ci = True)
        cmp.estimate_all_containment_ani()
        idA_sc_hashes = len(cmp.mh1_cmp)
        idB_sc_hashes = len(cmp.mh2_cmp)
        contain1 = cmp.mh1_containment_in_mh2
        cANI_1 = cmp.ani_from_mh2_containment_in_mh1
        cANI_1_low = cmp.ani_from_mh1_containment_in_mh2_low
        cANI_1_high = cmp.ani_from_mh1_containment_in_mh2_high

        contain2 = cmp.mh2_containment_in_mh1
        cANI_2 = cmp.ani_from_mh2_containment_in_mh1
        cANI_2_low = cmp.ani_from_mh2_containment_in_mh1_low
        cANI_2_high = cmp.ani_from_mh2_containment_in_mh1_high

        res = CompareResult(comparison_name, idA, idB, ksize, scaled, \
                             cmp.avg_containment_ani, cmp.jaccard, contain1, \
                             cANI_1, cANI_1_low, cANI_1_high,
                             contain2, cANI_2, cANI_2_low, cANI_2_high,
                             idA_sc_hashes, idB_sc_hashes, len(cmp.intersect_mh))
        results.append(res)
        if cmp.avg_containment_ani >= ani_thresh:
            ani_res = ANIResult(comparison_name, idA, idB, cmp.avg_containment_ani)
            ani_results.append(ani_res)

    # convert comparison info to pandas dataframe
    comparisonDF = pd.DataFrame.from_records(results, columns = CompareResult._fields)
    aniDF= pd.DataFrame.from_records(ani_results, columns = ANIResult._fields)


    # print to csv
    comparisonDF.to_csv(args.output_csv, index=False)
    aniDF.to_csv(args.ani_csv, index=False)
    print(f"done! Detailed comparison info written to {args.output_csv}. Abbreviated ANI info written to {args.ani_csv}")

def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("db")
    p.add_argument('-c', "--comparison-csv")
    p.add_argument('-k', "--ksize", default=21, type=int)
    p.add_argument('-s', "--scaled", default=1000, type=int)
    p.add_argument("-t", "--ani-threshold", default=0.0, type=float)
    p.add_argument("-a", "--ani-csv", required=True)
    p.add_argument("-o", "--output-csv", required=True)
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)
