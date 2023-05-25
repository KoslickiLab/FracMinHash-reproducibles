import os
import sys
import argparse
#import glob
#import pprint

import numpy as np
import pandas as pd

from collections import defaultdict, namedtuple
#from itertools import product

# sourmash for tax/lca utils
#from sourmash.tax import tax_utils
#from sourmash.lca import lca_utils
from sourmash.logging import notify

fastani_res = namedtuple('fastani_res',
                           'comparison_name, identA, identB, avg_fastani_ident, avg_fastani_alignment_fraction')


def main(args):

    # read in and format fastani results
    fastani = pd.read_csv(args.fastani, sep = "\t", header=None, names=['a_fn','b_fn','fastani_ident','count_bidirectional_frag_mappings','total_query_frags'])
    fastani["idA"] = fastani["a_fn"].str.rsplit(pat="/", n=1, expand=True)[1].str.rsplit(pat="_genomic.fna", n=1, expand=True)[0]
    fastani["idB"] = fastani["b_fn"].str.rsplit(pat="/", n=1, expand=True)[1].str.rsplit(pat="_genomic.fna", n=1, expand=True)[0]
    fastani["comparison_name"] = fastani["idA"] + "_x_" + fastani["idB"]
    fastani.drop_duplicates(inplace=True)
    fastani.set_index("comparison_name",inplace=True)
    fastani_comparisons = fastani.index.to_list()

    # read in list of comparisons we actually want
    compareInfo = pd.read_csv(args.comparison_csv)
    # make list of comparisons
    comparisons = list(zip(compareInfo.identA, compareInfo.identB))

    # select comparison info from fastani dataframe
    results = []
    num_empty=0
    for n, (idA, idB) in enumerate(comparisons):
        fwd_name = f"{idA}_x_{idB}"
        rev_name = f"{idB}_x_{idA}"
        # check if both comparison directions are available (>= 80% ani)
        if fwd_name in fastani_comparisons and rev_name in fastani_comparisons:
            # now grab directional fastani values. calc average ANI
            fwd_ani = fastani.at[fwd_name, "fastani_ident"]
            rev_ani = fastani.at[rev_name, "fastani_ident"]
            avg_ani = np.mean([fwd_ani, rev_ani])

            fwd_fmap = fastani.at[fwd_name, "count_bidirectional_frag_mappings"]
            rev_fmap = fastani.at[rev_name, "count_bidirectional_frag_mappings"]

            fwd_qfrag = fastani.at[fwd_name, "total_query_frags"]
            rev_qfrag = fastani.at[rev_name, "total_query_frags"]

            # get directional alignment fractions
            # alignment fraction = count_bidirectional_frag_mappings/total_query_frags
            fwd_af = float(fwd_fmap)/float(fwd_qfrag)
            rev_af = float(rev_fmap)/float(rev_qfrag)
            avg_af = np.mean([fwd_af, rev_af])
            this_info = fastani_res(fwd_name, idA, idB, avg_ani, avg_af)
            results.append(this_info)
        else:
            if fwd_name in fastani_comparisons:
                notify(f"only {fwd_name} directional comparison available")

            elif rev_name in fastani_comparisons:
                notify(f"only {rev_name} directional comparison available")
            else:
                num_empty +=1

    notify(f"No FastANI results available for {num_empty} desired comparisons.")

    # convert path ANI comparison info to pandas dataframe. to do: just write this with csv dictwriter to save on dict conversion
    aniDF = pd.DataFrame.from_records(results, columns = fastani_res._fields)

    # print to csv
    aniDF.to_csv(args.output_csv, index=False)
    print(f"done! Comparison info written to {args.output_csv}")

def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("-f", "--fastani", required=True)
    p.add_argument("-c", "--comparison-csv", required=True)
    p.add_argument("-o", "--output-csv", required=True)
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)

