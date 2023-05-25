import os
import sys
import argparse
import glob
import pprint

import numpy as np
import pandas as pd

from collections import defaultdict, namedtuple

from sourmash.logging import notify

pyani_res = namedtuple('pyani_res',
                           'comparison_name, identA, identB, pyani_ident, pyani_coverage, pyani_aln_length, pyani_sim_errors, pyani_hadamard')

def main(args):


     # read in list of comparisons we actually want
    compareInfo = pd.read_csv(args.comparison_csv)
    # make list of comparisons
    comparisons = list(zip(compareInfo.identA, compareInfo.identB))

    # find ANI files to pull from
    if "3" in str(args.pyani_version):
        len_fn= "matrix_aln_lengths_1.tab"
        cov_fn = "matrix_coverage_1.tab"
        had_fn = "matrix_hadamard_1.tab"
        id_fn = "matrix_identity_1.tab"
        se_fn = "matrix_sim_errors_1.tab"
    elif "2" in str(args.pyani_version):
        len_fn= "ANIb_alignment_lengths.tab"
        cov_fn = "ANIb_alignment_coverage.tab"
        had_fn = "ANIb_hadamard.tab"
        id_fn = "ANIb_percentage_identity.tab"
        se_fn = "ANIb_similarity_errors.tab"

    results_dir = args.pyani_results_dir
    lenF = os.path.join(results_dir, len_fn)
    covF = os.path.join(results_dir, cov_fn)
    hadamardF = os.path.join(results_dir, had_fn)
    identF = os.path.join(results_dir, id_fn)
    simerrF = os.path.join(results_dir, se_fn)

    # read in all matrices
    lenD = pd.read_csv(lenF, sep="\t", header=0, index_col=0)
    covD = pd.read_csv(covF, sep="\t", header=0, index_col=0)
    hadD = pd.read_csv(hadamardF, sep="\t", header=0, index_col=0)
    idD = pd.read_csv(identF, sep="\t", header=0, index_col=0)
    seD = pd.read_csv(simerrF, sep="\t", header=0, index_col=0)

    # use headers on one file to get full column names, in case they're not just the accessions:
    names = lenD.columns.tolist()

    # select comparison info from pyani results files
    results = []
    num_empty=0
    for n, (idA, idB) in enumerate(comparisons):
        fwd_name = f"{idA}_x_{idB}"
        rev_name = f"{idB}_x_{idA}"

        # now grab pyani values
        pyani_ident, pyani_coverage, pyani_aln_length, pyani_sim_errors, pyani_hadamard = np.nan, np.nan, np.nan, np.nan, np.nan

        try:
            a_label = [x for x in names if x.startswith(idA)][0]
            b_label = [x for x in names if x.startswith(idB)][0]
        except IndexError:
            num_empty+=1
            #print(f"skipping comparison {fwd_name}")
            continue # skip this comparison


        # matrix will not be symmetric. Average values for each direction.
        pyani_identA = idD.at[a_label, b_label]
        pyani_identB = idD.at[b_label, a_label]
        pyani_ident = np.mean([pyani_identA, pyani_identB])

        pyani_coverageA = covD.at[a_label, b_label]
        pyani_coverageB = covD.at[b_label, a_label]
        pyani_coverage = np.mean([pyani_coverageA, pyani_coverageB])

        pyani_aln_lengthA = lenD.at[a_label, b_label]
        pyani_aln_lengthB = lenD.at[b_label, a_label]
        pyani_aln_length = np.mean([pyani_aln_lengthA, pyani_aln_lengthB])

        pyani_sim_errorsA = seD.at[a_label, b_label]
        pyani_sim_errorsB = seD.at[b_label, a_label]
        pyani_sim_errors = np.mean([pyani_sim_errorsA, pyani_sim_errorsB])

        pyani_hadamardA = hadD.at[a_label, b_label]
        pyani_hadamardB = hadD.at[b_label, a_label]
        pyani_hadamard = np.mean([pyani_hadamardA, pyani_hadamardB])

        this_info = pyani_res(fwd_name, idA, idB, pyani_ident, pyani_coverage, pyani_aln_length, pyani_sim_errors, pyani_hadamard)
        results.append(this_info)

    notify(f"No pyANI results available for {num_empty} desired comparisons.")
    # convert path ANI comparison info to pandas dataframe. to do: just write this with csv dictwriter to save on dict conversion
    aniDF = pd.DataFrame.from_records(results, columns = pyani_res._fields)

    # print to csv
    aniDF.to_csv(args.output_csv, index=False)
    print(f"done! path comparison info written to {args.output_csv}")

def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("pyani_results_dir")
    p.add_argument("--pyani-version", default = 'v0.3')
    p.add_argument("-c", "--comparison-csv", required=True)
    p.add_argument("-o", "--output-csv", required=True)
    #p.add_argument("--taxonomy", default="gtdb-rs207.taxonomy.csv")
    #p.add_argument("--rank") # comparison rank
    #p.add_argument("--ranktax-name") # just do the one ranktax
    #p.add_argument("--labels")
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)

