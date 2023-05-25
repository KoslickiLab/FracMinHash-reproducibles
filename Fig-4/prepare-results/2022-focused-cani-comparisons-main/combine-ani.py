import os
import sys
import argparse
import pandas as pd
from functools import partial, reduce

from sourmash.tax import tax_utils
from sourmash.lca import lca_utils
from sourmash.logging import notify



def get_lineage(name, tax_assign):
    # handle GCA/GCF issues
    ident = tax_utils.get_ident(name) #, keep_identifier_versions=True)
    try:
        lineage = tax_assign[ident]
    except KeyError:
        if "GCF" in ident:
            ident = ident.replace("GCF", "GCA")
        elif "GCA" in ident:
            ident = ident.replace("GCA", "GCF")
        lineage = tax_assign.get(ident, None)
    return lineage

def find_lca(row, tax_assign):
    linA = get_lineage(row['identA'], tax_assign)
    linB = get_lineage(row['identB'], tax_assign)
    lintree = lca_utils.build_tree([linA, linB])
    lca, node_len = lca_utils.find_lca(lintree)
    if lca is not None:
        row['lca_lineage'] = lca_utils.display_lineage(lca)
        row['lca_rank'] = lca[-1].rank
    return row

def main(args):

    pyani_res = pd.read_csv(args.pyani_csv)
    fastani_res = pd.read_csv(args.fastani_csv)
    mash_res = pd.read_csv(args.mash_csv)
    orthoani_res = pd.read_csv(args.orthoani_csv)
    sourmash_res = [pd.read_csv(f) for f in args.sourmash_csv]

    ani_dfs = [pyani_res, fastani_res, mash_res, orthoani_res] + sourmash_res

    #id_cols = ["comparison_name", "identA", "identB"]
    outer_merge = partial(pd.merge, how='outer')
    aniDF = reduce(outer_merge, ani_dfs)

    if args.taxonomy:
        # could add LCA info here if we want it.
        tax_assign = tax_utils.MultiLineageDB.load([args.taxonomy],
                                                   keep_full_identifiers=False,
                                                   keep_identifier_versions=False)
        aniDF = aniDF.apply(find_lca, axis=1, tax_assign=tax_assign)

    # print to csv
    aniDF.to_csv(args.output_csv, index=False)
    print(f"done! Combined comparison info written to {args.output_csv}")

def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("-p", "--pyani-csv")
    p.add_argument("-f", "--fastani-csv")
    p.add_argument("-m", "--mash-csv")
    p.add_argument("--orthoani-csv")
    p.add_argument("-s", "--sourmash-csv", nargs='+')
    p.add_argument("-o", "--output-csv", required=True)
    p.add_argument("-t", "--taxonomy", default="/group/ctbrowngrp/sourmash-db/gtdb-rs207/gtdb-rs207.taxonomy.csv")
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)

