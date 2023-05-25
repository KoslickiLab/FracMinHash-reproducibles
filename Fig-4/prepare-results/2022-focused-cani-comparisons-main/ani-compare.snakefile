import os
import csv
import pandas as pd

# compare different ANI methods on a series of genomes

out_dir = "output.ani-commonsp10-evolpath"
logs_dir = os.path.join(out_dir, "logs")

comparison_file = "gtdb-rs207.common-sp10-evolpaths.csv"
comparisons = pd.read_csv(comparison_file)
basename= "comparisons"
identA= comparisons['identA'].tolist()
identB= comparisons['identB'].tolist()
IDENTS = list(dict.fromkeys(identA + identB)) # rm duplicates but preserve order to make sure copy_genomes rule doesn't get rerun
PAIRS = zip(identA, identB)
COMPARISONS = [f"{a}_x_{b}" for a,b in PAIRS]
KSIZE = [21]
#KSIZE = [31] # ran separately, not included in combinedANI.csv since other paper comparisons are k21
SCALED = [1, 10, 100, 1000]

# GTDB genomes are here:
original_genomedir = "/home/ntpierce/2021-rank-compare/genbank/genomes"

rule all:
    input: 
        os.path.join(out_dir, "combinedANI.csv"),

### split into genome folders + unzip fna files and generate classes/labels, then run pyANI index and compare
# make a folder with all genomes
rule copy_and_unzip_genomes:
    input:
        expand(os.path.join(original_genomedir, "{acc}_genomic.fna.gz"),acc=IDENTS)
    output:
        classes=os.path.join(out_dir, "genomes", "py.classes.txt"),
        labels=os.path.join(out_dir, "genomes", "py.labels.txt"),
        genomes=expand(os.path.join(out_dir, "genomes", "{acc}_genomic.fna"), acc=IDENTS),
    params:
        acc_list = IDENTS,
        genome_dir = os.path.join(out_dir, 'genomes'),
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        time=1200,
        partition="bmh",
    run:
        import hashlib
        os.makedirs(params.genome_dir, exist_ok=True)
        with open(output.classes, 'w') as out_classes, open(output.labels, 'w') as out_labels:
            for acc in params.acc_list:
                fn = os.path.join(original_genomedir, f"{acc}_genomic.fna.gz")
                #dest = os.path.join(params.genome_dir, f"{acc}_genomic.fna.gz")
                dest_unz = os.path.join(params.genome_dir, f"{acc}_genomic.fna")
                md5_file = os.path.join(params.genome_dir, f"{acc}_genomic.md5")
                #shell("cp {fn} {dest}") # do we need both gzipped and unzipped versions??
                shell("gunzip -c {fn} > {dest_unz}")
                # get md5 of unzipped fna
                with open(dest_unz, "rb") as f:
                    bytes = f.read()
                    md5 = hashlib.md5(bytes).hexdigest()
                # write to md5 file
                with open(md5_file, 'w') as mfile:
                    mfile.write(f"{md5}\t{dest_unz}\n")
                fna_base = os.path.basename(dest_unz).rsplit('.fna')[0]
                out_classes.write(f"{md5}\t{fna_base}\t{acc}\n")
                out_labels.write(f"{md5}\t{fna_base}\t{acc}\n")

# we copy all the fasta files to the genome_dir for pyani, so we can just use those paths here
rule build_genome_info_files:
    input:
        expand(os.path.join(original_genomedir, "{acc}_genomic.fna.gz"),acc=IDENTS),
    output:
        filepaths=os.path.join(out_dir, "genomes", "genome-filepaths.txt"),
        fromfile_csv=os.path.join(out_dir, "genomes", "fromfile.csv")
    params:
        genome_dir=os.path.join(out_dir, "genomes"),
        original_genomedir = original_genomedir,
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        time=1200,
        partition="bmh",
    run:
        with open(str(output.filepaths), "w") as out, open(str(output.fromfile_csv), 'w') as ff_out:
            acc_list = IDENTS
            ff_out.write('name,genome_filename,protein_filename\n')
            for acc in acc_list:
                new_fn = os.path.join(params.genome_dir, f"{acc}_genomic.fna")
                out.write(f"{new_fn}\n")
                orig_fn = os.path.join(params.original_genomedir, f"{acc}_genomic.fna.gz")
                ff_out.write(f"{acc},{orig_fn},\n")

# sourmash
rule sourmash_sketch:
    input: os.path.join(out_dir, "genomes", "fromfile.csv")
    output: os.path.join(out_dir, "sourmash", "signatures.zip")
    conda: "conf/env/sourmash.yml"
    params:
        scaled=min(SCALED),
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        time=240,
        partition="high2", #"low2",
    log: os.path.join(logs_dir, "sourmash", "sketch.log")
    benchmark: os.path.join(logs_dir, "sourmash", "sketch.benchmark")
    shell:
        """
        sourmash sketch fromfile {input} -p dna,k=21,k=31,k=51,scaled={params.scaled} -o {output} > {log}
        """

# just use python api -- not very many comparisons and we want to keep some extra info (nhashes, ncommon, etc)
rule sourmash_api_compare:
    input: 
        sigs=os.path.join(out_dir, "sourmash", "signatures.zip"),
        comparisons=comparison_file,#COMPARISON_FILES
    output: 
        compare_csv=os.path.join(out_dir, "sourmash", "{basename}.k{ksize}-sc{scaled}.comparison-info.csv"),
        ani_csv=os.path.join(out_dir, "sourmash", "{basename}.k{ksize}-sc{scaled}.cANI.csv")
    conda: "conf/env/sourmash.yml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        time=240,
        partition="high2", #"low2",
    log: os.path.join(logs_dir, "sourmash", "{basename}.k{ksize}-sc{scaled}.api-compare.log")
    benchmark: os.path.join(logs_dir, "sourmash", "{basename}.k{ksize}-sc{scaled}.api-compare.benchmark")
    shell:
        """
        python sourmash-api-compare.py {input.sigs} -c {input.comparisons} \
                                       -o {output.compare_csv}  --ani-threshold 0 \
                                       -s {wildcards.scaled} --ani-csv {output.ani_csv} 2> {log}
        """

rule pyani_index_and_createdb:
    input:
        classes=os.path.join(out_dir, "genomes", "py.classes.txt"),
        labels=os.path.join(out_dir, "genomes", "py.labels.txt")
    output:
        classes=os.path.join(out_dir, "genomes", "classes.txt"),
        labels=os.path.join(out_dir, "genomes", "labels.txt"),
        db=os.path.join(out_dir, "genomes", ".pyanidb"),
    params:
        genome_dir = os.path.join(out_dir, 'genomes'),
        pyanidb = os.path.join(out_dir, 'genomes',f".pyanidb"),
        classes_basename = "classes.txt",
        labels_basename = "labels.txt"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *6000,
        time=1200,
        partition="bmh",
    log: os.path.join(logs_dir, "pyani", "index-and-createdb.log")
    benchmark: os.path.join(logs_dir, "pyani", "index-and-createdb.benchmark")
    conda: "conf/env/pyani0.3.yml"
    shell:
        """
        pyani index -i {params.genome_dir} --classes {params.classes_basename} --labels {params.labels_basename}
        pyani createdb --dbpath {params.pyanidb} -v -l {log}
        """

rule pyANI_ANIb:
    input:
        classes=os.path.join(out_dir, "genomes", "py.classes.txt"),
        labels=os.path.join(out_dir, "genomes", "py.labels.txt"),
        idx_classes=os.path.join(out_dir, "genomes", "classes.txt"),
        idx_labels=os.path.join(out_dir, "genomes", "labels.txt"),
    output:
        covF= os.path.join(out_dir, "pyani/ANIb_results","ANIb_alignment_coverage.tab"),
        lenF= os.path.join(out_dir, "pyani/ANIb_results","ANIb_alignment_lengths.tab"),
        hadF= os.path.join(out_dir, "pyani/ANIb_results","ANIb_hadamard.tab"),
        idF=  os.path.join(out_dir, "pyani/ANIb_results","ANIb_percentage_identity.tab"),
        seF=  os.path.join(out_dir, "pyani/ANIb_results","ANIb_similarity_errors.tab"),
        bn =  os.path.join(out_dir, "pyani/ANIb_results","blastn_output.tar.gz"),
    threads: 30
    resources:
        mem_mb=lambda wildcards, attempt: attempt *50000,
        time=120000,
        partition="bmh",#"med2"
    params:
        genome_dir = lambda w: os.path.join(out_dir, 'genomes'),
        output_dir = lambda w: os.path.join(out_dir, 'pyani', "ANIb_results"),
    log: os.path.join(logs_dir, "pyani_anib", "pyANI-anib.log")
    benchmark: os.path.join(logs_dir, "pyani_anib", "pyANI-anib.benchmark")
    conda: "conf/env/pyani0.2.yml"
    shell:
        """
        average_nucleotide_identity.py -i {params.genome_dir} \
             -o {params.output_dir} -m ANIb -v --workers {threads} \
             --labels {input.labels} --classes {input.classes} \
             --force > {log}
        """

localrules: aggregate_anib
rule aggregate_anib:
    input:
        comparisons=comparison_file,
        covF= os.path.join(out_dir, "pyani/ANIb_results", "ANIb_alignment_coverage.tab"),
        lenF= os.path.join(out_dir, "pyani/ANIb_results","ANIb_alignment_lengths.tab"),
        hadF= os.path.join(out_dir, "pyani/ANIb_results","ANIb_hadamard.tab"),
        idF=  os.path.join(out_dir, "pyani/ANIb_results","ANIb_percentage_identity.tab"),
        seF=  os.path.join(out_dir, "pyani/ANIb_results","ANIb_similarity_errors.tab"),
    output:
        os.path.join(out_dir, "pyani", "pyani.ANIb.csv"),
    params:
        results_dir =  os.path.join(out_dir, "pyani/ANIb_results"),
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        time=1200,
    shell:
        """
        python aggregate-pyani-results.py {params.results_dir} --pyani-version 0.2 -c {input.comparisons} -o {output}
        """

rule fastANI:
    input: os.path.join(out_dir, "genomes", "genome-filepaths.txt")
    output: os.path.join(out_dir, "fastani","fastani.tsv"),
    threads: 12
    resources:
        mem_mb=lambda wildcards, attempt: attempt *20000,
        time=12000,
        partition="high2",
    log: os.path.join(logs_dir, "fastani", "fastani.log")
    benchmark: os.path.join(logs_dir, "fastani", "fastani.benchmark")
    conda: "conf/env/fastani.yml"
    shell:
        """
        fastANI --threads {threads} --ql {input} --rl {input} -o {output} 2> {log}
        """


localrules: parse_fastani
rule parse_fastani:
    input: 
        fastani=os.path.join(out_dir, "fastani", "fastani.tsv"),
        comparisons=comparison_file,
    output: os.path.join(out_dir, "fastani", "fastani.ANI.csv")
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        time=12000,
    log: os.path.join(logs_dir, "fastani", "parse_fastani.log")
    benchmark: os.path.join(logs_dir, "fastani", "parse_fastani.benchmark")
    shell:
        """
        python parse-fastani-results.py --fastani {input.fastani} -c {input.comparisons} \
                                        --output-csv {output} 2> {log}
        """

localrules: download_usearch
rule download_usearch:
    output: "conf/env/usearch11.0.667_i86linux32",
    params:
        link="http://www.drive5.com/downloads/usearch11.0.667_i86linux32.gz",
        dl_name="usearch11.0.667_i86linux32.gz",
    shell:
        """
        curl -JLO {params.link}
        gunzip -c {params.dl_name} > {output}
        chmod a+rx {output}
        rm {params.dl_name}
        """

rule mash_sketch:
    input: os.path.join(out_dir, "genomes", "{genome}_genomic.fna"),
    output:
        os.path.join(out_dir, "mash", "{genome}_genomic.msh"),
    conda: "conf/env/mash.yml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        time=240,
        partition="high2",#"low2",
    log: os.path.join(logs_dir, "mash", "{genome}.sketch.log")
    benchmark: os.path.join(logs_dir, "sourmash", "{genome}.sketch.benchmark")
    shell:
        """
        mash sketch {input} -o {output} 2> {log}
        """

rule mash_dist:
    input:
        gA= os.path.join(out_dir, "mash", "{gA}_genomic.msh"),
        gB= os.path.join(out_dir, "mash", "{gB}_genomic.msh"),
    output: os.path.join(out_dir, "mash", "{gA}_x_{gB}.mashdist.tsv") 
    conda: "conf/env/mash.yml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        time=240,
        partition="high2",#"low2",
    log: os.path.join(logs_dir, "mash", "{gA}_x_{gB}.dist.log")
    benchmark: os.path.join(logs_dir, "mash", "{gA}_x_{gB}.dist.benchmark")
    shell:
        """
        mash dist {input.gA} {input.gB} > {output} 2> {log}
        """

localrules: aggregate_pairwise_mash
rule aggregate_pairwise_mash:
    input: expand(os.path.join(out_dir, "mash", "{comp}.mashdist.tsv"), comp=COMPARISONS),
    output: os.path.join(out_dir, "mash", "mash.ANI.csv")
    run:
        # Mash columns: Reference-ID, Query-ID, Mash-distance, P-value, and Matching-hashes
        results = []
        header = ["comparison_name", "identA", "identB", "MashDistance", "Mash_pvalue", "Mash_matchinghashes", "MashANI"]
        for inF in input:
            inF = str(inF)
            comparison_name = os.path.basename(inF).rsplit(".mashdist.tsv")[0]
            identA, identB = comparison_name.split("_x_")
            with open(str(inF), 'r') as res:
                result = [comparison_name, identA, identB] + res.readlines()[-1].rsplit('\n')[0].split('\t')[2:]
                mashdist = float(result[3])
                est_ani = 1 - mashdist
                result.append(est_ani)
                results.append(result)
        with open(str(output), 'w') as csv_out:
            csvwriter = csv.writer(csv_out, delimiter=',')
            csvwriter.writerow(header)
            for res in results:
                csvwriter.writerow(res)


# need to dl OAU manually and scp to hpc. download link: http://www.ezbiocloud.net/download2/download_oau
rule orthoani_pairwise:
    # NOTE: do not run more than one at a time -- tmp files may collide (I assume) -->  get "-1" as an ANI value.
    input: 
        gA= ancient(os.path.join(out_dir, "genomes", "{gA}_genomic.fna")),
        gB= ancient(os.path.join(out_dir, "genomes", "{gB}_genomic.fna")),
        usearch="conf/env/usearch11.0.667_i86linux32",
    output: os.path.join(out_dir, "orthoani", "{gA}_x_{gB}.orthoani.txt")
    params:
        oau=os.path.abspath("conf/env/OAU.jar"),
        usearch=os.path.abspath("conf/env/usearch11.0.667_i86linux32"),
    threads: 1
    conda: "conf/env/orthoani.yml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        time=240,
        partition="high2", #"low2",#"med2",
    log: os.path.join(logs_dir, "orthoani", "{gA}_x_{gB}.orthoani.log")
    benchmark: os.path.join(logs_dir, "orthoani", "{gA}_x_{gB}.orthoani.benchmark")
    shell:
        """
        java -jar {params.oau:q} --upath {params.usearch} -f1 {input.gA} -f2 {input.gB} \
                   --format list --threads {threads} > {output} 2> {log}
        """
#--tmp {params.tmp} 
# -fmt,--format <arg>    Output format (optional, default: list, possible options: list, matrix, json) 

localrules: aggregate_pairwise_orthoani
rule aggregate_pairwise_orthoani:
    input: expand(os.path.join(out_dir, "orthoani", "{comp}.orthoani.txt"), comp=COMPARISONS),
    output: os.path.join(out_dir, "orthoani", "orthoani.ANI.csv")
    run:
        results = []
        header = ["comparison_name", "identA", "identB", "orthoANI_value", "orthoANI_avg_aligned_length", "orthoANI_query_coverage", "orthoANI_subject_coverage", "orthoANI_query_length", "orthoANI_subject_length"]
        for inF in input:
            inF = str(inF)
            comparison_name = os.path.basename(inF).rsplit(".orthoani.txt")[0]
            identA, identB = comparison_name.split("_x_")
            with open(str(inF), 'r') as res:
                result = [comparison_name, identA, identB] + res.readlines()[-1].rsplit('\n')[0].split('\t')[1:]
                results.append(result)
        with open(str(output), 'w') as csv_out:
            csvwriter = csv.writer(csv_out, delimiter=',')
            csvwriter.writerow(header)
            for res in results:
                csvwriter.writerow(res)


localrules: combine_ani
rule combine_ani:
    input: 
        fastani= os.path.join(out_dir, "fastani", "fastani.ANI.csv"),
        pyani=os.path.join(out_dir, "pyani", "pyani.ANIb.csv"),
        orthoani=os.path.join(out_dir, "orthoani", "orthoani.ANI.csv"),
        mash=os.path.join(out_dir, "mash", "mash.ANI.csv"),
        sourmash=expand(os.path.join(out_dir, "sourmash", f"{basename}.k{{ksize}}-sc{{scaled}}.cANI.csv"), ksize=KSIZE, scaled=SCALED),
    output: os.path.join(out_dir, "combinedANI.csv")
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        time=12000,
    log: os.path.join(logs_dir, "combine-ani", "combine-ani.log")
    benchmark: os.path.join(logs_dir, "combine-ani", "combine-ani.benchmark")
    shell:
        """
        python combine-ani.py --fastani {input.fastani} --pyani {input.pyani} \
                              --sourmash {input.sourmash} --orthoani {input.orthoani} \
                              --mash {input.mash} --output-csv {output} 2> {log}
        """


## unused (but working) rules ##

#rule orthoani:
#    input: 
#        genomes=ancient(expand(os.path.join(out_dir, "genomes", "{acc}_genomic.fna"), acc=IDENTS)),
#        usearch="conf/env/usearch11.0.667_i86linux32",
#    output: os.path.join(out_dir, "orthoani", "orthoani.ANI.tsv")
#    params:
#        genome_dir=os.path.abspath(os.path.join(out_dir, "genomes")),
#        outdir=os.path.abspath(os.path.join(out_dir, "orthoani")),
#        oau=os.path.abspath("conf/env/OAU.jar"),
#        usearch=os.path.abspath("conf/env/usearch11.0.667_i86linux32"),
#        #tmp="/scratch/ntpierce"
#        #tmp="/tmp"
#    threads: 6
#    conda: "conf/env/orthoani.yml"
#    resources:
#        mem_mb=lambda wildcards, attempt: attempt *3000,
#        time=12000,
#        partition="med2",
#    log: os.path.join(logs_dir, "orthoani", "orthoani.log")
#    benchmark: os.path.join(logs_dir, "orthoani", "orthoani.benchmark")
#    shell:
#        """
#        java -jar {params.oau:q} --upath {params.usearch} --fastadir {params.genome_dir} \
#                   --format list --threads {threads} -o {params.outdir} 2> {log}
#        """
#--tmp {params.tmp} 
# -fmt,--format <arg>    Output format (optional, default: list, possible options: list, matrix, json) 


