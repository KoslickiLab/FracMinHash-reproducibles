# Select a subset of GTDB genomes and compare with a series of ANI tools

Final genome selection process:

Using GTDB-rs207, I selected species with the largest number of genomes (7 bacteria, 3 archaea). I chose the representative species genome from each of those species as the "anchor" to build an evolutionary path. 

Those 10 genomes are: 
- 'GCF_000742135.1 Klebsiella pneumoniae strain=ATCC 13883, ASM74213v1'
- 'GCF_000195955.2 Mycobacterium tuberculosis H37Rv strain=H37Rv, ASM19595v2'
- 'GCF_001457615.1 Pseudomonas aeruginosa strain=NCTC10332, NCTC10332'
- 'GCF_001457635.1 Streptococcus pneumoniae strain=NCTC7465, NCTC7465'
- 'GCF_001027105.1 Staphylococcus aureus subsp. aureus DSM 20231 strain=DSM 20231, ASM102710v1'
- 'GCA_003162175.1 Euryarchaeota archaeon, 20100900_E2D'
- 'GCF_003697165.2 Escherichia coli DSM 30083 = JCM 1649 = ATCC 11775 strain=ATCC 11775, ASM369716v2'
- 'GCF_000970205.1 Methanosarcina mazei S-6 strain=S-6, ASM97020v1'
- 'GCF_000006945.2 Salmonella enterica subsp. enterica serovar Typhimurium str. LT2, ASM694v2'
- 'GCF_000012285.1 Sulfolobus acidocaldarius DSM 639 strain=DSM 639, ASM1228v1'

Then I selected three non-representative genomes at each taxonomic rank in the path, e.g. three genomes in the same genus but different species, three in same family but different genera, etc. 
This selection was done here: https://github.com/bluegenes/2022-focused-cani-comparisons/blob/main/05.select-comparisons-by-common-species.ipynb
Initial plots generated here: https://github.com/bluegenes/2022-focused-cani-comparisons/blob/main/06.plot-comparisons.ipynb

Here's the snakefile used to run comparisons using all programs: https://github.com/bluegenes/2022-focused-cani-comparisons/blob/main/ani-compare.snakefile

I used the sourmash python api to run sourmash comparisons (so I could store # hashes, # intersected hashes, etc); script here: https://github.com/bluegenes/2022-focused-cani-comparisons/blob/main/sourmash-api-compare.py

You can also find all intermediate files here: https://osf.io/d3c85/ in the output.ani-commonsp10-evolpath folder
