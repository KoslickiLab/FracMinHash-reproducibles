## Requirements

Python requirements: `screed`, `mmh3`, `scipy`, `numpy`, `java`

System requirements: `mash` needs to be installed

### Select pairs of genomes randomly

```
cd result-generation
cd data
gunzip *
cd ..
python generate_representative_pairs.py
```

After running for a large number of iterations, the generated pairs were taken and hardcoded into the code that actually compares pairwise mutation distance

### Pairwise distance calculation between selected pairs
```
cd result-generation
python analyze_representative_pairs.py
```

Output will be stored in file `results`

### Plotting

```
cd plotting
cp ../result-generation/results results
manually delete the file names in the results file
python plotter-pairwise-genomes-expt.py
```
