## Requirements

Python requirements: `screed`, `mmh3`, `scipy`, `numpy`

System requirements: `mash` needs to be installed

### Result generation

```
cd result-generation
python simulated_mutation_experiment.py
```

The results will be in the file named `stats`

### Plotting

```
cd plotting
cp ../result-generation/stats stats
python plotter-mutated-staphylo-expt.py
```
