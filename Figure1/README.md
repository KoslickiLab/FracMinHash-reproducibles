## Requirements

Python requirements: `screed`, `mmh3`

System requirements: `mash` needs to be installed


### Result generation

```
cd result-generation
python experiment_containment.py > containment_results.txt
```

### Plotting

```
cd plotting
cp ../result-generation/containment_results.txt containment_results.txt
python plotter-containment-expt.py
```