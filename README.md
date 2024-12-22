# Random UniProtKB sampler

Decimate the UniRef50 database by a factor of 1/100000 and randomly sample EC (Enzyme commission) IDs, 1000 times

`python mc_sample.py .5 1000 10000`

Prequisite:

- `pip install tqdm SPARQLWrapper`
