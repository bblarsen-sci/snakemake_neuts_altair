# Snakemake_neut
Plot neutralization curves with Altair

This snakemake pipeline uses [neut_curve](https://jbloomlab.github.io/neutcurve/) to fit neutralization curves to data and [Altair](https://altair-viz.github.io/) to plot the results.

## Installation and Usage
clone the repository
```bash
git clone https://github.com/bblarsen-sci/luciferase_neut_snakemake.git
```

Then run the pipeline with snakemake
```bash
./run_pipeline.bash
```

The pipeline works by reading any .csv files in /data and generates plots and a .csv of the fit paramaters in /results.

**Important**
The pipeline runs by parsing the name of the file placed in /data. Valid strings are:
- 'serum' will assume you are working with serum and adjust labels
- 'antibody' will assume you are working with antibodies and adjust to be in µg/mL scale
- 'receptor' will assume you are working with receptors and adjust to be in µM scale
- 'facet' if facet is in the name, it will facet the plots 

## Data format
The data should be in a csv file with the following columns. Note 'serum', is used generally here, it could be an antibody or receptor, or anything that neutralizes virus.
- serum: the serum name
- virus: the virus name
- replicate: the replicate number
- concentration: the concentration of serum/antibody/etc.
- fraction infectivity: the fraction of infectivity at that concentration

## Examples
Examples of data and output are in /data and /results respectively.