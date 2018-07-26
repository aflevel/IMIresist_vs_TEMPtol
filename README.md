<h1 align="center">IMIresist_vs_TEMPtol</h1>
This suite of R scripts allows to replicate the set of analysis proposed in Fournier-Level et al. (2018, submitted)

## Usage and Dependencies:
The scripts require to first source and unzip the data archive that can be found here: 
The following R packages need to be install to exectute the scripts
- lme4
- vegan
- fields
- raster
- plotrix
- FactoMineR

The scripts can then be executed using the Rscript command from the directory that was unpacked

The GNU parallel package was extensively used to parallelise computation [O. Tange (2011): GNU Parallel - The Command-Line Power Tool,login: The USENIX Magazine, February 2011:42-47.]

## Content:

- **GISanalysis.r**:
  R script used to design the fly sampling using Bioclim data (Hijmans et al 2006) and assess the level of pesticide pollution using water analysis data (Vörösmarty et al. 2010).
  It returns a plot the sampling zones, a plot of the pesticide load Worldwide (normalised unit) and a table of the pesticide load data for the population sampled
  
- **Phenotype_Analysis.r**:
  Rscript used to perform the linear modelling of the phenotype variation in longevity in response to imidacloprid and temperature.
  It returns a serie of tests for the effect of the continent, temperature and precipitation regime at the location of origin on longevity and poymorphism. It also plots the cumulative distribution of time-to-death for each assay (stored in the Pheno_data directory), the distribution of time-to-death for the entire experiment and the population and line effect on time-to-death.
  Alternatively, lines 46 64-66 and 72 can be commented out to perform a maximum likelihood estimation of mean, std dev. minimun and maximum for the insecticide data only (no controls) to return the pheno.py design files required tfor he GWAlpha association analysis.
