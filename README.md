<h1 align="center">IMIresist_vs_TEMPtol</h1>
This suite of R scripts allows to replicate the analysis proposed in Fournier-Level et al. (2018, submitted)
If you require some futher support or information please email Alex Fournier-Level (afournier@unimelb.edu.au)

## Installation Usage and Dependencies:
The scripts require to first source and unzip the data archive that can be found here: **<not public yet, please email afournier@unimelb.edu.au>**

The following R packages need to be installed to execute the scripts:
- lme4
- vegan
- fields
- raster
- plotrix
- FactoMineR

The scripts can then be executed using an Rscript command from the directory that was unpacked or through the R GUI having set the working directory to were the archive was unpacked.

The GNU parallel package was extensively used to parallelise computation [O. Tange (2011): GNU Parallel - The Command-Line Power Tool,login: The USENIX Magazine, February 2011:42-47.]

## Content:

- **GISanalysis.r**:
  R script used to determine the fly sampling locations using Bioclim data (Hijmans et al 2006) and assess the level of pesticide pollution using water analysis data (Vörösmarty et al. 2010).
  It returns a plot of the potential the sampling zones, a plot of the pesticide load Worldwide (normalised unit) and a table of the pesticide load data for the population sampled.
  
- **Phenotype_Analysis.r**:
  Rscript used to perform the linear modelling of the phenotype variation for time-to-death in response to imidacloprid and temperature.
  It returns a serie of tests for the effect of the continent, temperature and precipitation regime at the location of origin on time-to-death and genomic poymorphism. It also plots the cumulative distribution of time-to-death for each assay (stored in the Pheno_data directory), the distribution of time-to-death for the entire experiment and the population and line effect on time-to-death.
  Alternatively, lines 46 64-66 and 72 can be commented out to perform a maximum likelihood estimation of mean, std dev. minimun and maximum time-to-death for the insecticide data only (no controls) and return the pheno.py design files required for the GWAlpha genome-wide association analysis (https://github.com/aflevel/GWAlpha).

- **CoverageEstimation.r**:
  Rscript used to estimate the distribution of coverage for every sequencing library. The coverage estimation was made by windows of 500bp using the pysamstats package (https://github.com/alimanfoo/pysamstats). 
  It returns a LibMedCov file; for the GWAS analysis using the GWAlpha package (https://github.com/aflevel/GWAlpha), this data was used to only retain the genomic regions that showed a coverage within the interquartile distribution of each pool (q25).

- **GWAS_Analysis.r**:
  Rscript used to visualise the GWAS analysis for time-to-death exposed to 1000ppm imidacloprid at either 20°C or 30°C performed on SNP+InDels, Chromosomal Inversion, Transposable Elements and Copy Number Variation.
  It prints-out the list of candidate genes for each GWAS and generates separate plots for the SNP/Inv/TE and for the CNV associations. It also generates an RData archive containing the list of candidate genes for all associations and a GWAlist_dm5.57.txt files for the candidate genes for the SNP/Inv/TE association only. If Graph_dia in line 17 is True, diagnostic plots for the distribution of the test statistics and effect of coverage on associations are generated.
  It also prints out a table of the association between the abundance of microbes or endosymbionts of Drosophila and time-to-death as well as a table reprorting the abundance of each microbe in each population (as a ratio of the number of microbe sequences relative to the number of drosophila sequences).
  Finally it generates a plot of the density of candidate genes along the genome.
 
 - **PopGen_Analysis.r**:
  Rscript used to analyse the pattern of population genetic diversity within and among populations.
  It prints a chi-squared test of whether the set of alleles at candidate genes associated with resistance in each GWAS are mostly ancestral or derived states.
  It returns plots showing the level of genetic differentiation along the genome for the set of 16 populations (Fst estimated using the Popoolation package: https://sourceforge.net/p/popoolation2/wiki/Main/) together with the contribution of georgraphic and climatic factor in this differentiation. It also returns the position and strength of putative selective sweeps (estimated using a Hidden Markov Model impelemented in Pool-HMM: https://forge-dga.jouy.inra.fr/projects/pool-hmm). Finally it plots the results of the MonteCralo simulations testing that the pattern of differentiation and selction at insecticide resistance candidate loci is not compatible with randomness or neutrality.
  
- **pHMM_MonteCarloAnalysis.r and pHMM_Ancestral_MonteCarloAnalysis.r**:
  Rscripts used to perform the MonteCarlo simulation testing the enrichment for candidate genes among the most intensely or least intensely (when resistance is ancestral) selected loci.
  It uses as argument the number of GWAS for which a gene should be shown as associated to be included in the candiadate gene list (eg: 3 if a given gene is to be identified in 3 GWAS) and returns an RData archive with the observed and the permuted (1000 distinct permutations by default that can be modified in line 53).

 - **Functional_Analysis.r**:
  Rscript used to analyse the results from the mutant/deletion vs wildtype comparisons. It generates barplots for the parwise comparison between pairs of lines and a box-and-whisker plot for the overall effect of a specific gene disruption.
  
