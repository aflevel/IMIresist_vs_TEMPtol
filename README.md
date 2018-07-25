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

The scripts can then be executed using the Rscript command from the directory that was unpacked

The GNU parallel package was extensively used to parallelise computation [O. Tange (2011): GNU Parallel - The Command-Line Power Tool,
;login: The USENIX Magazine, February 2011:42-47.]

## Content:

- **GISanalysis.r**
  R script used to design the fly sampling using Bioclim data (Hijmans et al 2006) and assess the level of pesticide pollution using water analysis data (Vörösmarty et al. 2010)
  
- ****
