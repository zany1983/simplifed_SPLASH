# simplifed_SPLASH

### requirments 

To completely run analysis from fastq to results in the paper, the following tools are required:

Microsoft Excel 2017, Python 3.6,R 3.5.2,  numpy 1.15.4, scipy 1.3.0,pandas 0.23.4, dplyr 0.8.0.1, Hyb, DeSeq2, ggplot2 3.1.0

### installation 
All the custome scripts could be downloaded and run directly in R or python without installation

### main pipeline
This repository contains whole scripts to dected specific interaction (viewpoint analysis), differential interractions, count suporting chimeras for secondary structure and also domain calling and 3D modeling from chimeras detected by Hyb pipeline.

the raw and processed data of SARS-CoV-2 simplified SPLASH data can be downloaded from GEO (accession number GEO: GSE164565). https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164565

#### calling chimeras
first chimeras could be detected as follow:

```{bash}
# defalt (stringent) pipeline
hyb analyse in=sample_R2.fq db=SARS-CoV-2_no_polyA format=fastq align=bowtie2 eval=0.001

# relax (gmax=20) pipeline
hyb analyse in= sample_R2.fq db=SARS-CoV-2_no_polyA format=fastq align=bowtie2 eval=0.001 gmax=20

```
For Hyb pipe, please refer to:(https://github.com/gkudla/hyb).

The processed files could also be found in https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164565 , as processed files

#### heatmaps for contact matrix
After calling chimeras, you can use merge_cdt.Rmd scripts to export contact matrix for heatmap, either globaly or for specific interactions.

For convienence, the processed cdt files containing bined contact matrix were uploaded to "cdtfiles" directory.
The annotation based on NCBI and Narry Kim were also included in annotations directory

#### detecting differential inteeractions
For detecting differential inteeractions between groups, DESeq2.Rmd could be used.

counting_chimeras_for_sencondary_structure.Rmd is designated to count supportive chimeras for each interaction.

For viewpoint analysis on specific loci, TRSl_VP.Rmd could be used. 

#### domain boudary analyses
To call domain boudaries, please enter DomainCalling subdirectory, cworld is implented: https://github.com/dekkerlab/cworld-dekker

To call domain boudaries, please enter 3DStructure subdirectory.
