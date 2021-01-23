# simplifed_SPLASH

These repository contains whole scripts to dected specific interaction (viewpoint analysis), differential interractions, count suporting chimeras for secondary structure and also domain calling and 3D modeling from chimeras detected by Hyb pipeline.

the raw and processed data of SARS-CoV-2 simplified SPLASH data can be downloaded from GEO (accession number GEO: GSE164565).

first chimeras could be detected as follow:

```{bash}
# defalt (stringent) pipeline
hyb analyse in=sample_R2.fq db=SARS-CoV-2_no_polyA format=fastq align=bowtie2 eval=0.001

# relax (gmax=20) pipeline
hyb analyse in= sample_R2.fq db=SARS-CoV-2_no_polyA format=fastq align=bowtie2 eval=0.001 gmax=20

```


After calling chimeras, you can use merge_cdt.Rmd scripts to export contact matrix for heatmap, either globaly or for specific interactions.

For detecting differential inteeractions between groups, DESeq2.Rmd could be used.

counting_chimeras_for_sencondary_structure.Rmd is designated to count supportive chimeras for each interaction.

For viewpoint analysis on specific loci, TRSl_VP.Rmd could be used. 

