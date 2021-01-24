#step1 Construction of single chromosome N*N matrix
python3 Condense_CisMatrix.py -s C -b ../Data/Raw_index.bed -m ../Data/C_raw.matrix -g SARS_COV2 -c chrSC2 -o ../Analysis
#step2 Insulation calling
perl matrix2insulation.pl -i ../Analysis/C.SARS_COV2_chrSC2.SARS_COV2_chrSC2.matrix --is 500 --ids 150 --nt 1  -o../Analysis
#step3 Filter Robust Domain boundaries
python3 correcting_boundaries_by_insulation.py -i ../Analysis/C_chrSC2--is501--nt1.0--ids161--ss1--immean.insulation -nt 1.0  -o ../Analysis/
#step4 Correction of boundary parameters
python3 tranform_bound.py -b C_chrSC2--is501--nt1.0--ids161--ss1--immean.insulation.boundaries -c chrSC2 -s C -r 10 -o ../Analysis -cs ../Data/chrom_SC2.sizes

