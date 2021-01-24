#step1 Construction of single chromosome N*N matrix
python3 iced_asis2cis_npy_matrix_filterZeros.py -m ../Data/C_100_Turn.matrix -b ../Data/chrSC2_100_abs.bed -o ../Analysis -s C -c chrSC2
#step2 3D modeling
export LD_LIBRARY_PATH='/home/xiedejian/miniconda3/lib':$LD_LIBRARY_PATH
/home/xiedejian/miniconda3/envs/Python2/bin/python /home/xiedejian/miniconda3/bin/pastis-pm1 ../Analysis/
