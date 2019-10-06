### This script runs sex interaction scan using R script located:
#     https://github.com/duytpm16/qtl2-HPC-pipeline/blob/master/R%20scripts/qtl2_scan1.R



#PBS -q batch
#PBS -l nodes=1:ppn=2
#PBS -l walltime=24:00:00

module load R/3.5.1

viewer_data='weinstock_doma_viewer_v1.Rdata'
dataset_expr='dataset.doma.phylum|data|rz'
num_cores='2'
int_name='sex'
chunk_number='NA'
chunk_size='NA'



Rscript qtl2_scan1.R viewer_data=$viewer_data dataset_expr=$dataset_expr num_cores=$num_cores int_name=$int_name chunk_number=$chunk_number chunk_size=$chunk_size
