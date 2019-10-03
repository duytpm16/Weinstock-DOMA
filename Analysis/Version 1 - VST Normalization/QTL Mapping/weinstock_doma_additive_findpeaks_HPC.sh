### This script finds additive peaks above 6 using R script located:
#     https://github.com/duytpm16/qtl2-HPC-pipeline/blob/master/R%20scripts/qtl2_findpeaks.R


#PBS -l nodes=1:ppn=8
#PBS -q batch
#PBS -l walltime=72:00:00
module load R/3.5.1


viewer_data='weinstock_doma_viewer_v1.Rdata'
scan1_mat='doma_otu_additive_scan.rds'
dataset_expr='dataset.doma.otu|data|rz'
thr='6'
num_cores='8'
type_scan='additive'
prob='.95'
cis_threshold='NA'
int_mat='NA'



Rscript qtl2_findpeaks.R "viewer_data=$viewer_data" "scan1_mat=$scan1_mat" "dataset_expr=$dataset_expr" "thr=$thr" "num_cores=$num_cores" "type_scan=$type_scan" "drop=$drop" "cis_threshold=$cis_threshold" "int_mat=$int_mat"
