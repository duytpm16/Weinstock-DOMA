### This shell script is used to with the R script located at
#      https://github.com/duytpm16/qtl2-HPC-pipeline/tree/master/R%20scripts/qtl2_scan1.R



#PBS -q batch
#PBS -l nodes=1:ppn=8
#PBS -l walltime=24:00:00

module load R/3.5.1


viewer_data='weinstock_doma_viewer_v2.Rdata'
dataset_expr='dataset.doma.otu|data|norm'
num_cores='8'
int_name='NA'
chunk_number='NA'
chunk_size='NA'



Rscript qtl2_scan1.R viewer_data=$viewer_data dataset_expr=$dataset_expr num_cores=$num_cores int_name=$int_name chunk_number=$chunk_number chunk_size=$chunk_size
