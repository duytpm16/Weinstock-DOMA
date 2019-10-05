### This script runs additive permutations in parallel using R script from:
#     https://github.com/duytpm16/qtl2-HPC-pipeline/blob/master/R%20scripts/qtl2_scan1perm.R


#PBS -q batch
#PBS -l nodes=1:ppn=1

module load R/3.5.1



job_name='doma_additive_perm'


for i in {1..95}
do
  echo "#PBS -l nodes=1:ppn=8
#PBS -q batch
#PBS -l walltime=72:00:00
module load R/3.5.1
viewer_data='weinstock_doma_viewer_v2.Rdata'
dataset_expr='dataset.doma.otu|data|norm'
num_cores='8'
int_name='NA'
perm_run='1000'
chunk_number=${i}
chunk_size='4'


Rscript qtl2_scan1perm.R viewer_data=\$viewer_data dataset_expr=\$dataset_expr num_cores=\$num_cores int_name=\$int_name perm_run=\$perm_run chunk_number=\$chunk_number chunk_size=\$chunk_size" >> "${job_name}_${i}.sh"
qsub "${job_name}_${i}.sh"
