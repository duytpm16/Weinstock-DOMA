### This script runs qtl mapping for both additive and sex interaction using R script located:
#     https://github.com/duytpm16/qtl2-HPC-pipeline/blob/master/R%20scripts/qtl2_scan1.R


#PBS -q batch
#PBS -l nodes=1:ppn=1

module load R/3.5.1


taxa=('otu' 'genus' 'family' 'order' 'class' 'phylum')
cov=('NA' 'sex', 'cohort.age')



for i in "${taxa[@]}"
do
  for k in "${cov[@]}"
  do
    job_name="doma_${i}_${k}"

    echo "#PBS -q batch
    #PBS -l nodes=1:ppn=2
    #PBS -l walltime=24:00:00
    module load R/3.5.1

    viewer_data='weinstock_doma_viewer_v2.Rdata'
    dataset_expr='dataset.doma.${i}|data|norm'
    num_cores='2'
    int_name=${k}
    chunk_number='NA'
    chunk_size='NA'
    
    Rscript qtl2_scan1.R viewer_data=\$viewer_data dataset_expr=\$dataset_expr num_cores=\$num_cores int_name=\$int_name chunk_number=\$chunk_number chunk_size=\$chunk_size" >> "${job_name}.sh"
    qsub "${job_name}.sh"
  done
done
