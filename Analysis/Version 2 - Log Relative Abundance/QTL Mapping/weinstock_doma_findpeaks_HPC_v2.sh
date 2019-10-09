#PBS -l nodes=1:ppn=1
#PBS -q batch
#PBS -l walltime=72:00:00
module load R/3.5.1

taxa=('phylum' 'class' 'order' 'family' 'genus' 'otu')
cov=('NA' 'sex' 'cohort.age')

for i in "${taxa[@]}"
do
  for k in "${cov[@]}"
  do
    job_name="${i}_${k}_findpeaks.sh"

     if [ ${k} == 'NA' ]
     then
         peaks_name='additive'
         int_scan='NA'
     elif [ ${k} == "cohort.age" ]
     then
         peaks_name='age_int'
         int_scan="doma_${i}_${k}_int_scan_v2.rds"
     else
         peaks_name="sex_int"
         int_scan="doma_${i}_${k}_int_scan_v2.rds"
     fi



      echo "#PBS -l nodes=1:ppn=1
      #PBS -q batch
      #PBS -l walltime=72:00:00 
      module load R/3.5.1

      viewer_data=weinstock_doma_viewer_v2.Rdata
      scan1_mat="doma_${i}_additive_scan_v2.rds"
      dataset_expr='dataset.doma.${i}|data|norm'
      thr='6'
      num_cores='1'
      type_scan="${peaks_name}"
      prob='.95'
      cis_threshold='NA'
      int_mat="${int_scan}"

      Rscript qtl2_findpeaks.R viewer_data=\$viewer_data scan1_mat=\$scan1_mat dataset_expr=\$dataset_expr thr=\$thr num_cores=\$num_cores type_scan=\$type_scan prob=\$prob cis_threshold=\$cis_threshold int_mat=\$int_mat" >> $job_name
      sleep 7m
      qsub $job_name
  done
done
~       
