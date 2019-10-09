#PBS -l nodes=1:ppn=1
#PBS -q batch
#PBS -l walltime=72:00:00
module load R/3.5.1



taxa=('phylum' 'class' 'order' 'family' 'genus' 'otu')
month=(6 12 18)
cov=('NA' 'sex')



for i in "${taxa[@]}"
do
  for j in "${month[@]}"
  do
    for k in "${cov[@]}"
    do
      job_name="${i}_m${j}_${k}_findpeaks.sh"
      
      
      if [ ${k} == 'NA' ]
      then
         peaks_name='additive'
         int_scan='NA'
      elif [ ${k} == "cohort.age" ]
      then
         peaks_name='age_int'
         int_scan="${i}_m${j}_cohort.age_int_scan.rds"
      else
         peaks_name="sex_int"
         int_scan="${i}_m${j}_sex_int_scan.rds"
      fi




      echo "#PBS -l nodes=1:ppn=1
      #PBS -q batch
      #PBS -l walltime=72:00:00 
      module load R/3.5.1

      viewer_data=weinstock_doma_by_months_viewer_v2.Rdata
      scan1_mat="${i}_m${j}_additive_scan.rds"
      dataset_expr='dataset.${i}.m${j}|data|norm'
      thr='6'
      num_cores='1'
      type_scan="${peaks_name}"
      prob='.95'
      cis_threshold='NA'
      int_mat="${int_scan}"



      Rscript qtl2_findpeaks.R viewer_data=\$viewer_data scan1_mat=\$scan1_mat dataset_expr=\$dataset_expr thr=\$thr num_cores=\$num_cores type_scan=\$type_scan prob=\$prob cis_threshold=\$cis_threshold int_mat=\$int_mat" >> $job_name
      sleep 5m
      qsub $job_name
      
      
      
    done
  done
done
