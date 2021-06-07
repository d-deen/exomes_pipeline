#! /bin/bash
#SBATCH -A rtmngs
#SBATCH -p bigmem
#SBATCH --mem=10G


In=$1
Out=$2
Job=$3


script_path="/mnt/nfs/home/ndd73/exome_pipeline"


JobID1="Set_up_folder and copy/concatenate files"
JobID2="Running_fastqQC"
JobID3="BWA alignment and duplicate removal"
JobID4="gVCF generation"
JobID5="Combine gVCF, annotate and make a report"
JobID6='Annotate'
JobID7='Run a QC report'
JobID8='Mito SNPs'
JobID9='Annotate mito'


#job1=$(srun --job-name="${JobID1}" ${script_path}/folder_setup_ex.sh "${In}" "${Out}" "${Job}" )

srun --job-name="${JobID1}" ${script_path}/folder_setup_ex.sh "${In}" "${Out}" "${Job}" 

#echo "setting up folder system"
#sleep 2
#echo ${job1}
#j1=${job1}

wait

no_of_lines=$(wc -l < $2/samplesheet.txt)
echo "${no_of_lines}"

job2=$(sbatch --parsable  --job-name="${JobID2}" --array=1-"$((${no_of_lines}*2))"%24 ${script_path}/fastqc_ex.sh "${Out}" "${Job}")

sleep 2
echo ${job2}
j2=${job2} 

job3=$(sbatch --parsable  --job-name="${JobID3}" --array=1-"${no_of_lines}"%24 ${script_path}/bwa_duplicate_removal.sh "${Out}" "${Job}")

sleep 2
echo ${job3}
j3=${job3} 

job4=$(sbatch --parsable --dependency=aftercorr:$j3 --job-name="${JobID4}" --array=1-"${no_of_lines}"%24 ${script_path}/GATK_gVCF.sh "${Out}" "${Job}")

sleep 2
echo ${job4}
j4=${job4} 

job5=$(sbatch --parsable --dependency=aftercorr:$j4 --job-name="${JobID5}"  ${script_path}/combine_gVCF.sh "${Out}" "${Job}")

sleep 2
echo ${job5}
j5=${job5}

job6=$(sbatch --parsable --dependency=aftercorr:$j5 --job-name="${JobID6}" ${script_path}/vep_ex.sh "${Out}" "${Job}" "${script_path}")

sleep 2
echo ${job6}

job7=$(sbatch --parsable --dependency=aftercorr:$j3 --job-name="${JobID7}" ${script_path}/report_ex.sh "${Out}" "${Job}")
sleep 2
echo ${job7}

job8=$(sbatch --parsable --dependency=aftercorr:$j3 --job-name="${JobID8}"  ${script_path}/mito_ex.sh "${Out}" "${Job}")
sleep 2
echo ${job8}
j8=${job8}

job9=$(sbatch --parsable --dependency=aftercorr:$j8 --job-name="${JobID9}"  ${script_path}/mito_annot_ex.sh "${Out}" "${Job}" "${script_path}" )

sleep 2
echo ${job9}
