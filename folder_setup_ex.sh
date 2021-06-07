#! /bin/bash
#SBATCH -A rtmngs
#SBATCH -p bigmem
#SBATCH --mem=10G


In=$1
Out=$2
Job=$3

cd $Out
mkdir ${Job}_logs ${Job}_report ${Job}_results ${Job}_rawfiles
mkdir ${Job}_report/${Job}_priority
mkdir ${Job}_report/mito

#Making separate folders for the results of the analysis
fastq_folder=${Out}/${Job}_results/${Job}_fastqc
mkdir ${fastq_folder}

bwa_folder=${Out}/${Job}_results/${Job}_bam
mkdir ${bwa_folder}

#mkdir ${bwa_folder}/mito_coverage_exons
mkdir ${bwa_folder}/mito_coverage_CDS

vcf_folder=${Out}/${Job}_results/${Job}_vcf
mkdir ${vcf_folder}
mkdir ${vcf_folder}/vcf_pivot

echo "Setup of the folders is done"


##Sample concatenation; assumes that if the files are ending in _L*_R[1-2]_001.fastq.gz, then concatenation is required, if the ending is "*_R[1-2].fastq.gz", then concatenation is not reqiored
cd $In

if ls $In/*R1.fastq.gz 2> /dev/null
then
    echo 'No concatentation needed, files are ending in R*.fastq.gz'
    cp $In/*fastq.gz ${Out}/${Job}_rawfiles
else 
    echo "Concatenation required"
    ls -1 *R*_001.fastq.gz | sed 's/_L[0-9]*_R[1-2]_001.fastq.gz//' | uniq > ${Out}/${Job}_rawfiles/fastq_combine.txt

    samplesheet=${Out}/${Job}_rawfiles/fastq_combine.txt
    END=$(wc -l < ${samplesheet})

    for i in $(seq 1 $END)
    do
    sample=$(sed -n "$i"p $samplesheet |  awk '{print $1}') ;
    cat $(ls ${sample}_L*_R1_001.fastq.gz) >  ${Out}/${Job}_rawfiles/${sample}_R1.fastq.gz
    cat $(ls ${sample}_L*_R2_001.fastq.gz) >  ${Out}/${Job}_rawfiles/${sample}_R2.fastq.gz
    done

    rm ${Out}/${Job}_rawfiles/fastq_combine.txt

    echo "The fastq files are concatenated and copied to the working directory"

fi


##Make a samplesheet matrix for bwa and fastqc

cd ${Out}/${Job}_rawfiles

(ls -d $PWD/*.fastq.gz) > ${Out}/sample_list.txt

for f in `ls -1 *_R1.fastq.gz | sed 's/_R1.fastq.gz//' `
do
echo ${f} ${f}_R1.fastq.gz ${f}_R2.fastq.gz >> ${Out}/samples.txt
done

#To prevent the duplicates when the pipeline is rerun
sort ${Out}/samples.txt | uniq > ${Out}/samplesheet.txt
rm ${Out}/samples.txt

echo "Sample sheet generation is done"


##test