#!/bin/bash
#SBATCH -p unlimitq
#SBATCH -t 72:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=24

#module load bioinfo/Porechop/0.2.4


DATABASE="/home/acousson/work/NanoASV/Benchmark/NR99_ref/NanoASV_format_silva_nr99_v138.1_train_set.fa"

DATASET="/home/acousson/work/NanoASV/Benchmark/Nygaard/Dataset"

#Filtering
#for barcode in ${DATASET}/*.fastq.gz; do
#   /usr/bin/time -v porechop --verbosity 0 -i ${barcode} -o "choped_$(basename "$barcode")"
#	echo $barcode
#done

#module purge

module load devel/Miniconda/Miniconda3
module load bioinfo/chopper/0.7.0
#Choping

for barcode in ./choped*.fastq.gz; do
    /usr/bin/time -v zcat $barcode | chopper -l 1300 --maxlength 1700 -q 9 > filtered_$(basename ${barcode}).fastq
done

module purge

module load bioinfo/LAST/1542

#Index
#/usr/bin/time -v lastdb Nr_test_99.db ${DATABASE}

#Alignment
for barcode in ./filtered*.fasta; do
   /usr/bin/time -v  lastal -v Nr_test_99.db ${barcode} -f BlastTab -P 8 -r 1 -q 1 -a 1 -b 1 > ${barcode}.blast.tab
done


