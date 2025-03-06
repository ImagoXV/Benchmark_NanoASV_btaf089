#!/bin/bash
#SBATCH -p unlimitq
#SBATCH -t 12:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=24

source ~/work/miniconda3/bin/activate

conda activate NanoASV

DATABASE="/home/acousson/work/NanoASV/Benchmark/NR99_ref/NanoASV_format_silva_nr99_$

/usr/bin/time -v bash ${NANOASV_PATH}/workflow/run.sh --dir ~/work/Data/Nygaard/ --$


