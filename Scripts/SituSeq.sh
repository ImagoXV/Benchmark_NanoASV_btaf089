#!/bin/bash
#SBATCH -p unlimitq
#SBATCH -t 24:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=24

module load statistics/R/4.4.0

/usr/bin/time -v Rscript ~/work/NanoASV/Benchmark/SituSeq/scripts/setup.r

/usr/bin/time -v Rscript ~/work/NanoASV/Benchmark/SituSeq/scripts/preprocessing.r

/usr/bin/time -v Rscript ~/work/NanoASV/Benchmark/SituSeq/scripts/Stream-1A.r

/usr/bin/time -v Rscript ~/work/NanoASV/Benchmark/SituSeq/scripts/Stream-1B.r

/usr/bin/time -v Rscript ~/work/NanoASV/Benchmark/SituSeq/scripts/Phylosequization.r
