#!/usr/bin/env bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 8
#SBATCH --mem 40G 
#SBATCH --time 0-2:00:00
#SBATCH --partition norm 

module load cellranger

cellranger count --id test2 \
    --fastqs data/GEX_library/ \
    --transcriptome=$CELLRANGER_REF/refdata-gex-GRCh38-2020-A \
    --localcores=8 \
    --chemistry=auto \
    --localmem=36