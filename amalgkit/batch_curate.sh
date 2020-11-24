#!/bin/bash
#SBATCH --ntasks-per-node=2
#SBATCH --time=1:00:00
#SBATCH --mem=5GB

# Run the analysis
Rscript ${SCRIPTPATH} ${QUANT} ${META} ${WORK} ${LEN} ${DIST} '0' ${CUT} ${INTER} ${TISSUES}

