#!/bin/bash
#SBATCH --ntasks-per-node=2
#SBATCH --time=1:00:00

# Run the analysis
#wait 1

Rscript ${SCRIPTPATH} ${QUANT} ${META} ${WORK} ${LEN} ${DIST} '0' ${CUT} ${INTER} ${TISSUES}

