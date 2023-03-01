#!/bin/bash
#
cd $PBS_O_WORKDIR
snakemake --jobname 's.{jobid}.{rulename}' -k --stats snakemake.stats --rerun-incomplete --restart-times 4 -j 300 --cluster 'qsub {params.batch}'  >&  snakemake.log
