#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=8:00:00
#SBATCH --mem=8GB
#SBATCH --job-name=download_data
#SBATCH --output="%x-%j.out"
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=kk4764@nyu.edu

# WGS
wget --timestamping ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/WGS/WGS_EA_N_1.bwa.dedup.bam -P .
wget --timestamping ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/WGS/WGS_EA_N_1.bwa.dedup.bai -P .

# WES
wget --timestamping ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/WES/WES_EA_N_1.bwa.dedup.bam -P .
wget --timestamping ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/WES/WES_EA_N_1.bwa.dedup.bai -P .
