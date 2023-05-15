#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --time=48:00:00
#SBATCH --mem=150GB
#SBATCH --job-name=pileup
#SBATCH --output="%x-%j.out"
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=kk4764@nyu.edu

module purge

module load samtools/intel/1.14

dir=/scratch/kk4764/varify/data/
bam=ENCFF283TLK.bam
ref=/scratch/kk4764/varify/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
#samtools flagstat ${dir}${bam} > ${dir}flag.txt

TMPDIR=/scratch/kk4764/varify/data/
cd $TMPDIR

samtools mpileup -f ${ref} -s ${dir}${bam} --output mpileup_out.txt --min-MQ 21 --min-BQ 21 --max-depth 4000