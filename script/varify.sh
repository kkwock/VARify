#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --mem=160GB
#SBATCH --job-name=varify
#SBATCH --output="%x-%j.out"
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=kk4764@nyu.edu

module purge

mp=/scratch/kk4764/varify/data/ENCFF283TLK_chr12.mpileup
snpeff=/scratch/kk4764/VARify/data/test/ENCFF541HLI.snp_effects.chr12_only_filtergenic_snps.no_intron.no_utr.tsv

module load anaconda3/2020.07

#python3 '/scratch/kk4764/VARify/test.py' --in_pileup $mp --in_snpeff $snpeff
python3 '/scratch/kk4764/VARify/test.py'
