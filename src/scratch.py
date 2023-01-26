# Libraries
import sys, re

# read in the test fasta from data
g = '../VARify/data/test/test_gene.fa'
v = '../VARify/data/test/test_var.fa'

#f=open('genes.fasta','r')

def geneParse(x):
        f = open(x, 'r')
        lines=f.readlines()

        hre=re.compile('>(\S+)')
        lre=re.compile('^(\S+)$')

        gene={}

        for line in lines:
                outh = hre.search(line)
                if outh:
                        id=outh.group(1)
                else:
                        outl=lre.search(line)
                        if(id in gene.keys()):
                                gene[id] += outl.group(1)
                        else:
                                gene[id]  =outl.group(1)
        return gene

# test seq match
# seq = 'ATC----GCTGTA--A-----T'
# matches = list(re.finditer('-+', seq))
#
# print 'Number of gaps =', len(matches)
# print
#
# for region_number, match in enumerate(matches, 1):
#     print 'Index Position of Gap region {} = {} to {}'.format(
#             region_number,
#             match.start(),
#             match.end() - 1)
#     print 'Length of Gap region {} = {}'.format(
#             region_number,
#             match.end() - match.start())
#     print


# Test
gene = geneParse(g)
var = geneParse(v)
# Quick Comparison
tuple(gene.values()) == tuple(var.values())

# Plan
# 1. Get FASTA
# 2. Align sequences --> Create VCF
# 3. Analyze Differences

from Bio import AlignIO