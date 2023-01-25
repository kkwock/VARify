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


# Test
gene = geneParse(g)
var = geneParse(v)
# Quick Comparison
tuple(gene.values()) == tuple(var.values())

# compare the two

# capture difference