# Libraries
import sys, re

# read in the test fasta from data
g = '../VARify/data/test/test_gene.fa'
v = '../VARify/data/test/test_var.fa'

f=open(g,'r')
gene =f.readlines()
f=open(v,'r')
var = f.readlines()

# compare the two

# capture difference