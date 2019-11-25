# -*- coding: cp1253 -*-
import random
import sys
from Bio import SeqIO


for seq_record1 in SeqIO.parse("Lysozyme-c.fasta", "fasta"):
	chromosome1 = seq_record1.seq
	
for seq_record2 in SeqIO.parse("á-lactalbumin.fasta", "fasta"):
	chromosome2 = seq_record2.seq
	
k = random.randrange(1,10)
print "k is ", k, "\n"
n=['A', 'C', 'G', 'T']
dna1=[]
dna2=[]

x1=0
y1=random.randrange(5,30)
while x1 < y1:
    rdna1 = random.choice(n)
    dna1.append(rdna1)
    x1=x1+1

x2=0
if y1-k>0:
    y2=y1-k
else:
    y2=y1-1

y2=random.randrange(5,30)   
while x2 < y2:
    rdna2 = random.choice(n)
    dna2.append(rdna2)
    x2=x2+1

dna1 = seq_record1.seq
dna2 = seq_record2.seq

n=len(dna1) 
m=len(dna2) 
print dna1
print dna2, "\n"

#D = [[0 for x in xrange(n)] for y in xrange(m)]
D=[]
for i in range(m+1):
    D.append([])
    for j in range(n+1):
        D[i].append(0)

        
#loops that print score matrix properly
print "   ",
for j in range(n):
    print dna1[j],
print
print " ",
for j in range(n+1):
    print D[0][j],
print
for i in range(1,m+1):
    print dna2[i-1],
    for j in range(n+1):
        print D[i][j],
    print
print "\n"

    
for j in range(1, n+1):
  D[0][j] = (-2)*j

for i in range(1,m+1):
  D[i][0] = (-2)*i


#loops that print score matrix properly
print "   ",
for j in range(n):
    print dna1[j],
print
print " ",
for j in range(n+1):
    print D[0][j],
print
for i in range(1,m+1):
    print dna2[i-1],
    for j in range(n+1):
        print D[i][j],
    print
print "\n"



sim_mat = [ 
[2, -1, 1, -1], 
[-1, 2, -1, 1],
[1, -1, 2, -1],
[-1, 1, -1, 2]]

n1 =0
n2 =0   
for i in range(1, m+1):
    for j in range(1, n+1):
        if dna1[j-1] == "A":
            n1 =0
        elif dna1[j-1] == "C":
            n1 =1
        elif dna1[j-1] == "G":
            n1 =2
        elif dna1[j-1] == "T":
            n1 =3
        else:
            n1=-1

        if dna2[i-1] == "A":
            n2 =0
        elif dna2[i-1] == "C":
            n2 =1
        elif dna2[i-1] == "G":
            n2 =2
        elif dna2[i-1] == "T":
            n2 =3
        else:
            n2=-1
            
        match = D[i-1][j-1] + sim_mat[n1][n2]
        gapdna1 = D[i][j-1] -2
        gapdna2 = D[i-1][j] -2
        D[i][j] = max(match,gapdna1,gapdna2)

#loops that print score matrix properly
print "   ",
for j in range(n):
    print dna1[j],
print
print " ",
for j in range(n+1):
    print D[0][j],
print
for i in range(1,m+1):
    print dna2[i-1],
    for j in range(n+1):
        print D[i][j],
    print
print "\n"

x=m
y=n

dna2_aln=""
dna1_aln=""

while x > 0 or y > 0:
    if dna1[y-1] == "A":
        n1 =0
    elif dna1[y-1] == "C":
        n1 =1
    elif dna1[y-1] == "G":
        n1 =2
    elif dna1[y-1] == "T":
        n1 =3
    else:
        n1=-1

    if dna2[x-1] == "A":
        n2 =0
    elif dna2[x-1] == "C":
        n2 =1
    elif dna2[x-1] == "G":
        n2 =2
    elif dna2[x-1] == "T":
        n2 =3
    else:
        n2=-1
        
    if D[x][y] - sim_mat[n1][n2] == D[x-1][y-1]:
        dna2_aln = dna2_aln + dna2[x-1]
        dna1_aln = dna1_aln + dna1[y-1]
        x = x-1
        y = y-1
    elif D[x][y] - (-2) == D[x][y-1]:
        dna1_aln=dna1_aln+dna1[y-1]
        dna2_aln=dna2_aln + "_"
        y = y-1
    elif D[x][y] - (-2) == D[x-1][y]:
        dna1_aln = dna1_aln + "_"
        dna2_aln = dna2_aln + dna2[x-1]
        x = x-1
    else:
        print('should not happen')


if y > 1:  
  while y > 1:
    dna1_aln = dna1_aln + dna1[y-1]
    dna2_aln=dna2_aln + "_"
    y = y-1
elif x > 1: 
  while x > 1:
    dna1_aln = dna1_aln + "_" 
    dna2_aln = dna2_aln + dna2[x-1] 
    x = x-1

dna1_aln = dna1_aln[::-1]
dna2_aln = dna2_aln[::-1]
print dna1_aln
print dna2_aln


