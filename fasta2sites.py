#!/usr/bin/python

import sys
from collections import defaultdict

# >HLA:HLA00132_B_07:02:01
# GATCAGGACGAAGTCCCAGGTCCCGGACGGGGCTCTCAGGGTC

infasta = open(sys.argv[1],'r')

sitesDict = defaultdict(list)
line = infasta.readline()
outfile = open(sys.argv[2], 'w')

labelList = []
while line != '':
	count = 1
	labelList.append(line.rstrip('\n').lstrip('>'))
	seq = infasta.readline().rstrip('\n')
	for char in seq:
		sitesDict[count].append(char)
		count += 1
	line = infasta.readline()

outfile.write("NAMES")
for label in labelList:
	outfile.write('\t' + label)

outfile.write('\n')
outfile.write("REGION\tchr\t1\t" + str(count-1) + '\n')

for site in sorted(sitesDict.keys()):
	siteSet = set(sitesDict[site])
	siteSet.discard('N')
	if siteSet and len(siteSet) != 1:
		outfile.write(str(site) + '\t' + ''.join(sitesDict[site]) + '\n')


	

	
