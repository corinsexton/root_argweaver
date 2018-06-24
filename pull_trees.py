#!/usr/bin/python

import sys, re, gzip

#364935778134):[&error=6.045031813464135E-4, support=0.9938888888888889]0.001345974672911322

SMC_file = gzip.open(sys.argv[1])

names = SMC_file.readline().strip().split('\t')[1:]

outfile_label = sys.argv[1].replace('smc.gz','')

namesDict = {}
count = 0
for name in names:
	namesDict[str(count) + ':'] = name.strip().replace(':','_')
	count += 1

for line in SMC_file:
	if "TREE" in line:
		#TREE    1       54      ((((59:1545.314509[&&NHX:age=0.000000],412:1545.314509[&&NHX:age=0.000000])824:16499.3	
		ll = line.strip().split('\t')
		tree = ll[3]
		outfile = open(outfile_label + ll[1] + '_' + ll[2] +'.newick', 'w')
		for num in namesDict:
			name = namesDict[num]
			if '(' + num in tree:
				tree = tree.replace('(' + num, '(' + name + ':')
			elif ')' + num in line:
				tree = tree.replace(')' + num, ')' + name + ':')
			elif ',' + num in line:
				tree= tree.replace(',' + num, ',' + name + ':')

		outfile.write(tree + '\n')
		outfile.close()
