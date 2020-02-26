''' Adds missing fields to list of taxonomy strings. Use: python fix_taxonomy_L7.py my_taxa.txt > new_taxa.txt '''

import sys

with open(sys.argv[1], "r") as f:
	for line in f:
		x = line.rstrip("\n")
		y = x.count(";")
		if y < 6:
			add = 6 - y
			new = (";" + x.split(";")[-1] + "_unknown") * add
			print(x + new)
		else:
			print(x)
