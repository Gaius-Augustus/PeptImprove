#!python3

#Author: Leonie Johanna Lorenz
#Last modified: 19th September 2019

### This is a script for parsing a gff3 file (output of AUGUSTUS) to a GFF format that can be used to compare and combine GeneMark's and AUGUSTUS' outputs.

import re #needed for regular expressions
import sys #needed for arguments given to the script

try:
	file_handle_1 = open(sys.argv[1],"r")#open the gff3 file
	file_handle_2 = open(sys.argv[2],"w") #direction to the output file
except IndexError: #the arguement needed ist not given
	print("Please give the following arguments: a gff3 file produced by AUGUSTUS and the name of the output file.")
except FileNotFoundError: #the file could not be found where indicated by the user
	print("The file you entered does not exist.")
else: #if the input file was found and opened without any errors:
	
	for x in file_handle_1:	#iterating over the lines of the gff3 file
		matching = re.search(r"^(.*)\t(AUGUSTUS)\t(gene)\t(\d+)\t(\d+)\t.+\t([+-])\t.+\t(.*)\n",x) #matching all lines with gene predictions
		if(matching):
			file_handle_2.write(matching.group(1) + "\t" + matching.group(2) + "\t" + matching.group(3) + "\t" + matching.group(4) + "\t" + matching.group(5) + "\t" + "." + "\t" + matching.group(6) + "\t" + "." + "\t" + matching.group(7) + "\n") #output in gff format

file_handle_1.close()
file_handle_2.close()
