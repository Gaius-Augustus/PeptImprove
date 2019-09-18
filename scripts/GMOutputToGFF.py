#!python3

import re #needed for regular expressions
import sys #needed for arguments given to the script

### This is a script for parsing a nucleotide list (output of GeneMark) to a GFF format that can be used by an Augustus script to create a genbank flat file format.

try:
	file_handle_1 = open(sys.argv[1],"r")#open nucleotides.lst
except IndexError: #the arguement needed ist not given
	print("Please give the following argument: the list of nucleotides produced by GeneMark.")
except FileNotFoundError: #one or both file could not be found where indicated by the user
	print("The file you entered does not exist.")
else: #if the input file was found and opened without any errors:
	if(len(sys.argv)==3): #if there is another argument
		file_handle_2 = open(sys.argv[2],"w") #if there is a direction to the output file
	else: #default direction to the output file
		file_handle_2 = open("GffFromGM.gff","w") #the file I am creating here (default)
	transcript_id = 0 #stores the number of the current transcript
	for x in file_handle_1: #iterates over the nucleotide list
		matching = re.search(r"^>\d+\s(.+)\s(\d+)\s(\d+)\s([+-])\s",x)
		if(matching):
			transcript_id +=1 #id is increased by one
			file_handle_2.write(matching.group(1) + "\t" + "GeneMarkS2"  + "\t" + "CDS" + "\t" + matching.group(2) + "\t" + matching.group(3) + "\t" + "." + "\t" + matching.group(4) + "\t" + "0" + "\t" + 'gene_id "g' + str(transcript_id) +  '.t1"; transcript_id "g' + str(transcript_id) + '";' + "\n") #output is written in GFF format

file_handle_1.close()
file_handle_2.close()
