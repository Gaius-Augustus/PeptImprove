#! python3

import re #needed for regular expressions
import sys #needed for arguments given to the script

# This is a script to collapse redundant hints of IdentiPy into one single line with a remark of the multiplicity of that hint.
# This script will be used instead of sort -u in bash but after sort was done.


try:
	ip_hints_red = open(sys.argv[1],"r")#open file with sorted but redundant ip hints
except IndexError: #not both needed arguments are given
	print("Please give the following two arguments: 1) the file containing the sorted but still redundant IdentiPy hints and 2) the output file into which the sorted and nonredundant hints should be written.")
except FileNotFoundError: #one or both file could not be found where indicated by the user
	print("One or both files you entered do not exist.")
else:
	ip_hints_nonred = open(sys.argv[2],"w") #output file for nonredundant hints

	previous_line = " "
	current_mult = 1
	for x in ip_hints_red:
		type_match = re.search(r"CDSpart",x)
		line_match = re.search(r"^(.+)\n",x)
		if(type_match):
			if(previous_line == line_match.group(1)): #current line is equal to the previous line
				current_mult += 1
			elif(previous_line != " "): #the previous line was not empty (== line zero)
				ip_hints_nonred.write(previous_line + ";mult=" + str(current_mult) + "\n")
				previous_line = line_match.group(1)
				current_mult = 1
			else: #the current line is the very first line of the file
				previous_line = line_match.group(1)
		else:
			ip_hints_nonred.write(x)
	ip_hints_nonred.write(previous_line  + ";mult=" + str(current_mult) + "\n")

ip_hints_red.close()
ip_hints_nonred.close()
