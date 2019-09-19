#! python3


#Author: Leonie Johanna Lorenz
#Last modified: 19th September 2019

#This is a script that assists the ThermoRawFileParser in converting more than one file.
#The input is the directory in which all the raw files that need to be converted are located.
#For the output this script will create a directory in the raw file directoy called "mzml".
#You will find your mzml and metadata.txt files there.

import re #needed for regular expressions
import sys #needed for arguments given to the script
import os #needed to read all the protein files given in the directory
import subprocess


try: #trying to assign the location of the files to two variables and to open the parameter file
	protein_data = os.listdir(sys.argv[1]) #a list containing all the file names of the protein data in raw file format
	path_to_files = sys.argv[1]
	path_to_parser = sys.argv[2]
	path_to_output = sys.argv[3]
except IndexError: #if one or both arguments are not given:
	print("Please enter the location of the raw data directory.")
except FileNotFoundError:
	print("The raw data directory cannot be found. Please check its existence.")
else: #if all arguments are given and the parameter can be found:

	for x in protein_data: #iterating over all found files
		filename = path_to_files + x #appending the path to the file name
		a,b = x.split(".") #splitting x by the . which splits the name from the file type (raw)
		command = 'mono ' + path_to_parser + '/ThermoRawFileParser.exe -i=' + filename + ' -b=' +path_to_output + "/" + a + '.mgf -f=0 -m=1' #this is the bash command (for mzml: .mzml and -f=1)
		correct=subprocess.run(command,shell=True) #creating the sub.process that runs the bash command
