#! python3

import re #needed for regular expressions
import sys #needed for arguments given to the script
import os #needed to read all the protein files given in the directory
import subprocess

#This is a script that assists IdentiPy in iterating over a directory with input files and analysing one after another.
#The input is the directory in which all the mgf files are located.
#For the output this script will create a directory that should be given as an argument to this script.
#This script also needs the location of the six-frame translation of course.

try: #trying to assign the location of the files to two variables and to open the parameter file
	protein_data = os.listdir(sys.argv[1]) #a list containing all the file names of the protein data in mgf file format
	path_to_files = sys.argv[1]
	path_to_6_frame = sys.argv[2] #path to 6 frame translation
	path_to_output = sys.argv[3] #path to the location of the output directory
except IndexError: #if one or both arguments are not given:
	print("Please enter the location of the mgf data directory and the location of the 6 frame translation and the path to the output directory.")
except FileNotFoundError:
	print("The raw data directory cannot be found. Please check its existence.")
else: #if all arguments are given and the parameter can be found:

	
	if(os.path.exists(path_to_output)):
		rm_directory = subprocess.run("rm -r " + path_to_output, shell = True)
	command1 = 'mkdir ' + path_to_output #this a bash command
	create_directory=subprocess.run(command1,shell=True) #creating the sub.process that runs the bash command

	for x in protein_data: #iterating over all found files
		filename = path_to_files + x #appending the path to the file name
		print(filename + "\n")
		#a,b = x.split(".") #splitting x by the . which splits the name from the file type (raw)
		command = 'identipy ' + filename + ' -out ' + path_to_output + ' -db ' + path_to_6_frame #this is the bash command
		run_ip=subprocess.run(command,shell=True) #creating the sub.process that runs the bash command
