#! python3


#Author: Leonie Johanna Lorenz
#Last modified: 19th September 2019

### This is a script for parsing an pep.xml (output of IdentiPy) to a GFF format that can be used by AUGUSTUS.

import re #needed for regular expressions
import sys #needed for arguments given to the script

try:
	file_handle_1 = open(sys.argv[1],"r")#open file with six frame translation
	file_handle_2 = open(sys.argv[2],"r") #Output of Identipy, list of matched peptides
except IndexError: #not both needed arguments are given
	print("Please give the following two arguments: 1) the file containing the six-frame-translation of your genome, 2) the location of pep.xml, an output file of IdentiPy, 3) the output file (IdentiPy's hints in GFF format) and if you like 4) the percentage at which you would like to filter (quality of prediction). Default will be 0.1%.")
except FileNotFoundError: #one or both file could not be found where indicated by the user
	print("One or both files you entered do not exist.")
else: #if both input files where found and opened without any errors:
	if(len(sys.argv)>=4): #if there is a direction to the output file
		file_handle_3 = open(sys.argv[3],"w")
	else: #create default output file
		a,b,c = sys.argv[2].split(".") #splitting path to pep xml by .
		file_handle_3 = open(a + ".gff","w") #the file I am creating here
	quality_filter = 0.001 #default filter
	if(len(sys.argv)>=5): #if a filter is wanted
		quality_filter=float(sys.argv[4])

	#In this approach I will iterate over the 6 frame translation first and then over the identipy results.
	
	cds_dict_strand = {} #dictionairy to know which CDS is on which strand
	cds_dict_start = {} #for plus strand first number, for minus second number 
	cds_dict_seq = {} #dictionairy for the sequences of all CDSs
	cds_name = "" #name of the current cds, e.g. CDS1
	new_cds = True #True if a new cds is starting
	current_cds = "" #seq of the current cds
	for x in file_handle_1: #6-frame translation
		name_match = re.search(r"^>(.+)\|(CDS\d+)\slength.*from\s([a-z\(\s]*)(\d+)\.\.(\d+)[\s\)]",x) #matches the descriptive part
		seq_match = re.search(r"^([A-Z]+)\n",x) #matches all sequence lines
		if(name_match):
			new_cds = True #found the beginning of a new CDS
			cds_dict_seq[cds_name] = current_cds #set the sequence of the PREVIOUS match to the CDS found before 
			cds_name = name_match.group(2) #set the variable cds_name to the current CDS
			if(name_match.group(3)=="complement("): #minus strand
				cds_dict_strand[name_match.group(2)]="-"
				cds_dict_start[name_match.group(2)]=name_match.group(5)
			else: #plus strand
				cds_dict_strand[name_match.group(2)]="+"
				cds_dict_start[name_match.group(2)]=name_match.group(4)
		if(seq_match): #matched a line containing a part of a CDS
			if(new_cds): #if it is the first part of a CDS
				current_cds = seq_match.group(1) #set it to the currently found part
				new_cds = False #not anymore at the beginning of the CDS
			else: #not the first part of the CDS
				current_cds = current_cds + seq_match.group(1) #append current part
	cds_dict_seq[cds_name] = current_cds #set the sequence of the PREVIOUS match to the CDS found before 
#entering file_handle_2:
	first_pep_dict = {} #dictionairy that saves for each CDS the most upstream lying peptide with a good quality
	current_line = ""
	for x in file_handle_2: #iterating over IdentiPy's peptides
		pep_match = re.search(r'<search_hit.+peptide="([A-Z]+)".+protein="(.+)\|(CDS\d+)".+',x) #matches the desctiption and the sequence of the peptides
		evalue_match = re.search(r'"expect"\svalue="(.+)"',x) #matches the e-value of a peptide
		if(pep_match):
			current_pep = pep_match.group(3) #name of the cds
			current_seq=pep_match.group(1) #seq of the current peptide
			local_index = cds_dict_seq[current_pep].find(current_seq) #finds out at which index of the CDS the peptide starts
			current_chr = pep_match.group(2) #saves which on which chromosome the peptide was found on
			if((cds_dict_strand.get(current_pep,-1))==("+")): #plus strand
				current_line = pep_match.group(2)+"\t"+"Proteomics"+"\t"+"CDSpart"+"\t"+str(int(cds_dict_start[current_pep])+(int(local_index)*3))+"\t"+str(int(cds_dict_start[current_pep])+(int(local_index)*3)+((len(current_seq))*3)-1)+"\t"+"."+"\t"+"+"+"\t"+".\t"+"src=P;len=" + str((len(current_seq))*3) +"\n" #saves the output of the currently found peptide
				
			elif((cds_dict_strand.get(current_pep,-1))==("-")): #minus strand
				current_line = pep_match.group(2)+"\t"+"Proteomics"+"\t"+"CDSpart"+"\t"+str(int(cds_dict_start[current_pep])-(int(local_index)*3)-((len(current_seq))*3)+1)+"\t"+str(int(cds_dict_start[current_pep])-(int(local_index)*3))+"\t"+"."+"\t"+"-"+"\t"+".\t"+"src=P;len=" + str((len(current_seq))*3) +"\n" #saves what should be written into the output file for the current peptide

			else: #meaning =-1
				print("There are peptides with CDS names not corresponding to any CDSs in the 6 frame translation. Please make sure that the 6 frame translation file you entered was the one used for creating the IdentiPy hints.\n")

		if(evalue_match): #if it matched a line containing an e-value
			if(float(evalue_match.group(1))<quality_filter): #if the e-value is below the chosen threshold
				file_handle_3.write(current_line) #write the information of the current peptide into the output file
				if((cds_dict_strand.get(current_pep,-1))==("+")): #if plus strand
					if((first_pep_dict.get(current_pep,[-1])[0]==(-1))or (first_pep_dict.get(current_pep,[-1])[0]>int(local_index))): #if for this CDS there hasn't been found a peptide yet or if the currently found one is more upstream
						if(current_seq[0]=="M"):  #if the current peptide starts with Methionin
							first_pep_dict[current_pep]=[int(local_index),int(cds_dict_start[current_pep])+(int(local_index)*3),"+",current_chr,1] #save the one which has its start the nearest to the beginning of the ORF
						else: #if it does not start with Methionin
							first_pep_dict[current_pep]=[int(local_index),0,"+",current_chr,1] #save 0 as the stop site
					elif(first_pep_dict.get(current_pep,[-1])[0]==int(local_index)): #if the just found start hint is equal to the one already noted in the dictionairy
						if(current_seq[0]=="M"): #only if it really is a start hint, otherwise multiplicity is of no interest
							first_pep_dict[current_pep][4] += 1 #increase multiplicity by one
				else: #minus strand
					if((first_pep_dict.get(current_pep,[-1])[0]==(-1))or (first_pep_dict.get(current_pep,[-1])[0]>int(local_index))): #if for this CDS there hasn't been found a peptide yet or if the currently found one is more upstream
						if(current_seq[0]=="M"): #if it starts with Methionin
							first_pep_dict[current_pep]=[int(local_index),int(cds_dict_start[current_pep])-(int(local_index)*3),"-",current_chr,1] #save the one which has its start the nearest to the beginning of the ORF
						else: #if it does not start with Methionin
							first_pep_dict[current_pep]=[int(local_index),0,"-",current_chr,1]
					elif(first_pep_dict.get(current_pep,[-1])[0]==int(local_index)):#if the just found start hint is equal to the one already noted in the dictionairy
						if(current_seq[0]=="M"):
							first_pep_dict[current_pep][4] += 1 



	for x in first_pep_dict:#iterating over the peptide dictionairy 
		if(first_pep_dict[x][1]!=0): #if a peptide starts with Methionin
			if(first_pep_dict[x][2]=="+"): #if plus strand
				file_handle_3.write(first_pep_dict[x][3] + "\t" + "Proteomics" + "\t" + "start" + "\t" + str(first_pep_dict[x][1]) + "\t" + str(first_pep_dict[x][1]+2) + "\t.\t" + "+" + "\t.\t" + "src=P;len=3;mult=" + str(first_pep_dict[x][4]) + "\n")
			else: #if minus strand
				file_handle_3.write(first_pep_dict[x][3] + "\t" + "Proteomics" + "\t" + "start" + "\t" + str(first_pep_dict[x][1]-2) + "\t" + str(first_pep_dict[x][1]) + "\t.\t" + "-" + "\t.\t" + "src=P;len=3;mult="+ str(first_pep_dict[x][4])+ "\n")

file_handle_1.close()
file_handle_2.close()
file_handle_3.close()
