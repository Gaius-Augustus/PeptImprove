#!/usr/bin/env python

#Author: Leonie Johanna Lorenz
#Last modified: 19th September 2019

#This is a script that combines the output of GeneMarkS-2 and AUGUSTUS.
#If GM and AUG differ in their predictions, AUGUSTUS' gene will be chosen if there is extrinsic evidence for this difference.
#If not, GM's gene will be chosen.

import os
import argparse #needed for parsing the arguments
import re #needed for regular expressions

parser = argparse.ArgumentParser(description='combines the output of GeneMarkS2 and AUGUSTUS depending on extrinsic evidence')
parser.add_argument('-gm','--genemark',type=str,required=True,
	help='Gene prediction of GeneMark in gff format (required)')
parser.add_argument('-aug','--augustus',type=str,required=True,
	help='Gene prediction of AUGUSTUS in gff format (required)')
parser.add_argument('-ext','--extrinsic_evidence',type=str,required=True,
	help='Hints based on MS/MS experiments in gff format')
parser.add_argument('-o','--output',type=str,required=True,
	help='location of the outut file')
args = parser.parse_args()

if(args.genemark):
	genemark_handle = open(args.genemark,"r")
if(args.augustus):
	augustus_handle = open(args.augustus,"r")
if(args.extrinsic_evidence):
	peptide_handle = open(args.extrinsic_evidence,"r")
if(args.output):
	output_handle = open(args.output,"w")

gm_dictionairy_plus = {} #dictionairy with GeneMark's predictions concerning the plus strand
gm_dictionairy_minus = {}#dictionairy with GeneMark's predictions concerning the minus strand
gm_list_start = [] #list of all first numbers predicted by GM (plus strand: start site, minus strand: stop site)
gm_list_stop = [] #list of all second number predicted by GM (pluss strand: stop site, minus strand: start site)
gm_gff_start = "" #string for saving the name of the chromosome
for x in genemark_handle:
	gm_match = re.search(r"^(.+\t\w+\t\w+\t)(\d+)\t(\d+)\t.\t([+-])\t",x) #matching all lines in the gff file
	if(gm_match):
		gm_gff_start = gm_match.group(1) 
		gm_list_start.append(gm_match.group(2)) #this is start for plus, stop for minus
		gm_list_stop.append(gm_match.group(3)) #this is stop for plus, start for minus
		if(gm_match.group(4)=="+"): # if it is a plus strand gene
			gm_dictionairy_plus[gm_match.group(3)]=gm_match.group(2) #key: stop site, value: start site
		else: #minus strand gene
			gm_dictionairy_minus[gm_match.group(2)]=gm_match.group(3) #key: stop site, value: start site
			
pep_list_plus_start = [] #this will be a list of all start sites of peptide evidence on the plus strand
pep_list_plus_stop = [] #this will be a list of all stop sites of peptide evidence on the plus strand
pep_list_minus_start = [] #this will be a list of all start sites of peptide evidence on the minus strand
pep_list_minus_stop = [] #this will be a list of all stop sites of peptide evidence on the minus strand
start_hints_dict_plus = {} # a dictionairy containing all start hints on the plus strand
start_hints_dict_minus ={} # a dictionairy containing all start hints on the minus strand
for x in peptide_handle:
	pep_match = re.search(r"^.+\t\w+\t(\w+)\t(\d+)\t(\d+)\t.\t([+-])\t",x) #matching all lines of the gff file
	if(pep_match):
		if(pep_match.group(4)=="+"): #plus strand
			if(pep_match.group(1)=="start"): #start hint
				start_hints_dict_plus[int(pep_match.group(2))]=pep_match.group(3) #add start hint to the dictionairy
			else: #normal hint (=no start hint)
				pep_list_plus_start.append(int(pep_match.group(2))) #append start of the hint to list of hints
				pep_list_plus_stop.append(int(pep_match.group(3))) #append end of the hint to list of hints
		else: #minus strand
			if(pep_match. group(1)=="start"): #start hint
				start_hints_dict_minus[int(pep_match.group(3))]=pep_match.group(2) #add start hint to the dictionairy
			else: #normal hint
				pep_list_minus_start.append(int(pep_match.group(3))) #append start of hint to list
				pep_list_minus_stop.append(int(pep_match.group(2))) #append stop of hint to list

aug_list_start = [] #list of the first number of each predicted gene, plus: start, minus: stop
aug_list_stop = [] #list of the second number of each predicted gene, plus: stop, minus:start
aug_list_strand = [] #list which saves for which strand the gene with this index was predicted
aug_gff_start = "" #saves the name of the chromosome
for x in augustus_handle:
	aug_match = re.search(r"^(.+\t\w+\t\w+\t)(\d+)\t(\d+)\t.\t([+-])\t",x) #matching all lines in the gff file
	if(aug_match):
		aug_gff_start = aug_match.group(1)
		aug_list_start.append(aug_match.group(2)) #appends number to list, plus: start, minus: stop
		aug_list_stop.append(aug_match.group(3)) #appends second number to list, plus: stop, minus:start
		if(aug_match.group(4)=="+"): #plus strand
			aug_list_strand.append("+") #appends strand
		else: #minus strand
			aug_list_strand.append("-") #appends strand

aug_position = -1 # the index of AUGUSTUS' lists
next = 0 #a variable for making sure that each of the genemark genes that is not found by Augustus will be printed into the output file
last_stop_plus = 0 #saves the position of the previous plus strand gene
last_stop_minus = 0 #saves the position of the previous minus strand gene
next_2=0 # an integer for printing the gene number into the last column of the output file
ext_count = 0 #counter for extrinsic evidence. only needed for new genes by AUG
for a in range(0,len(aug_list_start)):
	x = int(aug_list_start[a]) #x=start site if it is a plus gene, stop if it is a minus strand gen
	z = int(aug_list_stop[a]) # z=stop site if on plus strand, start site otherwise
	aug_position += 1 #increasing the position in AUGUSTUS' list by one at each iteration
	if((next<(len(gm_list_start))) and (last_stop_plus == gm_list_stop[next] or last_stop_minus == gm_list_start[next])): #if we haven't reached the end of gm_list_start and the next gene in the list was already written into the output file
		next += 1
	
	while((next<(len(gm_list_start)))and(int(gm_list_start[next])<int(x))and(int(gm_list_stop[next])<int(z))): #while we haven't reached the end of the list and both the start and the stop site have a lower index than the one we are looking at from AUGUSTUS
		if(gm_dictionairy_plus.get(gm_list_stop[next],-1)!=-1): #if the next gene is a plus strand gene
			output_handle.write(gm_gff_start+gm_list_start[next]+"\t"+gm_list_stop[next]+ "\t" + "." + "\t" +  "+" + "\t" + "." + "\t"+str(next_2+1)+"\n") #write it into the output file
		else: #next gene is on minus strand
			output_handle.write(gm_gff_start + gm_list_start[next]+"\t" + gm_list_stop[next]+"\t"+ "." + "\t" +  "-" + "\t" + "." + "\t"+str(next_2+1)+"\n") #write it into the output file
			
		next += 1 #gm's index is increased by one
		next_2 +=1 #gene number for the output file is increased by one
	ext_evidence = False #boolean variable for existence of extrinsic evidence
	print_gm = False #boolean variable whether GM's genes should be printed

	if(aug_list_strand[aug_position]=="+"): #if AUGUSTUS' next gene is on plus strand
		gm_start = int(gm_dictionairy_plus.get(aug_list_stop[aug_position],-1)) #returns -1 if the stop site does not exist in GM
			
		if(gm_start==x): #meaning: gene existis in GM AND starts at the same position as AUGUSTUS' gene
			output_handle.write(gm_gff_start + str(x) + "\t" + aug_list_stop[aug_position] +"\t" + "." + "\t" + aug_list_strand[aug_position] + "\t" + "." + "\t" + str(next_2+1) + "\n")
			next_2 += 1
			last_stop_plus = aug_list_stop[aug_position] #sets the variable for the previous gene to the current gene
		elif(gm_start!=-1): #meaning: GeneMark gene exists but does not have the same start	
			print_gm = True
			if(x>gm_start): #if AUGUSTUS' gene starts downstream of GM's
				if(start_hints_dict_plus.get(x,-1)!=(-1)): #if I have a start hint that supports this site, please keep
					ext_evidence = True
					print_gm = False	
				else: #no start hint to support AUGUSTUS' gene
					print_gm = True #print GM's gene instead
			else: #AUGUSTUS' gene start lies upstream of GM's
				if(start_hints_dict_plus.get(x,-1)!=(-1)): #checking for a start hint at this site
					ext_evidence = True
					print_gm = False
				else: #checking for some other hint
					for y in pep_list_plus_start: #iterating over all non-start hints
						if((x<=y<gm_start)and((y-x)%3==0)): #if there is a hint with its beginning lying between AUG's and GM's start sites and that is also in-frame
							ext_evidence = True
							print_gm = False
				
		else: #meaning: the searched STOP site does not exist in GeneMark
			ext_count = 0
			if(start_hints_dict_plus.get(x,-1)!=(-1)): #checking for a start hint at this site
					ext_count = 1		
#					ext_evidence = True
			pep_number = -1
			for y in pep_list_plus_start:
				pep_number += 1 #simply stores the position in the list
				if((x<y<int(aug_list_stop[aug_position]))and((y-x)%3==0)and(pep_list_plus_stop[pep_number]<=int(aug_list_stop[aug_position]))): #if there is a hint that begins behind the start site of AUG's gene and end in front of its stop codon and is in the same frame
					ext_count += 1	
			if(ext_count>=2):			
				ext_evidence = True
					#ext_evidence = True
#			for y in pep_list_plus_stop:
#				if((x<=y<int(aug_list_stop[aug_position]))and((int(aug_list_stop[aug_position])-y)%3==0)):
#					ext_evidence = True
			#if(ext_evidence): #if I found extrinsic evidence
			#	output_handle.write(aug_gff_start + str(x) + "\t" + aug_list_stop[aug_position] +"\t" + "." + "\t" + aug_list_strand[aug_position] + "\t" + "." + "\t" + str(next_2+1) + "\n")
			#	next_2 += 1			
			#ext_evidence = False
			
		if(print_gm): #if it was decided that GM's gene is more reliable than AUGUSTUS'
			output_handle.write(gm_gff_start + gm_dictionairy_plus[aug_list_stop[aug_position]] + "\t" + aug_list_stop[aug_position] + "\t"+ "." +"\t" + aug_list_strand[aug_position] + "\t" + "." + "\t" + str(next_2+1) + "\n")
			next_2 += 1
			last_stop_plus = aug_list_stop[aug_position]
		if(ext_evidence): #if I found extrinsic evidence for the cases needed
			output_handle.write(aug_gff_start + str(x) + "\t" + aug_list_stop[aug_position] +"\t" + "." + "\t" + aug_list_strand[aug_position] + "\t" + "." + "\t" + str(next_2+1) + "\n")
			next_2 += 1
			last_stop_plus = aug_list_stop[aug_position]

	else: #current gene is on the minus strand
		gm_start = int(gm_dictionairy_minus.get(aug_list_start[aug_position],-1)) #returns -1 if the stop site does not exist in gm
		if(gm_start==z): #meaning: gene existis in GM AND starts at the same position as AUGUSTUS' gene
			output_handle.write(gm_gff_start + str(x) + "\t" + aug_list_stop[aug_position] +"\t" + "." + "\t" + aug_list_strand[aug_position] + "\t" + "." + "\t" + str(next_2+1) + "\n")
			next_2 += 1
			last_stop_minus = aug_list_start[aug_position]
		elif(gm_start!=-1): #meaning: GM's gene exists but does not have the same start
			print_gm = True	
			if(z<gm_start): #AUGUSTUS' gene starts downstream of GM's (on the minus strand)
				if(start_hints_dict_minus.get(z,-1)!=(-1)): #if I have a start hint for z
					ext_evidence = True
					print_gm = False
			else: #AUGUSTUS' gene starts upstream of GM's
				if(start_hints_dict_minus.get(z,-1)!=(-1)): #if I have a start hint for z
					ext_evidence = True
					print_gm = False
				else: #no start hint so I need to search for different hints
					for y in pep_list_minus_start: #searching for a peptide hint
						if((z>=y>gm_start)and((z-y)%3==0)): #if I find an in-framt hint that starts between AUGUSTUS' start and GM'S
							ext_evidence = True
							print_gm = False
		else: #meaning: the searched stop site does not exist in GM
			ext_count = 0
			if(start_hints_dict_minus.get(z,-1)!=(-1)): #if I have a start hint for z
					ext_count = 1
					#ext_evidence = True
			pep_num = -1			
			for y in pep_list_minus_start:
				pep_num += 1
				if((x<=y<z)and((z-y)%3==0)and(pep_list_minus_stop[pep_num]>=x)): #if there is any evidence between start and stop site of the augustus gene
					ext_count +=1
			if(ext_count>=2):
				ext_evidence = True
					#ext_evidence = True
				#if(ext_evidence): #if I found extrinsic evidence for the cases needed
				#	output_handle.write(aug_gff_start + str(x) + "\t" + aug_list_stop[aug_position] +"\t" + "." + "\t" + aug_list_strand[aug_position] + "\t" + "." + "\t" + str(next_2+1) + "\n")
				#	next_2 += 1
				#ext_evidence = False
		if(print_gm): #print GM's gene because of missing evidence for AUGUSTUS' gene
			output_handle.write(gm_gff_start + aug_list_start[aug_position] +"\t" + gm_dictionairy_minus[aug_list_start[aug_position]] + "\t"+ "." +"\t" + aug_list_strand[aug_position] + "\t" + "." + "\t" + str(next_2+1) + "\n")
			next_2 += 1
			last_stop_minus = aug_list_start[aug_position]
		if(ext_evidence): #if I found extrinsic evidence for the cases needed
			output_handle.write(aug_gff_start + str(x) + "\t" + aug_list_stop[aug_position] +"\t" + "." + "\t" + aug_list_strand[aug_position] + "\t" + "." + "\t" + str(next_2+1) + "\n")
			next_2 += 1
			last_stop_minus = aug_list_start[aug_position]

#after iterating over all of AUGUSTUS' genes, check for missed GM genes:
if((next<(len(gm_list_start))) and (last_stop_plus == gm_list_stop[next] or last_stop_minus == gm_list_start[next])): 
		next += 1
while(next<len(gm_list_start)):
	if(gm_dictionairy_plus.get(gm_list_stop[next],-1)!=-1):
		output_handle.write(gm_gff_start+gm_list_start[next]+"\t"+gm_list_stop[next]+ "\t" + "." + "\t" +  "+" + "\t" + "." + "\t"+str(next_2+1)+"\n")
	else:
		output_handle.write(gm_gff_start + gm_list_start[next]+"\t" + gm_list_stop[next]+"\t"+ "." + "\t" +  "-" + "\t" + "." + "\t"+str(next_2+1)+"\n")			
	next += 1
	next_2 +=1

print(aug_position)
genemark_handle.close()
augustus_handle.close()
peptide_handle.close()
