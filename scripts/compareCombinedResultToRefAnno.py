#!/usr/bin/env python

import os
import argparse #needed for parsing the arguments
import re #needed for regular expressions

#Please not: This script says that a predicted gene is better than in the ref anno if there is extrinsic evidence in the part where those to genes differ.

parser = argparse.ArgumentParser(description='combines the output of GeneMarkS2 and AUGUSTUS depending on extrinsic evidence')
parser.add_argument('-p','--prediction',type=str,required=True,
	help='Joined gene prediction of GeneMark and AUGUSTUS (required)')
parser.add_argument('-r','--reference_annotation',type=str,required=True,
	help='Reference annotation (required)')
parser.add_argument('-ext','--extrinsic_evidence',type=str,required=True,
	help='Hints based on MS/MS experiments in gff format')
parser.add_argument('-o','--output',type=str,required=True,
	help='location of the outut file') 
parser.add_argument('-c','--correctly_predicted',type=str,required=False,
	help='Outputfile with all correctly predicted genes')
parser.add_argument('-u','--unfound_genes',type=str,required=False,
	help='Outputfile with all genes of the ref anno that were not found by the prediction')
args = parser.parse_args()

if(args.prediction):
	prediction_handle = open(args.prediction,"r")
if(args.reference_annotation):
	ref_handle = open(args.reference_annotation,"r")
if(args.extrinsic_evidence):
	peptide_handle = open(args.extrinsic_evidence,"r")
if(args.output):
	output_handle = open(args.output,"w")
if(args.correctly_predicted):
	correct_handle = open(args.correctly_predicted,"w")
if(args.unfound_genes):
	unfound_handle = open(args.unfound_genes,"w")

#THIS IS A SCRIPT THAT IS MEANT TO COMPARE THE RESULTS OF MY GENE PREDICTION TO THE REFERENCE ANNOTATION.
#THE FIRST THING TO DO IS FIND ALL THE GENES THAT GM AND AUG WERE ABLE TO PREDICT CORRECTLY (ACCORDING TO THE REFERENCE ANNOTATION).
#THE SECOND STEP WILL BE LOOK INTO THOSE GENES THAT WERE PREDICTED BY THOSE PROGRAMS IN A DIFFERENT WAY THAN IN THE REFERENCE ANNOTATION.
#"DIFFERENT WAY" MEANS EITHER THAT GM AND AUG PREDICTED A START SITE THAT LIES UPSTREAM OF THE ONE IN THE REFERENCE ANNOTATION (AND THE STOP SITE
#IS THE SAME) OR THAT THE TWO PROGRAMS PREDICT A COMPLETELY NEW GENE (NOT INCLUDED IN THE REFERENCE ANNOTATION). IN BOTH CASES THERE HAS TO BE 
#EXTRINSIC EVIDENCE IN ORDER TO JUDGE THE PREDICTION CORRECT (AND THEREBY THE REFERENCE ANNOTATION INCOMPLETE).
#ALL OTHER CASES OF DIFFERENT PREDICTION (e.g. A GENE START THAT LIES DOWNSTREAM) WILL BE JUDGED AS WRONG.

ref_dictionairy_plus = {} #the idea is to use the stop sites as keys and the start sites as values.
ref_dictionairy_minus = {} #again, stop sites as keys and start sites as values for the minus strand
ref_dictionairy_found_plus = {} #keys: stop sites, value: whole line if not found, True if found
ref_dictionairy_found_minus = {} #same
ref_anno_genes = 0
for x in ref_handle: #iterating over the reference annotation
	ref_match = re.search(r"^(.+\t\w+\t\w+\t)(\d+)\t(\d+)\t.\t([+-])\t",x) 
	if(ref_match):
		ref_anno_genes += 1
		if(ref_match.group(4)=="+"): #plus strand
			if(ref_dictionairy_plus.get(ref_match.group(3),-1)==(-1)):
				ref_dictionairy_plus[ref_match.group(3)]=ref_match.group(2)
				ref_dictionairy_found_plus[ref_match.group(3)]=x #saves the whole line as value
			else:
				if(int(ref_dictionairy_plus[ref_match.group(3)])<(-1)): #I alreade have more than one gene with this stop site
					ref_dictionairy_plus[ref_match.group(3)]=int(ref_dictionairy_plus[ref_match.group(3)])-1
					ref_dictionairy_plus[ref_match.group(3)+"_all"]=(ref_dictionairy_plus[ref_match.group(3)+"_all"]).append(ref_match.group(2))
					ref_dictionairy_found_plus[ref_match.group(3)]=ref_dictionairy_found_plus[ref_match.group(3)].append(x)
				else: #this is the second gene I found with this stop site
					ref_dictionairy_plus[ref_match.group(3)+"_all"]=[ref_dictionairy_plus[ref_match.group(3)],ref_match.group(2)]
					ref_dictionairy_plus[ref_match.group(3)]=-2
					ref_dictionairy_found_plus[ref_match.group(3)]=[ref_dictionairy_found_plus[ref_match.group(3)],x]
		else: #minus strand
			if(ref_dictionairy_minus.get(ref_match.group(2),-1)==(-1)):
				ref_dictionairy_minus[ref_match.group(2)]=ref_match.group(3)
				ref_dictionairy_found_minus[ref_match.group(2)]=x #saves the whole line as value
			else:
				if(int(ref_dictionairy_minus[ref_match.group(2)])<(-1)): #I alreade have more than one gene with this stop site
					ref_dictionairy_minus[ref_match.group(2)]=int(ref_dictionairy_minus[ref_match.group(2)])-1
					ref_dictionairy_minus[ref_match.group(2)+"_all"]=(ref_dictionairy_minus[ref_match.group(2)+"_all"]).append(ref_match.group(3))
					ref_dictionairy_found_minus[ref_match.group(2)]=ref_dictionairy_found_minus[ref_match.group(2)].append(x)
				else: #this is the second gene I found with this stop site
					ref_dictionairy_minus[ref_match.group(2)+"_all"]=[ref_dictionairy_minus[ref_match.group(2)],ref_match.group(3)]
					ref_dictionairy_minus[ref_match.group(2)]=-2
					ref_dictionairy_found_minus[ref_match.group(2)]=[ref_dictionairy_found_minus[ref_match.group(2)],x]
pep_list_plus = [] #this will be a LIST of all start sites of peptide evidence on the plus strand
pep_list_minus = [] #same for the minus strand
for x in peptide_handle: #iterating over the peptides hints
	pep_match = re.search(r"^.+\t\w+\t\w+\t(\d+)\t(\d+)\t.\t([+-])\t",x)
	if(pep_match):
		if(pep_match.group(3)=="+"): #plus strand
			pep_list_plus.append(int(pep_match.group(1))) 
		else: #minus strand
			pep_list_minus.append(int(pep_match.group(2)))

#At first I will implement counting the different types of genes: correct, more correct than ref (meaning there is good evidence for the difference) , in part and wrong
correct = 0 #int for counting the number of completely correct found genes (according to the ref anno)
better = 0 #int for counting the number of genes that were in part or completely different than in the ref anno but have ext evidence
in_part = 0#int for counting the genes that were predicted partly correct
wrong = 0 #int for counting the genes that were predicted in a way not existing in the ref anno at all and have no ext evidence whatsoever
pred_genes = 0 #number of all predicted genes
correct_plus_evidence = 0 #number of genes that were correct and had extrinsic evidence 
same_stop = 0 #number of genes that have the same stop as some reference annotation's gene
for x in prediction_handle: #iterating over the gene prediction file
	pred_match = re.search(r"^(.+\t\w+\t\w+\t)(\d+)\t(\d+)\t.\t([+-])\t",x)
	if(pred_match):
		pred_genes += 1 # counts the number of predicted genes
		ext_evidence = False #boolean variable for storing the information whether there is extrinsic evidence for this gene prediction
		if(pred_match.group(4)=="+"): #plus strand
			ref_start = int(ref_dictionairy_plus.get(pred_match.group(3),-1)) #returns -1 if the stop site does not exist in the ref anno
			if(ref_start==int(pred_match.group(2))): #same start site
				correct += 1
				if(ref_dictionairy_found_plus[pred_match.group(3)]!=True):
					same_stop += 1
				ref_dictionairy_found_plus[pred_match.group(3)] =True #sets value from whole line to True
				ext_evidence_correct = False
				for y in pep_list_plus:
						if((int(pred_match.group(2))<=y<int(pred_match.group(3)))and((y-int(pred_match.group(2)))%3==0)):
							ext_evidence_correct = True
				if(ext_evidence_correct): #I found evidence for the predicted difference
					correct_plus_evidence += 1
				#del ref_dictionairy_plus[pred_match.group(3)]
			
			elif(ref_start < (-1)):
				same_stop += 1
				ref_dictionairy_found_plus[pred_match.group(3)] =True #sets value from whole line to True
				for i in ref_dictionairy_plus[pred_match.group(3)+"_all"]:
					if(i==int(pred_match.group(2))):
						correct += 1
						ext_evidence_correct = False
						for y in pep_list_plus:
							if((int(pred_match.group(2))<=y<int(pred_match.group(3)))and((y-int(pred_match.group(2)))%3==0)):
								ext_evidence_correct = True
					if(ext_evidence_correct): #I found evidence for the predicted difference
						correct_plus_evidence += 1
			elif(ref_start > -1):#meaning the stop site does exist in the ref anno but apparently the start site is different
				#now I need to test whether I have extrinsic evidence for this difference
				if(args.correctly_predicted): #if parameter -c was chosen
					correct_handle.write(x) #write all genes into it that have the correct stop site but not the correct start site
				if(ref_dictionairy_found_plus[pred_match.group(3)]!=True):
					same_stop += 1
				ref_dictionairy_found_plus[pred_match.group(3)] =True
				if(ref_start < int(pred_match.group(2))): #meaning: the reference start site is upstream of the one in the prediction
					in_part += 1
				else: #the reference start site is downstream of the one in the prediction
					for y in pep_list_plus: #search for peptide evidence
						if((int(pred_match.group(2))<=y<ref_start)and((y-int(pred_match.group(2)))%3==0)):
							ext_evidence = True
					if(ext_evidence): #I found evidence for the predicted difference
						better += 1
						output_handle.write(x)
					else:
						wrong +=1
			else: #meaning this is a gene that does not exist at all in the ref anno
				for y in pep_list_plus: #search for evidence
					if((int(pred_match.group(2))<=y<int(pred_match.group(3)))and((y-int(pred_match.group(2)))%3==0)): #is this enough? no
					#only testing whether there is some peptide STARTING in the range between start and stop site
					#seems like light-weighted evidence to me
						ext_evidence = True
				if(ext_evidence): #I found evidence for the predicted difference
					better += 1
				else:
					wrong +=1
					
		else: #minus strand
			ref_start = int(ref_dictionairy_minus.get(pred_match.group(2),-1)) #returns -1 if the stop site does not exist in gm
			if(ref_start==int(pred_match.group(3))):
				correct +=1
				if(ref_dictionairy_found_minus[pred_match.group(2)]!=True):
					same_stop += 1				
				ref_dictionairy_found_minus[pred_match.group(2)] =True
				ext_evidence_correct = False
				for y in pep_list_plus:
						if((int(pred_match.group(2))<=y<int(pred_match.group(3)))and((int(pred_match.group(3))-y)%3==0)):
							ext_evidence_correct = True
				if(ext_evidence_correct):
					correct_plus_evidence += 1
				#del ref_dictionairy_minus[pred_match.group(3)]
			elif(ref_start < (-1)):
				same_stop += 1
				ref_dictionairy_found_minus[pred_match.group(2)] =True #sets value from whole line to True
				for i in ref_dictionairy_minus[pred_match.group(2)+"_all"]:
					if(i==int(pred_match.group(3))):
						correct += 1
						ext_evidence_correct = False
						for y in pep_list_plus:
							if((int(pred_match.group(2))<=y<int(pred_match.group(3)))and((int(pred_match.group(3))-y)%3==0)):
								ext_evidence_correct = True
					if(ext_evidence_correct): #I found evidence for the predicted difference
						correct_plus_evidence += 1
			elif(ref_start > -1):#meaning the stop site does exist in the ref anno but apparently the start site is different
				#now I need to test whether I have extrinsic evidence for this difference
				if(args.correctly_predicted):
					correct_handle.write(x)
				if(ref_dictionairy_found_minus[pred_match.group(2)]!=True):
					same_stop += 1				
				ref_dictionairy_found_minus[pred_match.group(2)] =True
				if(ref_start > int(pred_match.group(3))): #meaning: the reference start site is upstream of the one in the prediction
					in_part += 1
				else:
					for y in pep_list_minus:
						if((int(pred_match.group(3))>=y>ref_start)and((int(pred_match.group(3))-y)%3==0)):
							ext_evidence = True
					if(ext_evidence): #I found evidence for the predicted difference
						better += 1
						output_handle.write(x)
					else:
						wrong +=1
			else: #meaning this is a gene that does not exist at all in the ref anno
				for y in pep_list_minus:
					if((int(pred_match.group(2))<=y<int(pred_match.group(3)))and((int(pred_match.group(3))-y)%3==0)):
						ext_evidence = True
				if(ext_evidence): #I found evidence for the predicted difference
					better += 1
				else:
					wrong +=1

#all genes of the reference annotation that were not found by the prediction will be written into the unfound file: 
unfound_handle.write("These are all the genes of the reference annotation that have no gene in the prediction : " + args.prediction + " that shares the same stop codon. \n")
for ref_gene in ref_dictionairy_found_plus:
	if(ref_dictionairy_found_plus[ref_gene]!=True):
		if(int(ref_dictionairy_plus[ref_gene])<(-1)):
			for i in ref_dictionairy_found_plus[ref_gene]:
				if(i!=True):
					unfound_handle_write(i)
		else:
			unfound_handle.write(ref_dictionairy_found_plus[ref_gene])
for ref_gene in ref_dictionairy_found_minus:
	if(ref_dictionairy_found_minus[ref_gene]!=True):
		if(int(ref_dictionairy_minus[ref_gene])<(-1)):
			for i in ref_dictionairy_found_minus[ref_gene]:
				if(i!=True):
					unfound_handle_write(i)
		else:
			unfound_handle.write(ref_dictionairy_found_minus[ref_gene])

#all other infos will be written into the output file:
output_handle.write("total number of genes in the reference annotation: " + str(ref_anno_genes) + "\n")
output_handle.write("total number of genes predicted by GeneMark and AUGUSTUS: " + str(pred_genes) + "\n")
#output_handle.write("Number of correctly predicted genes: " + str(correct) + "\n")
#output_handle.write("Number of genes that were predicted more accurately than in the ref anno: " + str(better) + "\n")
#output_handle.write("Number of genes that were predicted only partly correct (meaning different start, same stop site): " +str(in_part) + "\n")
#output_handle.write("Number of predicted genes that were wrong: " +str(wrong) + "\n")
output_handle.write("\n")
#output_handle.write("This leads to the following numbers: \n")
output_handle.write("True Positives(start&stop): " + str(correct) + "\n")
output_handle.write("False Positives(start&stop): " + str(better)+"+"+str(in_part)+"+"+str(wrong)+"=\t" +str(better+in_part+wrong) + "\n")
output_handle.write("False Negatives(start&stop): " + str(int(str(len(ref_dictionairy_plus)+len(ref_dictionairy_minus)) + "\n")-correct) +"\n")
output_handle.write("Sensitivitiy = TruePositives/#RefAnno:\t" + str(correct) + "/" + str(ref_anno_genes) + " =\t" + str(int(correct)/int(ref_anno_genes)) + "\n")
output_handle.write("Specificity = #TruePositives/#Predictions :\t" + str(correct) + "/" + str(pred_genes) + "=\t" + str(correct/pred_genes) + "\n")
output_handle.write("Harmonic mean: " + str(2/(1/(int(correct)/int(ref_anno_genes))+1/(correct/pred_genes))) + "\n")
output_handle.write("\n")
#output_handle.write("Number of correctly predicted genes that additionally have extrinsic evidence: " + str(correct_plus_evidence) + "\n")
output_handle.write("True Positives(stop): " + str(same_stop) + "\n")
output_handle.write("Sensitivity(stop): " + str(same_stop/(int(ref_anno_genes))) + "\n")
output_handle.write("Specificity(stop): " + str(same_stop/pred_genes) + "\n")
output_handle.write("Harmonic mean: " + str(2/(1/(int(same_stop)/int(ref_anno_genes))+1/(same_stop/pred_genes))) + "\n")
prediction_handle.close()
ref_handle.close()
peptide_handle.close()
output_handle.close()
correct_handle.close()
unfound_handle.close()
