#!/usr/bin/env python3


#Author: Leonie Johanna Lorenz
#Last modified: 18th September 2019

# PeptImprove is a pipeline for prokaryotic genome annotation.
# Given a genome FASTA file, the species name and extrinsic evidence
# in form of mass-spectrometric data in raw file format,
# PeptImprove will fully automated make a gene prediction using GeneMarkS-2 and AUGUSTUS.
# If desired, PeptImprove can instead be used to improve a given reference annotation.


import re #needed for regular expressions
import sys #needed for arguments given to the script
import os #needed to read all the protein files given in the directory
import argparse
import subprocess
import multiprocessing
from functools import partial
import shutil

### ARGUMENTS ###
parser = argparse.ArgumentParser(
    description='A pipeline for prokaryotic gene prediction using IdentiPy, GeneMarkS-2 and AUGUSTUS.')
parser.add_argument('-g', '--genome_file', type=str,required=True,
                    help="Genome input file in FASTA format")
parser.add_argument('-e','--extrinsic_evidence',type=str, nargs='+',required=True,
	help="The location of the directory containing all MS data. You may name multiple directories.")
parser.add_argument('-n','--name',type=str,required=True,
	help="Please give the name of your species (avoid white spaces as it will also function as the name of a directory).")
parser.add_argument('--cores',type=int,required=False,
	help="If you have information about the number of available cores you can enter it here. Default will be 4.")
parser.add_argument('--re_run', action='store_true',
	help="Option to leave out GeneMark prediction and AUGUSTUS training. Only makes sense when the pipeline was used before on the same species but with different extrinsic evidence.")
parser.add_argument('-r','--reference_annotation',type=str,required=False,
	help="Reference annotation is gff format. Can be used for comparing the results of for an different usage of the pipeline.")
parser.add_argument('--use_ref',action='store_true',
	help="Option to use reference annotation for training AUGUSTUS. GeneMarkS2 won't be run in this version.")
parser.add_argument('--no_optimization',required=False, action='store_true',
	help="Option to omit optimize_augustus. Only advisable if there is not enough time or computational capacity.")
parser.add_argument('--output_directory',required=False,type=str,
	help="Option to name the output directory and to define where it should be. If this argument is not given the default directory 'output' will be created where the pipeline was called.")
parser.add_argument('--augustus_config_path',required=False,type=str,
	help="Option give the path to the AUGUSTUS config directory.")
parser.add_argument('--path_to_rawfileparser',required=False,type=str,
	help="Option give the path to the directory of the ThermoRawFileParser.")
parser.add_argument('--delete_intermediate',required=False, action='store_true',
	help="Option to delete all results except the final one.")
parser.add_argument('--archea',required=False, action='store_true',
	help="Use option if species is of type Archea.")
parser.add_argument('--translation_table',type=int,required=False,
	help="If the translation table of the species' genome is not 11, please enter the number her. Default will be 11.")
#The following arguments are for experts:
parser.add_argument('--start',type=str,required=False,
	help="You may give an integer (incl. 0) that indicates how many steps should be left out at the beginning of the pipeline.")
parser.add_argument('--stop',type=str,required=False,
	help="You may give an integer <=9 that indicates how many steps should be left out at the end of the pipeline.")
parser.add_argument('-c','--compare_to_ref',action='store_true',
	help="Option to compare results to a reference annotation. Only possible when the reference annotation is given via -r.")
parser.add_argument('-f','--filter_ip',type=float,required=False,
	help="Value at which the hints by IdentiPy should be filtered. 0.01 stands for 1 percent. Default is 0.001.")
parser.add_argument('--extrinsic_cfg',required=False,type=str,
	help="Option to use a specific Cfg file that is derived from the original one. Elsewise the pipeline's default Cfg file is used.")
args = parser.parse_args()

if(args.genome_file):
	#genome_handle = open(args.genome_file,"r")
	genome_file_name = args.genome_file
if(args.extrinsic_evidence):
	#ext_handle = open(args.extrinsic_evidence,"r")
	ext_file_directory = args.extrinsic_evidence
start_at = 0
if(args.start):
	start_at = int(args.start)
stop_at = 20
if(args.stop):
	stop_at = int(args.stop)
if(args.name):
	species = args.name
cores_available = 4
if(args.cores):
	cores_available = int(args.cores)
if(args.reference_annotation):
	reference = os.path.abspath(args.reference_annotation)
	print(reference)
else:
	if(args.compare_to_ref or args.use_ref):
		print("Warning: You chose an option that needs the reference annotation as input via -r!\n")
hints_filter = 0.001
if(args.filter_ip):
	hints_filter = args.filter_ip
ext_cfg_file = "extrinsic_prok.MP.cfg"
if(args.extrinsic_cfg):
	ext_cfg_file=args.extrinsic_cfg
optimize = True
if(args.no_optimization):
	optimize = False
path_to_output = "./output"
if(args.output_directory):
	path_to_output = args.output_directory
path_to_output = os.path.abspath(path_to_output)
delete_results = False
if(args.delete_intermediate):
	delete_results = True
transl_table = 11
if(args.translation_table):
	transl_table = int(args.translation_table)

### ARGUMENTS ###

# To make sure no more cores are accessed than defined by the user:
command1 = "taskset -p -c 1"
if(cores_available > 1):
	for x in range(2,cores_available+1):
		command1 = command1 + "," + str(x)
command1 = command1 + " %d"
os.system(command1 % os.getpid())


# Store path to PeptImprove
path_to_pipeline, filename1 = os.path.split(str(os.path.abspath(__file__)))
#path_to_call_location = os.getcwd()
#if(not re.search(r"^/",path_to_output)):
#	path_to_output = path_to_call_location + "/" + path_to_output
### BEGIN FUNCTIONS ###

# Find augustus_config_path
def find_config(general_augustus):
	augustus_config=""
	if (args.augustus_config_path):
		augustus_config = args.augustus_config_path
	else:
		if os.environ.get('AUGUSTUS_CONFIG_PATH') is not None:
			test_augustus_config = os.environ.get('AUGUSTUS_CONFIG_PATH')
			if (os.path.exists(test_augustus_config)):
				augustus_config = test_augustus_config
	if augustus_config == "":
		test_augustus_config = general_augustus + "/config/"
		if (os.path.exists(test_augustus_config)):
			augustus_config = test_augustus_config
		else:
			print("Unable to find the AUGUSTUS_CONFIG_PATH. Please check whether you installed AUGUSTUS properly and if so use --augustus_config_path to enter the path manually.")
			exit(1)
	return(augustus_config)

# Find path to ThermoRawFileParser
def locate_rawfileparser():
	path = ""
	if(args.path_to_rawfileparser):
		path = args.path_to_rawfileparser
	else:
		if(os.path.exists(path_to_pipeline + "/ThermoRawFileParser")):
			path = path_to_pipeline + "/ThermoRawFileParser"
		else:
			print("Uanble to find the location of the ThermoRawFileParser. Either store it as a subdirectory of your PeptImprove's directory or enter its path via --path_to_rawfileparser.")
			exit(1)
	return(path)
# Find path to ThermoRawFileParser
def locate_getOrfsOrCDSs():
	path = ""
	if(not(shutil.which("get_orfs_or_cdss.py")== None)):
		path = shutil.which("get_orfs_or_cdss.py")
	else:
		if(os.path.exists(path_to_pipeline + "/get_orfs_or_cdss.py")):
			path = path_to_pipeline
		else:
			print("Unable to find the location of the script get_orfs_or_cdss.py. Please either add it to your $PATH and make it executable or store it in PeptImprove's directory.")
			exit(1)
	return(path)

# Find path to ThermoRawFileParser
def locate_augustus_scripts(augustus_general_directory):
	path = ""
	if(not(shutil.which("optimize_augustus.pl")== None)):
		path = shutil.which("optimize_augustus.pl")
	else:
		if(os.path.exists(augustus_general_directory + "/scripts")):
			path = augustus_general_directory + "/scripts"
		else:
			print("Unable to find the location of Augustus/scripts. Please either add it to your $PATH or it next to the Augustus/bin.")
			exit(1)
	return(path)

# Function to locate all necessary binaries
def locate_binaries(tool):
	if(shutil.which(tool) == None):
		print("I could not locate " + tool + ". Please make it available in you $PATH.")
		exit(1)
	return(shutil.which(tool))

# Function to create the 6-frame translation. It probably would make sense to return the name of the file created here. (local variable)
def create6frame(genome, sixframe):
	print("Creating the six-frame translation\n")
	try:
		six_frame = subprocess.run("python3 " + path_to_getOrfsOrCDSs + " -i" + genome + " -t CDS --table " + str(transl_table) + " --op " + sixframe, shell=True) #subprocess that is created to run get_orfs_or_cdss.py. Output: sixframe_name. Later on: create a handle corresponding to that file
	except IOError:
		print("The six-frame translation could not be created.\n")
#function for parsing the raw data to mgf format using HelpForRawFileParser.py and in it ThermoRawFilerParser
def parseRawData(ext_file_directory, mgf,path_to_parser):
	print("Parsing RAW data to MGF format\n")
	try:
		for x in ext_file_directory:
			raw_to_mgf = subprocess.run("python3 " + path_to_pipeline + "/HelpForRawFileParserMgf.py " + x + " " + path_to_parser + " " + mgf, shell=True)
		rm_txt_files = subprocess.run("rm "+ mgf + "/*.txt", shell=True)
	except IOError:
		print("The raw files could not be parsed.\n")

#function for using IdentiPy (with HelpForIdentiPy.py)
def run_ip(path_to_ip_output, path_to_6_frame, mgf_direc, mgf_file):			
	return subprocess.run('identipy ' + mgf_direc + "/" + mgf_file + ' -out ' + path_to_ip_output + ' -db ' + path_to_6_frame, shell=True)
def runningIdentiPy(sixframe, mgf, ip_direc, cores):
	print("IdentiPy is identifying peptides\n")
	try:
		if(os.path.exists(ip_direc)):
			rm_directory = subprocess.run("rm -r " + ip_direc, shell = True)
		create_directory=subprocess.run('mkdir ' + ip_direc,shell=True) #creating the sub.process that runs the bash command
		
		mgf_list = os.listdir(mgf)
		if __name__ == "__main__":
			with multiprocessing.Pool(processes = cores) as pool:
				function_to_run_ip = partial(run_ip, ip_direc, sixframe, mgf)
				results = pool.map(function_to_run_ip, mgf_list)
	except IOError:
		print("There is a problem with running IdentiPy. Maybe it is not in the path variable? \n")
		exit(1)

#Function for parsing IP's output to GFF format
def parseIdentiPyOutput(ip_direc, joinedfile, outfile, sixframe,filter_hints):
	print("Parsing IdentiPy's output\n")
	try:		
		createJoinedFile = subprocess.run("cat " + ip_direc + "/*.pep.xml > " + path_to_output + "/" + joinedfile, shell = True)
		moveJoinedFile = subprocess.run("mv " + path_to_output + "/" + joinedfile + " " + ip_direc + "/",shell=True)
		parseFile = subprocess.run("python3 " + path_to_pipeline + "/identipyToGffAndFilterAndstartHints.py " + sixframe + " " + ip_direc + "/" + joinedfile + " " + outfile + " " + str(filter_hints), shell=True)
		sortHints = subprocess.run("sort " + outfile + "> Hints_sorted.gff",shell=True) #sort hints
		collapseMult = subprocess.run("python3 " + path_to_pipeline + "/collapseRedundantHints.py Hints_sorted.gff Hints_collapsed.gff",shell=True)
		mvNonred = subprocess.run("mv Hints_collapsed.gff "+ outfile, shell=True)
	except IOError:
		print("I cannot join and parse IdentiPy's output files.\n")

#Function for running GeneMarkS2
def runningGMS2(genome, gm_output):
	print("Predicting genes with GeneMarkS-2\n")
	try:
		if(not args.archea):
				runGM = subprocess.run("gms2.pl --genome-type bacteria --fnn " + gm_output +" --seq " + genome, shell=True)
		else:
				runGM = subprocess.run("gms2.pl --genome-type archea --fnn " + gm_output +" --seq " + genome, shell=True)
	except:
		print("Running GeneMarkS2 does not work properly.\n")

def parsingGMS2Output(gm_output,gm_gff):
	print("Parsing GeneMarkS2's output\n")
	try:
		parseGM = subprocess.run("python3 " + path_to_pipeline + "/GMOutputToGFF.py " + gm_output + " " + gm_gff, shell=True)
	except IOError:
		print("Parsing GeneMarkS2's output does not work properly.\n")

def trainingAUGUSTUS(genome, training_all,cores, species):
	print("Training AUGUSTUS\n")
	try:
		flankingRegion = subprocess.check_output("perl " + augustus_path + "/computeFlankingRegion.pl " + training_all, shell=True) #computes possible length of flanking regions
		minflanking = re.search(r"The\sflanking\_DNA\svalue\sis:\s(\d+)\s",str(flankingRegion))#in minflanking the minimal possible flanking region is stored correctly!
		createGenbankFile = subprocess.run("perl " + augustus_path + "/gff2gbSmallDNA.pl --overlap " + training_all + " " + genome + " " + str(int(minflanking.group(1))) + " " + path_to_output + "/training_all.gb", shell=True) #creates genbank flat file containg GM's genes plus flanking regions
		createSoftLink = subprocess.run("ln -s " + training_all + " " + path_to_output +"/training_all.f.gtf",shell=True) #creates soft link to the gff file
		convertToaa = subprocess.run("python3 " + augustus_path + "/getAnnoFastaFromJoingenes.py -g " + genome + " -t " + str(transl_table) + " -f " + path_to_output+ "/training_all.f.gtf -o " + path_to_output + "/protSeq.fa -s True",shell=True)
		mvaa = subprocess.run("mv " + path_to_output + "/protSeq.fa.aa " + path_to_output + "/prot.aa",shell=True)
		blastaa = subprocess.run("perl " + augustus_path + "/aa2nonred.pl " + path_to_output + "/prot.aa " + path_to_output+ "/prot.nr.aa --cores=" + str(cores), shell=True) #makes a blast run for filtering for genes that are not identical in more than 80% on amino acid level
		getNonred = subprocess.run("grep '>' " + path_to_output + "/prot.nr.aa | perl -pe 's/>//' > " + path_to_output + "/nonred.lst", shell=True) #filters for non-redundant genes
		partOne = "cat " + path_to_output + "/training_all.gb | perl -ne ' "
		partTwo = 'if($_=~ m/LOCUS\s+(\S+)\s/){$txLocus = $1;} elsif ($_ =~ m/\/gene=\\"(\S+)\\"/){$txInGb3{$1}=$txLocus} if(eof()) {foreach (keys %txInGb3) {print "$_\\t$txInGb3{$_}\\n";}}'
		partThree= "' > " + path_to_output + "/loci.lst"
		createLociList = subprocess.run(partOne + partTwo + partThree, shell=True)	
		makeLociListnonred = subprocess.run("grep -wf " + path_to_output + "/nonred.lst " + path_to_output + "/loci.lst | cut -f2 > " + path_to_output + "/nonred.loci.lst", shell=True) #loci list is now nonred
		filterTrainingGenes = subprocess.run("perl " + augustus_path + "/filterGenesIn.pl " + path_to_output +"/nonred.loci.lst " + path_to_output + "/training_all.gb > " + path_to_output + "/GbFromGM.f.gb",shell=True)
		mvFilteredGb = subprocess.run("mv " + path_to_output + "/GbFromGM.f.gb " + path_to_output + "/training_all.gb", shell=True) #saves the filtered version as the normal one
		createSpeciesDirec = subprocess.run("perl " + augustus_path + "/new_species.pl --AUGUSTUS_CONFIG_PATH=" + augustus_config_path + " --prokaryotic --species=" + species, shell=True) #creates a directory into which AUGUSTUS will store all its species specific parameter files
		createTestSet = subprocess.run("perl " + augustus_path + "/randomSplit.pl " + path_to_output + "/training_all.gb 200", shell=True) #creates a test and a traing set
		nameTestSet = subprocess.run("mv " + path_to_output + "/training_all.gb.test " + path_to_output + "/test.gb",shell=True)
		nameTrainingSet = subprocess.run("mv " + path_to_output + "/training_all.gb.train " + path_to_output + "/train.gb",shell=True)
		etrain1 = subprocess.run("etraining --species=" + species + " " + path_to_output + "/train.gb > " + path_to_output + "/etrain.out",shell=True) #run first etraining
		tagFreq = 0.34
		taaFreq = 0.33
		tgaFreq = 0.33
		etrainHandle = open(path_to_output + "/etrain.out","r")
		for x in etrainHandle:
			matchTAG = re.search(r"^tag.\s+\d+\s+\((.+)\)",x) #matching the frequency of tag as a stop codon
			matchTAA = re.search(r"^taa.\s+\d+\s+\((.+)\)",x)#matching the frequency of taa as a stop codon
			matchTGA = re.search(r"^tga.\s+\d+\s+\((.+)\)",x)#matching the frequency of tga as a stop codon
			if(matchTAG):
				tagFreq = float(matchTAG.group(1))
			if(matchTAA):
				taaFreq = float(matchTAA.group(1))
			if(matchTGA):
				tgaFreq = float(matchTGA.group(1))
		etrainHandle.close()
		
		paramFileName = species + "_parameters.cfg"
		helpParamHandle = open(path_to_output + "/" + paramFileName,"w")
		paramFileHandle = open(augustus_config_path + "/species/" + species + "/" + paramFileName, "r")
		
		for x in paramFileHandle:
			tagMatch = re.search(r"(\/Constant\/amberprob\s+)0\.\d+(\s.+)\n",x)
			taaMatch = re.search(r"(\/Constant\/ochreprob\s+)0\.\d+(\s.+)\n",x)
			tgaMatch = re.search(r"(\/Constant\/opalprob\s+)0\.\d+(\s.+)\n",x)
			if(tagMatch):
				helpParamHandle.write(tagMatch.group(1) + str(tagFreq) + tagMatch.group(2) + "\n")
			elif(taaMatch):
				helpParamHandle.write(taaMatch.group(1) + str(taaFreq) + taaMatch.group(2) + "\n")
			elif(tgaMatch):
				helpParamHandle.write(tgaMatch.group(1) + str(tgaFreq) + tgaMatch.group(2) + "\n")
			else:
				helpParamHandle.write(x)
		helpParamHandle.close()
		paramFileHandle.close()
		mvParamFile = subprocess.run("mv " + path_to_output + "/" +paramFileName + " " + augustus_config_path + "/species/" + species + "/",shell=True)
	except IOError:
		print("Training AUGUSTUS did not work properly.\n")

def optimizingAUGUSTUS(train,species,cores):
	print("Further improving AUGUSTUS' model (this may take a while)\n")
	if(int(cores)<=8):
		optAUG = subprocess.run("perl " + augustus_path + "/optimize_augustus.pl --AUGUSTUS_CONFIG_PATH=" + augustus_config_path + " --species=" + species + " --kfold=" + str(8) + " --cpus=" + str(cores) + " " + train + " > optimize.out",shell=True)
	else: #many cores 
		lineCount = 0
		countTrainGenesHandle = open(train,"r")
		for x in countTrainGenesHandle:
			lineCount += 1
		countTrainGenesHandle.close()
		if((int(lineCount)/int(cores))>1): #if there is at least one gene in each training file
			optAUG = subprocess.run("perl " + augustus_path + "/optimize_augustus.pl --AUGUSTUS_CONFIG_PATH=" + augustus_config_path + " --species=" + species + " --kfold=" + str(int(cores)) + " --cpus=" + str(cores) + " " + train + " > optimize.out",shell=True)
		elif((2*int(lineCount)/(int(cores)))>1): 
			optAUG = subprocess.run("perl " + augustus_path + "/optimize_augustus.pl --AUGUSTUS_CONFIG_PATH=" + augustus_config_path + " --species=" + species + " --kfold=" + str(int(cores)/2) + " --cpus=" + str(cores) + " " + train + " > optimize.out",shell=True)
		else:
			optAUG = subprocess.run("perl " + augustus_path + "/optimize_augustus.pl --AUGUSTUS_CONFIG_PATH=" + augustus_config_path + " --species=" + species + " --kfold=" + str(8) + " --cpus=" + str(cores) + " " + train + " > optimize.out",shell=True)
	finalEtrain = subprocess.run("etraining --species=" + species + " " + train, shell=True)

def predictWithAUGUSTUS(genome, species, hints, outfile,ext_cfg):
	print("Predicting genes with AUGUSTUS\n")
	runAug = subprocess.run("augustus --species=" + species + " --extrinsicCfgFile=" + path_to_pipeline + "/"+ ext_cfg + " --hintsfile=" + hints + " " + genome + " --gff3=on --outfile=" + path_to_output + "/augustus.gff3 --errfile=" + path_to_output + "/augustus.err", shell=True)
	#print("augustus --species=" + species + " --extrinsicCfgFile=" + "/home/leonie/Augustus/config/species/" + species + "/" + species + "_parameters.cfg" + " --hintsfile=" +hints + " --gff3=on --outfile=augustus.gff3 --errfile=augustus.err\n")
	parseAug = subprocess.run("python3 " + path_to_pipeline + "/AugOutputToGff.py " + path_to_output+ "/augustus.gff3 " + outfile, shell=True)

def combineGMandAUGUSTUS(genemark, augustus, hints, outfile):
	print("Combining the GeneMarkS2 and the AUGUSTUS predictions\n")
	combine = subprocess.run("python3 " + path_to_pipeline + "/combineGMandAUGUSTUS.py -gm " + genemark + " -aug " + augustus + " -ext " + hints + " -o " + outfile, shell=True)

def compareToRef(prediction, reference, hints, outfile):
	print("Comparing results to the reference annotation\n")
	compare = subprocess.run("python3 " + path_to_pipeline + "/compareCombinedResultToRefAnno.py -p " + prediction + " -r " + reference + " -ext " + hints + " -o " + outfile + " -c " + path_to_output + "/correct_stop.gff -u " + path_to_output + "/unfound_genes.gff", shell=True)

### END FUNCTIONS ###
	
### CHECKING PROGRAMS AND DIRECTORIES ###
locate_binaries("python3")
locate_binaries("perl")
locate_binaries("mono")
if(not(args.use_ref)):
	locate_binaries("gms2.pl")
locate_binaries("blastp")
locate_binaries("identipy")
path_to_parser = locate_rawfileparser()
path_to_getOrfsOrCDSs = locate_getOrfsOrCDSs()
augustus_bin = locate_binaries("augustus")
augustus_path_general, binariesDirectory = augustus_bin.split("/bin")
augustus_path = locate_augustus_scripts(augustus_path_general)
augustus_config_path = find_config(augustus_path_general)
### Checking the output directory
if(start_at == 0):
	##Create output directory
	if (os.path.isdir(path_to_output)):
		if (os.listdir(path_to_output)):
			print("The directory you entered already exists and is not empty. Please name another directory or empty this one.\n")
			exit(1)
	path = path_to_output.rsplit('/',1)
	if(not os.access(path[0], os.W_OK)):
		print("Do not have the permission to create output directory. Make path writeable or choose different path.\n")
		exit(1)
	else:
		create_directory = subprocess.run(["mkdir",path_to_output], shell=False)
else: #check whether the given directory is writeable
	if (os.path.isdir(path_to_output)):
		if(not(os.access(path_to_output,os.W_OK))):
			print("Do not have the permission to write into the output directory. Please name a different one.")
	else:
		print("You started PeptImprove at step " + str(start_at) + ". If doing so, the output files of all former steps must already be given in the output directory but your output directory does not exist as of yet. Please either start PeptImprove from the beginning or name a different output directory.")


### FUNCTIONS CALLS ###
#Calling the function to create the 6 frame translation
sixframe_name = path_to_output + "/6_frame_translation_CDS_prot.fa" #this is the name of the file with the six-frame translation
if(start_at == 0 and stop_at >0):
	create6frame(genome_file_name, sixframe_name)
#Calling the function to parse the raw data
mgf_directory = path_to_output + "/mgf/"
if(start_at <=1 and stop_at > 1):
	if(not(os.path.exists(mgf_directory))):
		createMGFDirectory = subprocess.run("mkdir " + mgf_directory, shell=True)
	parseRawData(ext_file_directory, mgf_directory, path_to_parser)
#Calling the function to run IdentiPy
ip_output_directory = path_to_output + "/ip_output"
if(start_at <= 2 and stop_at>2):
	runningIdentiPy(sixframe_name, mgf_directory, ip_output_directory, cores_available)
#Calling the function to parse IdentiPy's output
joined_ip_xml_hints = "allIPpeptides.pep.xml"
ip_hints = path_to_output + "/HintsByIP.gff"
if(start_at <=3 and stop_at>3):
	parseIdentiPyOutput(ip_output_directory,joined_ip_xml_hints, ip_hints, sixframe_name,hints_filter)
#Defining the first gene set (either GeneMarkS-2 prediction or reference annotation)
training_genes = path_to_output + "/GffFromGMS2.gff"
if(not (args.use_ref)):
	gm_output_list = path_to_output + "/nucleotides.lst"
	if(start_at <= 4 and not(args.re_run) and stop_at>4): #if pipeline not only runs because of new evidence
		runningGMS2(genome_file_name, gm_output_list)
	gm_gff_file = path_to_output + "/GffFromGMS2.gff"
	if(start_at <= 5 and not(args.re_run) and stop_at>5): #if pipeline not only runs because of new evidence
		parsingGMS2Output(gm_output_list, gm_gff_file)
	training_genes = gm_gff_file
else: #user defined that the reference should be used instead of GM
	training_genes = reference
#Calling function to train AUGUSTUS
train_file = path_to_output + "/train.gb"
test_file = path_to_output + "/test.gb"
if(start_at <= 6 and not(args.re_run) and stop_at>6):#if pipeline not only run because of new evidence
	trainingAUGUSTUS(genome_file_name,training_genes, cores_available, species)
#Calling function to optimize AUGUSTUS
if(start_at <=7 and not(args.re_run) and stop_at>7 and optimize):#if pipeline not only run because of new evidence
	optimizingAUGUSTUS(train_file,species,cores_available)
#Calling function to make the AUGUSTUS prediction (gene set 2)
augPred = path_to_output + "/augustus.gff"
if(start_at <=8 and stop_at>8):
	predictWithAUGUSTUS(genome_file_name,species,ip_hints, augPred,ext_cfg_file)
#Calling function to join the two gene sets
joinedPred = path_to_output + "/GMandAugCombined.gff"
if(start_at <=9 and stop_at>9):
	combineGMandAUGUSTUS(training_genes, augPred, ip_hints, joinedPred)
#Calling the function to compare the results to the given reference annotation
if(args.compare_to_ref and start_at <=10 and stop_at >10):
	comparison = path_to_output + "/ComparisonToReferenceAnnotation.txt"
	compareToRef(joinedPred, reference, ip_hints, comparison)
### deleting all temporary results
if(delete_results):
	delete_intermediate = subprocess.run("find " + path_to_output + "/ -type f -not -name 'GMandAugCombined.gff' -delete", shell=True) #deletes all files except the final output file
	delete_intermediate_folders = subprocess.run("find " + path_to_output + "/ -maxdepth 1 -mindepth 1 -type d -exec rm -rf {} \;", shell=True) #deletes all subdirectories that formerly contained intermediate results 
